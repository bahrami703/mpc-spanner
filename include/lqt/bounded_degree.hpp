#pragma once

#include "spanner.hpp"

#include <limits>
#include <numbers>

namespace lqt {

struct BoundedDegreeParams {
  double t = 2.0;
  double split = 0.5;
  u32 kappa = 0;
};

struct BoundedDegreeSpanner {
  std::vector<Edge> edges;
  double base_stretch = 0;
  double sink_stretch = 0;
  double base_eps = 0;
  u32 kappa = 0;
  u64 base_edges = 0;
  u32 base_max_out_degree = 0;
  u32 max_sink_in_degree = 0;
  u32 max_degree = 0;
  u64 degree_bound = 0;
};

inline u32 smallest_cone_count(double sink_stretch) {
  for (u32 k = 9; k < 1u << 20; ++k) {
    const double th = 2 * std::numbers::pi / k, c = std::cos(th) - std::sin(th);
    if (c > 0 && 1 / c <= sink_stretch) return k;
  }
  throw std::invalid_argument("lqt: the sink stretch is too close to 1");
}

namespace detail {

inline void last_leaf_representatives(const Tree<2>& T, std::vector<u32>& rep) {
  std::vector<u32> last(T.node_count);
  rep.resize(T.node_count);
  for (u32 v = T.node_count; v-- > 0;) {
    const u32 b = T.child_begin[v], e = T.child_end[v];
    if (b == e) { rep[v] = last[v] = (u32)T.point[v]; continue; }
    last[v] = last[e - 1];
    rep[v] = e - b == 1 ? rep[b] : last[b];
  }
}

struct EdgeBlocks {
  static constexpr size_t block_size = 1 << 20;
  std::vector<std::vector<u64>> blocks;
  u64 count = 0;
  void push(u64 x) {
    if (blocks.empty() || blocks.back().size() == block_size) { blocks.emplace_back(); blocks.back().reserve(block_size); }
    blocks.back().push_back(x);
    ++count;
  }
};

struct BaseEdgeVisitor {
  const u32* rep;
  EdgeBlocks edges;
  std::vector<u32> in_degree;
  void operator()(u32 a, u32 b) {
    const u32 p = rep[a], q = rep[b];
    ++in_degree[q];
    edges.push((u64)p | (u64)q << 32);
  }
};

struct Cones {
  u32 kappa = 0, L = 0;
  double scale = 0;
  std::vector<double> b, ux, uy;
  std::vector<u32> lut;
  static double pseudo_angle(double x, double y) {
    const double ax = std::fabs(x), ay = std::fabs(y), r = ay / (ax + ay);
    return y >= 0 ? (x >= 0 ? r : 2 - r) : (x < 0 ? 2 + r : 4 - r);
  }
  explicit Cones(u32 k) : kappa(k) {
    const double th = 2 * std::numbers::pi / k;
    b.resize(k + 1); ux.resize(k); uy.resize(k);
    for (u32 c = 0; c < k; ++c) {
      b[c] = c ? pseudo_angle(std::cos(c * th), std::sin(c * th)) : 0.0;
      ux[c] = std::cos((c + 0.5) * th); uy[c] = std::sin((c + 0.5) * th);
    }
    b[k] = std::numeric_limits<double>::infinity();
    L = 1;
    while (L < 64 * k) L <<= 1;
    scale = L / 4.0;
    lut.resize(L);
    for (u32 i = 0, c = 0; i < L; ++i) {
      const double s = i / scale;
      while (b[c + 1] <= s) ++c;
      lut[i] = c;
    }
  }
  u32 operator()(double x, double y) const {
    const double p = pseudo_angle(x, y);
    u32 c = lut[std::min<u32>(L - 1, (u32)(p * scale))];
    c += p >= b[c + 1];
    c -= p < b[c];
    return c;
  }
};

inline void replace_in_stars_by_sink_spanners(BaseEdgeVisitor& base, const std::vector<Point<2>>& pts, u32 kappa, BoundedDegreeSpanner& S) {
  const u32 n = (u32)pts.size();
  const u64 E = base.edges.count;
  std::vector<u64> off(n + 1);
  off[0] = 0;
  u32 maxlen = 0;
  for (u32 q = 0; q < n; ++q) { off[q + 1] = off[q] + base.in_degree[q]; maxlen = std::max(maxlen, base.in_degree[q]); }
  std::vector<u32> W(E);
  std::vector<u64> pos(off.begin(), off.end() - 1);
  for (const auto& block : base.edges.blocks)
    for (u64 x : block) W[pos[x >> 32]++] = (u32)x;
  base.edges.blocks.clear();
  base.edges.blocks.shrink_to_fit();
  base.in_degree.clear();
  base.in_degree.shrink_to_fit();

  std::vector<u64> A;
  A.reserve(E);
  const Cones C(kappa);
  std::vector<u32> outd(n, 0);
  std::vector<u32> I(maxlen), Pq(maxlen), ids(maxlen), cone(maxlen), cnt(kappa, 0), start(kappa), touched(kappa + 1);
  std::vector<u64> lx(maxlen), ly(maxlen), tx(maxlen), ty(maxlen);
  std::vector<double> fx(maxlen), fy(maxlen), proj(maxlen);
  struct Frame { u32 b, e, apex; };
  std::vector<Frame> st;
  const int nbits = bit_length(n), npass = (nbits + 10) / 11, dbits = (nbits + npass - 1) / npass;
  const u32 dmask = (1u << dbits) - 1;
  for (u32 q = 0; q < n; ++q) {
    u32 Lq = (u32)(off[q + 1] - off[q]);
    std::copy(W.data() + off[q], W.data() + off[q] + Lq, I.data());
    if (Lq < 128) std::sort(I.begin(), I.begin() + Lq);
    else {
      u32 *a = I.data(), *b = ids.data();
      for (int sh = 0; sh < nbits; sh += dbits) {
        u32 c[2049];
        std::fill(c, c + dmask + 2, 0u);
        for (u32 j = 0; j < Lq; ++j) ++c[((a[j] >> sh) & dmask) + 1];
        for (u32 d = 0; d < dmask; ++d) c[d + 1] += c[d];
        for (u32 j = 0; j < Lq; ++j) b[c[(a[j] >> sh) & dmask]++] = a[j];
        std::swap(a, b);
      }
      if (a != I.data()) std::copy(a, a + Lq, I.data());
    }
    Lq = (u32)(std::unique(I.begin(), I.begin() + Lq) - I.begin());
    for (u32 j = 0; j < Lq; ++j) ++outd[I[j]];
    if (!Lq) continue;
    for (u32 j = 0; j < Lq; ++j) { const Point<2>& p = pts[I[j]]; lx[j] = p[0]; ly[j] = p[1]; }
    st.clear();
    st.push_back({0, Lq, q});
    while (!st.empty()) {
      const Frame f = st.back();
      st.pop_back();
      const u32 len = f.e - f.b, m = len + 1;
      if (len <= 2) {
        Pq[f.b] = f.apex;
        if (len == 2) Pq[f.b + 1] = f.apex;
        S.max_sink_in_degree = std::max(S.max_sink_in_degree, len);
        continue;
      }
      u32* Ib = I.data() + f.b;
      u64 *X = lx.data() + f.b, *Y = ly.data() + f.b;
      const Point<2> a = f.b ? Point<2>{X[-1], Y[-1]} : pts[q];
      for (u32 j = 0; j < len; ++j) {
        ids[j] = Ib[j]; tx[j] = X[j]; ty[j] = Y[j];
        const double dx = (double)(i64)(X[j] - a[0]), dy = (double)(i64)(Y[j] - a[1]);
        fx[j] = dx; fy[j] = dy;
        cone[j] = C(dx, dy);
      }
      u32 nt = 0;
      for (u32 j = 0; j < len; ++j) { const u32 c = cone[j]; touched[nt] = c; nt += cnt[c] == 0; ++cnt[c]; }
      for (u32 i = 0, o = 0; i < nt; ++i) { const u32 c = touched[i]; start[c] = o; o += cnt[c]; }
      for (u32 j = 0; j < len; ++j) {
        const u32 c = cone[j], k = start[c]++;
        Ib[k] = ids[j]; X[k] = tx[j]; Y[k] = ty[j];
        proj[k] = fx[j] * C.ux[c] + fy[j] * C.uy[c];
      }
      u32 subsets = 0;
      auto emit = [&](u32 sb, u32 se) {
        u32 rk = sb;
        double bp = proj[sb];
        for (u32 k = sb + 1; k < se; ++k) { const bool lt = proj[k] < bp; bp = lt ? proj[k] : bp; rk = lt ? k : rk; }
        std::swap(Ib[sb], Ib[rk]); std::swap(X[sb], X[rk]); std::swap(Y[sb], Y[rk]);
        Pq[f.b + sb] = f.apex;
        ++subsets;
        if (se > sb + 1) st.push_back({f.b + sb + 1, f.b + se, Ib[sb]});
      };
      for (u32 i = 0; i < nt; ++i) {
        const u32 c = touched[i], s = cnt[c], ge = start[c], g = ge - s;
        cnt[c] = 0;
        if (2 * s <= m) emit(g, ge);
        else { emit(g, g + s / 2); emit(g + s / 2, ge); }
      }
      S.max_sink_in_degree = std::max(S.max_sink_in_degree, subsets);
    }
    for (u32 j = 0; j < Lq; ++j) A.push_back((u64)I[j] | (u64)Pq[j] << 32);
  }
  const u64 w = A.size();
  for (u32 x = 0; x < n; ++x) S.base_max_out_degree = std::max(S.base_max_out_degree, outd[x]);

  std::vector<u64> h(n + 1, 0);
  for (u64 j = 0; j < w; ++j) {
    const u32 x = (u32)A[j], y = (u32)(A[j] >> 32);
    ++h[std::max(x, y) + 1];
  }
  for (u32 x = 0; x < n; ++x) h[x + 1] += h[x];
  pos.assign(h.begin(), h.end() - 1);
  for (u64 j = 0; j < w; ++j) {
    const u32 x = (u32)A[j], y = (u32)(A[j] >> 32);
    W[pos[std::max(x, y)]++] = std::min(x, y);
  }
  A.clear();
  A.shrink_to_fit();
  std::vector<u32> deg(n, 0);
  std::vector<u32> mark(n, ~0u);
  u64 k = 0;
  for (u32 x = 0; x < n; ++x) {
    const u64 b = h[x], e = h[x + 1];
    h[x] = k;
    for (u64 j = b; j < e; ++j) {
      const u32 y = W[j];
      if (mark[y] == x) continue;
      mark[y] = x; W[k++] = y; ++deg[x]; ++deg[y];
    }
  }
  h[n] = k;
  for (u32 x = 0; x < n; ++x) S.max_degree = std::max(S.max_degree, deg[x]);
  S.degree_bound = S.base_max_out_degree + (u64)(S.base_max_out_degree + 1) * (kappa + 1);
  S.edges.resize(k);
  for (u32 x = 0; x < n; ++x)
    for (u64 j = h[x]; j < h[x + 1]; ++j) S.edges[j] = {W[j], x};
}

}

inline BoundedDegreeSpanner build_bounded_degree_spanner(const PointSet<2>& points, const Tree<2>& tree, const BoundedDegreeParams& params = {}) {
  if (!(params.t > 1) || !std::isfinite(params.t)) throw std::invalid_argument("lqt: the stretch t must be finite and greater than 1");
  if (!(params.split > 0 && params.split < 1)) throw std::invalid_argument("lqt: split must be in (0, 1)");
  if (tree.point_count != points.size()) throw std::invalid_argument("lqt: the tree was not built on these points");
  BoundedDegreeSpanner S;
  const double t = params.t, tp = std::pow(t, params.split);
  if (!(tp > 1)) throw std::invalid_argument("lqt: the base stretch t^split must be greater than 1");
  const u32 kappa = params.kappa ? params.kappa : smallest_cone_count(t / tp);
  const double th = 2 * std::numbers::pi / kappa, tsink = 1 / (std::cos(th) - std::sin(th));
  if (kappa < 9 || tsink <= 0 || tp * tsink > t * (1 + 1e-12))
    throw std::invalid_argument("lqt: kappa too small for the stretch t");
  S.base_stretch = tp;
  S.sink_stretch = tsink;
  S.base_eps = (tp - 1) / (2 * (tp + 1));
  S.kappa = kappa;
  std::vector<u32> rep;
  detail::last_leaf_representatives(tree, rep);
  detail::BaseEdgeVisitor base{rep.data(), {}, std::vector<u32>(points.size(), 0)};
  {
    auto nodes = detail::make_wspd_nodes<2>(tree, SeparationTest::BoundingBox);
    S.base_edges = detail::run_wspd<2>(nodes, S.base_eps, base);
  }
  rep.clear();
  rep.shrink_to_fit();
  detail::replace_in_stars_by_sink_spanners(base, points.points, kappa, S);
  return S;
}

}

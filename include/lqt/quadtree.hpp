#pragma once

#include "point_set.hpp"

#include <span>

namespace lqt {

enum class NodeKind : u8 { Grid = 0, EmptyRegion = 1, LeafCellRegion = 2, Donut = 3, Compressed = 4, Point = 5 };

struct LevelQuadtreeParams {
  double beta = 0;
  double lambda = 0.5;
};

struct LevelQuadtreeInfo {
  int levels = 0;
  u64 calls = 0;
  u64 base_calls = 0;
  u64 forced_base_calls = 0;
  u64 donuts = 0;
};

template <int D> struct Tree {
  u32 node_count = 0;
  u32 point_count = 0;
  int bits = 0;
  u32 height = 0;
  std::vector<u32> parent, child_begin, child_end, size, min_input_index, leaf_begin, depth;
  std::vector<u32> light_representative, first_representative;
  std::vector<u16> call_depth;
  std::vector<u8> cell_level, call_level, hole_level;
  std::vector<NodeKind> kind;
  std::vector<i32> point;
  std::vector<Point<D>> cell_corner, box_lo, box_hi, hole_corner;
  std::vector<u32> leaf_points;

  bool is_leaf(u32 v) const { return child_begin[v] == child_end[v]; }
  std::span<const u32> points_of(u32 v) const { return {leaf_points.data() + leaf_begin[v], size[v]}; }
};

namespace detail {

template <int D> struct Arena {
  std::vector<i32> parent;
  std::vector<u8> level, kind, call_level;
  std::vector<Point<D>> corner;
  std::vector<i32> point;
  std::vector<u8> hole_level;
  std::vector<Point<D>> hole_corner;
  u32 size() const { return (u32)parent.size(); }
  void reserve(size_t m) {
    parent.reserve(m); level.reserve(m); kind.reserve(m); call_level.reserve(m);
    corner.reserve(m); point.reserve(m); hole_level.reserve(m); hole_corner.reserve(m);
  }
  u32 add(i32 p, int cell_level, NodeKind k, int lev, const Point<D>& c, i32 pt = -1) {
    parent.push_back(p); level.push_back((u8)cell_level); kind.push_back((u8)k); call_level.push_back((u8)lev);
    corner.push_back(c); point.push_back(pt); hole_level.push_back(0xff); hole_corner.push_back(Point<D>{});
    return (u32)parent.size() - 1;
  }
};

template <int D> struct Builder {
  const PointSet<D>& P;
  Arena<D> A;
  std::vector<u32> stack, leaves, scratch;

  explicit Builder(const PointSet<D>& p) : P(p) { A.reserve((size_t)P.size() * 3 + 16); }

  u32 compressed(const u32* Q, u32 m, const Point<D>& C, int e, int lev, i32 parent) {
    if (m == 1) return A.add(parent, e, NodeKind::Point, lev, C, (i32)Q[0]);
    stack.clear();
    leaves.clear();
    u32 cur = A.add(-1, 0, NodeKind::Point, lev, Point<D>{}, (i32)Q[0]);
    leaves.push_back(cur);
    for (u32 i = 1; i < m; ++i) {
      int h = common_cell_level<D>(P.keys[Q[i - 1]], P.keys[Q[i]]);
      u32 x = cur;
      while (!stack.empty() && A.level[stack.back()] < h) { A.parent[x] = (i32)stack.back(); x = stack.back(); stack.pop_back(); }
      if (!stack.empty() && A.level[stack.back()] == h) {
        A.parent[x] = (i32)stack.back();
      } else {
        u32 node = A.add(-1, h, NodeKind::Compressed, lev, cell_corner_of<D>(P.points[Q[i]], h));
        A.parent[x] = (i32)node;
        stack.push_back(node);
      }
      cur = A.add(-1, 0, NodeKind::Point, lev, Point<D>{}, (i32)Q[i]);
      leaves.push_back(cur);
    }
    u32 x = cur;
    while (!stack.empty()) { A.parent[x] = (i32)stack.back(); x = stack.back(); stack.pop_back(); }
    for (u32 l : leaves) {
      int le = A.level[A.parent[l]] - 1;
      A.level[l] = (u8)le;
      A.corner[l] = cell_corner_of<D>(P.points[A.point[l]], le);
    }
    if (A.level[x] < e) {
      u32 w = A.add(parent, e, NodeKind::Compressed, lev, C);
      A.parent[x] = (i32)w;
      return w;
    }
    A.parent[x] = parent;
    return x;
  }

  u32 build_compressed() {
    std::vector<u32> all(P.size());
    for (u32 i = 0; i < P.size(); ++i) all[i] = i;
    return compressed(all.data(), P.size(), Point<D>{}, P.bits, 0, -1);
  }

  struct Call {
    u32 offset, m, parent_m;
    Point<D> C;
    u8 e;
    i32 attach;
  };
  struct Region {
    u32 node;
    Point<D> corner;
    u8 e;
    u32 a1, b1, a2, b2;
    u32 count;
  };
  struct CellItem {
    Point<D> corner;
    u8 e;
    u32 zlo[D], zhi[D];
    u32 qa, qb;
    i32 parent;
  };

  std::vector<u32> morton_order[2];
  std::vector<u32> axis_order[2][D];
  std::vector<i32> region_of;
  std::vector<u64> grid[D];
  std::vector<Region> regions;
  std::vector<CellItem> items;
  std::vector<u32> slot_offset, cursor;

  u32 lower_bound(const u32* Q, u32 a, u32 b, u128 k) const {
    while (a < b) {
      u32 mid = a + (b - a) / 2;
      if (P.keys[Q[mid]] < k) a = mid + 1; else b = mid;
    }
    return a;
  }

  u32 walk_grid_tree(const Call& c, int lev, const u32* Q, LevelQuadtreeInfo& info) {
    items.clear();
    CellItem first;
    first.corner = c.C; first.e = c.e; first.qa = c.offset; first.qb = c.offset + c.m; first.parent = c.attach;
    for (int i = 0; i < D; ++i) { first.zlo[i] = 0; first.zhi[i] = (u32)grid[i].size(); }
    items.push_back(first);
    i32 call_root = -1;
    bool at_root = true;
    while (!items.empty()) {
      CellItem it = items.back();
      items.pop_back();
      bool root = at_root;
      at_root = false;
      if (it.qa == it.qb) continue;
      u64 g = 1;
      for (int i = 0; i < D; ++i) g *= (u64)(it.zhi[i] - it.zlo[i]);
      if (g <= 1) {
        u32 v = A.add(it.parent, it.e, g == 0 ? NodeKind::EmptyRegion : NodeKind::LeafCellRegion, lev, it.corner);
        regions.push_back({v, it.corner, it.e, it.qa, it.qb, it.qb, it.qb, it.qb - it.qa});
        if (root) call_root = (i32)v;
        continue;
      }
      int h = 0;
      for (int i = 0; i < D; ++i) h = std::max(h, bit_length(grid[i][it.zlo[i]] ^ grid[i][it.zhi[i] - 1]));
      Point<D> cc;
      for (int i = 0; i < D; ++i) cc[i] = (grid[i][it.zlo[i]] >> h) << h;
      const u128 cstart = morton_key<D>(cc);
      i32 par = it.parent;
      u32 qa2 = it.qa, qb2 = it.qb;
      if (h < it.e) {
        if (root) {
          u32 w = A.add(par, it.e, NodeKind::Grid, lev, it.corner);
          call_root = (i32)w;
          par = (i32)w;
        }
        qa2 = lower_bound(Q, it.qa, it.qb, cstart);
        qb2 = lower_bound(Q, qa2, it.qb, cstart + ((u128)1 << (D * h)));
        u32 donut_count = (qa2 - it.qa) + (it.qb - qb2);
        if (donut_count) {
          u32 dn = A.add(par, it.e, NodeKind::Donut, lev, it.corner);
          A.hole_level[dn] = (u8)h;
          A.hole_corner[dn] = cc;
          ++info.donuts;
          regions.push_back({dn, it.corner, it.e, it.qa, qa2, qb2, it.qb, donut_count});
        }
      }
      if (qa2 == qb2) continue;
      u32 B = A.add(par, h, NodeKind::Grid, lev, cc);
      if (root && h == it.e) call_root = (i32)B;
      constexpr int child_cells = 1 << D;
      u32 qs[child_cells + 1];
      qs[0] = qa2;
      qs[child_cells] = qb2;
      const u128 step = (u128)1 << (D * (h - 1));
      for (int t = 1; t < child_cells; ++t) qs[t] = lower_bound(Q, qs[t - 1], qb2, cstart + step * (u128)t);
      u32 zs[D];
      for (int i = 0; i < D; ++i) {
        const u64* z = grid[i].data();
        u32 a = it.zlo[i], b = it.zhi[i];
        while (a < b) {
          u32 mid = a + (b - a) / 2;
          if (((z[mid] >> (h - 1)) & 1) == 0) a = mid + 1; else b = mid;
        }
        zs[i] = a;
      }
      for (int t = child_cells - 1; t >= 0; --t) {
        if (qs[t] == qs[t + 1]) continue;
        CellItem ch;
        ch.e = (u8)(h - 1);
        ch.qa = qs[t];
        ch.qb = qs[t + 1];
        ch.parent = (i32)B;
        for (int i = 0; i < D; ++i) {
          bool up = (t >> i) & 1;
          ch.corner[i] = cc[i] + ((u64)up << (h - 1));
          ch.zlo[i] = up ? zs[i] : it.zlo[i];
          ch.zhi[i] = up ? it.zhi[i] : zs[i];
        }
        items.push_back(ch);
      }
    }
    return (u32)call_root;
  }

  u32 build_level(double beta, double lambda, LevelQuadtreeInfo& info) {
    const u32 n = P.size();
    const double threshold = std::pow((double)n, lambda);
    for (int b = 0; b < 2; ++b) {
      morton_order[b].resize(n);
      for (int i = 0; i < D; ++i) axis_order[b][i].resize(n);
    }
    for (u32 i = 0; i < n; ++i) morton_order[0][i] = i;
    {
      struct CoordinateRecord { u64 c; u32 id; };
      std::vector<CoordinateRecord> r(n);
      for (int i = 0; i < D; ++i) {
        for (u32 j = 0; j < n; ++j) r[j] = {P.points[j][i], j};
        if (n < 4096)
          std::stable_sort(r.begin(), r.end(), [](const CoordinateRecord& a, const CoordinateRecord& b) { return a.c < b.c; });
        else
          radix_sort(r, P.bits, [](const CoordinateRecord& x, int p) -> u32 { return (u32)((x.c >> (16 * p)) & 0xffff); });
        for (u32 j = 0; j < n; ++j) axis_order[0][i][j] = r[j].id;
      }
    }
    region_of.assign(n, -1);
    std::vector<Call> current, next;
    current.push_back({0, n, 0, Point<D>{}, (u8)P.bits, -1});
    i32 root = -1;
    int cb = 0;
    for (int lev = 0; !current.empty(); ++lev) {
      next.clear();
      u32 next_offset = 0;
      u32* Qm = morton_order[cb].data();
      for (const Call& c : current) {
        info.levels = std::max(info.levels, lev + 1);
        bool base = (double)c.m <= threshold || (c.parent_m && c.m >= c.parent_m);
        if (base) {
          if ((double)c.m > threshold) ++info.forced_base_calls;
          ++info.base_calls;
          u32 r = compressed(Qm + c.offset, c.m, c.C, c.e, lev, c.attach);
          if (c.attach < 0) root = (i32)r;
          continue;
        }
        ++info.calls;
        const u64 k = ceil_power(c.m, 1 - beta), s = c.m / k;
        for (int i = 0; i < D; ++i) {
          grid[i].clear();
          const u32* a = axis_order[cb][i].data() + c.offset;
          for (u64 j = 1; j <= s; ++j) {
            u64 v = P.points[a[j * k - 1]][i];
            if (grid[i].empty() || grid[i].back() != v) grid[i].push_back(v);
          }
        }
        regions.clear();
        u32 r = walk_grid_tree(c, lev, Qm, info);
        if (c.attach < 0) root = (i32)r;
        slot_offset.clear();
        bool any_recursive = false;
        for (const Region& R : regions) {
          bool region_base = (double)R.count <= threshold || R.count >= c.m;
          i32 tag;
          if (region_base) {
            if ((double)R.count > threshold) ++info.forced_base_calls;
            info.levels = std::max(info.levels, lev + 2);
            ++info.base_calls;
            scratch.assign(Qm + R.a1, Qm + R.b1);
            scratch.insert(scratch.end(), Qm + R.a2, Qm + R.b2);
            compressed(scratch.data(), R.count, R.corner, R.e, lev + 1, (i32)R.node);
            tag = -1;
          } else {
            u32 o = next_offset;
            next_offset += R.count;
            u32* dst = morton_order[cb ^ 1].data() + o;
            std::memcpy(dst, Qm + R.a1, (R.b1 - R.a1) * sizeof(u32));
            std::memcpy(dst + (R.b1 - R.a1), Qm + R.a2, (R.b2 - R.a2) * sizeof(u32));
            tag = (i32)slot_offset.size();
            slot_offset.push_back(o);
            next.push_back({o, R.count, c.m, R.corner, R.e, (i32)R.node});
            any_recursive = true;
          }
          for (u32 t = R.a1; t < R.b1; ++t) region_of[Qm[t]] = tag;
          for (u32 t = R.a2; t < R.b2; ++t) region_of[Qm[t]] = tag;
        }
        if (any_recursive) {
          for (int i = 0; i < D; ++i) {
            cursor = slot_offset;
            const u32* a = axis_order[cb][i].data() + c.offset;
            u32* out = axis_order[cb ^ 1][i].data();
            for (u32 t = 0; t < c.m; ++t) {
              u32 p = a[t];
              i32 rr = region_of[p];
              if (rr >= 0) out[cursor[rr]++] = p;
            }
          }
        }
      }
      std::swap(current, next);
      cb ^= 1;
    }
    return (u32)root;
  }
};

template <int D> Tree<D> prune_and_number(const Arena<D>& A, u32 root, const PointSet<D>& P) {
  const u32 M = A.size();
  std::vector<u32> coff(M + 1, 0), children(M);
  for (u32 v = 0; v < M; ++v) if (A.parent[v] >= 0) ++coff[A.parent[v] + 1];
  for (u32 v = 0; v < M; ++v) coff[v + 1] += coff[v];
  {
    std::vector<u32> pos(coff.begin(), coff.end() - 1);
    for (u32 v = 0; v < M; ++v) if (A.parent[v] >= 0) children[pos[A.parent[v]]++] = v;
  }
  std::vector<u32> order;
  order.reserve(M);
  order.push_back(root);
  for (size_t h = 0; h < order.size(); ++h) {
    u32 v = order[h];
    for (u32 j = coff[v]; j < coff[v + 1]; ++j) order.push_back(children[j]);
  }
  std::vector<u32> count(M, 0);
  for (size_t h = order.size(); h-- > 0;) {
    u32 v = order[h];
    if (A.point[v] >= 0) count[v] = 1;
    if (A.parent[v] >= 0 && v != root) count[A.parent[v]] += count[v];
  }
  Tree<D> T;
  T.point_count = P.size();
  T.bits = P.bits;
  u32 N = 0;
  for (u32 v : order) N += count[v] > 0;
  order.assign(N, 0);
  std::vector<u32> cb(N), ce(N);
  {
    std::vector<u32> new_id(M, ~0u), st;
    order[0] = root;
    new_id[root] = 0;
    u32 next = 1;
    st.push_back(root);
    while (!st.empty()) {
      u32 v = st.back();
      st.pop_back();
      u32 first = next;
      for (u32 j = coff[v]; j < coff[v + 1]; ++j) {
        u32 w = children[j];
        if (count[w]) { new_id[w] = next; order[next++] = w; }
      }
      cb[new_id[v]] = first;
      ce[new_id[v]] = next;
      for (u32 k = next; k-- > first;) st.push_back(order[k]);
    }
  }
  T.node_count = N;
  T.child_begin = std::move(cb);
  T.child_end = std::move(ce);
  std::vector<u32> new_id(M, ~0u);
  for (u32 v = 0; v < N; ++v) new_id[order[v]] = v;
  T.parent.resize(N); T.size.resize(N); T.min_input_index.resize(N); T.leaf_begin.resize(N); T.depth.resize(N);
  T.light_representative.resize(N); T.first_representative.resize(N); T.call_depth.resize(N);
  T.cell_level.resize(N); T.kind.resize(N); T.call_level.resize(N); T.hole_level.resize(N); T.point.resize(N);
  T.cell_corner.resize(N); T.box_lo.resize(N); T.box_hi.resize(N); T.hole_corner.resize(N);
  for (u32 v = 0; v < N; ++v) {
    u32 o = order[v];
    T.parent[v] = v == 0 ? ~0u : new_id[A.parent[o]];
    T.cell_level[v] = A.level[o]; T.kind[v] = (NodeKind)A.kind[o]; T.call_level[v] = A.call_level[o]; T.point[v] = A.point[o];
    T.cell_corner[v] = A.corner[o]; T.hole_level[v] = A.hole_level[o]; T.hole_corner[v] = A.hole_corner[o];
  }
  for (u32 v = N; v-- > 0;) {
    if (T.child_begin[v] == T.child_end[v]) {
      u32 p = (u32)T.point[v];
      T.size[v] = 1; T.box_lo[v] = T.box_hi[v] = P.points[p]; T.min_input_index[v] = P.input_index[p];
      T.light_representative[v] = T.first_representative[v] = p;
      continue;
    }
    u32 b = T.child_begin[v], en = T.child_end[v];
    u32 c = 0, min_index = ~0u, lightest = b;
    Point<D> l = T.box_lo[b], h = T.box_hi[b];
    for (u32 w = b; w < en; ++w) {
      c += T.size[w];
      min_index = std::min(min_index, T.min_input_index[w]);
      for (int i = 0; i < D; ++i) { l[i] = std::min(l[i], T.box_lo[w][i]); h[i] = std::max(h[i], T.box_hi[w][i]); }
      if (T.size[w] < T.size[lightest] || (T.size[w] == T.size[lightest] && T.min_input_index[w] < T.min_input_index[lightest])) lightest = w;
    }
    T.size[v] = c; T.min_input_index[v] = min_index; T.box_lo[v] = l; T.box_hi[v] = h;
    T.light_representative[v] = T.light_representative[lightest];
    T.first_representative[v] = T.first_representative[b];
  }
  T.leaf_points.resize(P.size());
  T.depth[0] = 0; T.call_depth[0] = 0; T.leaf_begin[0] = 0;
  u32 height = 0;
  for (u32 v = 0; v < N; ++v) {
    u32 off = T.leaf_begin[v];
    for (u32 w = T.child_begin[v]; w < T.child_end[v]; ++w) {
      T.depth[w] = T.depth[v] + 1;
      T.call_depth[w] = T.call_level[w] == T.call_level[v] ? (u16)std::min<u32>(T.call_depth[v] + 1, 65535) : 0;
      T.leaf_begin[w] = off;
      off += T.size[w];
    }
    if (T.child_begin[v] == T.child_end[v]) T.leaf_points[T.leaf_begin[v]] = (u32)T.point[v];
    height = std::max(height, T.depth[v]);
  }
  T.height = height;
  return T;
}

template <int D> void require_points(const PointSet<D>& P) {
  if (P.size() == 0) throw std::invalid_argument("lqt: the point set is empty");
}

}

template <int D>
Tree<D> build_level_quadtree(const PointSet<D>& points, LevelQuadtreeParams params = {}, LevelQuadtreeInfo* info = nullptr) {
  detail::require_points(points);
  const double beta = params.beta > 0 ? params.beta : 1.0 / D;
  if (!(beta > 0 && beta <= 1.0 / D + 1e-12)) throw std::invalid_argument("lqt: beta must be in (0, 1/d]");
  if (!(params.lambda > 0 && params.lambda <= 1)) throw std::invalid_argument("lqt: lambda must be in (0, 1]");
  LevelQuadtreeInfo local;
  u32 root;
  detail::Arena<D> arena;
  {
    detail::Builder<D> builder(points);
    root = builder.build_level(beta, params.lambda, info ? *info : local);
    arena = std::move(builder.A);
  }
  return detail::prune_and_number<D>(arena, root, points);
}

template <int D> Tree<D> build_compressed_quadtree(const PointSet<D>& points) {
  detail::require_points(points);
  detail::Builder<D> builder(points);
  u32 root = builder.build_compressed();
  return detail::prune_and_number<D>(builder.A, root, points);
}

}

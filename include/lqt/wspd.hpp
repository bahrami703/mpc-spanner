#pragma once

#include "quadtree.hpp"

#include <type_traits>

namespace lqt {

enum class SeparationTest { BoundingBox, Cell };

struct CountOnly {
  static constexpr bool counts_only = true;
  void operator()(u32, u32) const {}
};

namespace detail {

template <int D> struct alignas(D == 2 ? 64 : 16) WspdNode {
  Point<D> lo, hi;
  double diameter;
  u64 tie_key;
  u32 child_begin, child_end;
  u32 light_representative, first_representative;
};

template <int D> std::vector<WspdNode<D>> make_wspd_nodes(const Tree<D>& T, SeparationTest test) {
  const u32 N = T.node_count;
  std::vector<WspdNode<D>> W(N);
  const double sqrt_d = std::sqrt((double)D);
  for (u32 v = 0; v < N; ++v) {
    WspdNode<D>& w = W[v];
    w.child_begin = T.child_begin[v]; w.child_end = T.child_end[v];
    w.light_representative = T.light_representative[v]; w.first_representative = T.first_representative[v];
    if (test == SeparationTest::BoundingBox || T.size[v] == 1) {
      w.lo = T.box_lo[v]; w.hi = T.box_hi[v];
    } else {
      w.lo = T.cell_corner[v];
      for (int i = 0; i < D; ++i) w.hi[i] = T.cell_corner[v][i] + (1ull << T.cell_level[v]);
    }
    if (T.size[v] == 1) w.diameter = 0;
    else if (test == SeparationTest::Cell) w.diameter = sqrt_d * std::ldexp(1.0, T.cell_level[v]);
    else {
      u128 s = 0;
      for (int i = 0; i < D; ++i) { u64 g = T.box_hi[v][i] - T.box_lo[v][i]; s += (u128)g * g; }
      w.diameter = std::sqrt(to_double(s));
    }
    w.tie_key = ~(((u64)T.call_level[v] << 56) | ((u64)std::min<u32>(T.call_depth[v], 0xffffff) << 32) | T.min_input_index[v]);
  }
  return W;
}

template <class Visitor> constexpr bool counts_only_v = requires { Visitor::counts_only; };

template <int D, class Visitor> struct WspdLoop {
  const WspdNode<D>* w;
  double eps, eps2;
  Visitor& visit;
  std::vector<u32>* split_counts = nullptr;
  u64 pairs = 0, examined = 0, exact_tests = 0;
  std::vector<u64> buffer;

  WspdLoop(const WspdNode<D>* nodes, double e, Visitor& v) : w(nodes), eps(e), eps2(e * e), visit(v) {
    buffer.resize(1 << 16);
  }

  inline bool examine(u32 a, u32 b, u64& packed, u64& exact) const {
    const WspdNode<D>& A = w[a];
    const WspdNode<D>& B = w[b];
    const bool a_larger = A.diameter > B.diameter || (A.diameter == B.diameter && A.tie_key > B.tie_key);
    const double dx = a_larger ? A.diameter : B.diameter;
    u64 g[D];
    double d2f = 0;
    for (int i = 0; i < D; ++i) {
      u64 g1 = A.lo[i] > B.hi[i] ? A.lo[i] - B.hi[i] : 0;
      u64 g2 = B.lo[i] > A.hi[i] ? B.lo[i] - A.hi[i] : 0;
      g[i] = g1 | g2;
      double gf = (double)g[i];
      d2f += gf * gf;
    }
    const double lhs = dx * dx, rhs = eps2 * d2f;
    bool separated = lhs < rhs * (1 - 1e-9);
    if (__builtin_expect(!separated && !(lhs > rhs * (1 + 1e-9)), 0)) {
      u128 d2 = 0;
      for (int i = 0; i < D; ++i) d2 += (u128)g[i] * g[i];
      separated = dx <= eps * std::sqrt(to_double(d2));
      ++exact;
    }
    const u32 x = a_larger ? a : b, y = a_larger ? b : a;
    if (!separated) {
      const WspdNode<D>* cx = w + (a_larger ? A.child_begin : B.child_begin);
      __builtin_prefetch(cx);
      __builtin_prefetch(cx + 2);
    }
    packed = (u64)x | ((u64)y << 32);
    return separated;
  }

  void run_from(u32 z) {
    const u32 b0 = w[z].child_begin, e0 = w[z].child_end;
    u64 need = (u64)(e0 - b0) * (e0 - b0) + 64;
    if (buffer.size() < need) buffer.resize(need);
    u64* st = buffer.data();
    u64 cap = buffer.size(), top = 0, np = 0, nc = 0, ex = 0;
    for (u32 i = b0; i < e0; ++i)
      for (u32 j = i + 1; j < e0; ++j) {
        u64 pk;
        bool sp = examine(i, j, pk, ex);
        ++nc; np += sp;
        if constexpr (!counts_only_v<Visitor>) if (sp) visit(i, j);
        st[top] = pk; top += !sp;
      }
    constexpr u32 max_children = 2u << D;
    while (top) {
      const u64 p = st[--top];
      const u32 x = (u32)p, y = (u32)(p >> 32);
      if (split_counts) ++(*split_counts)[x];
      const u32 cb = w[x].child_begin, ce = w[x].child_end;
      if (__builtin_expect(top + max_children > cap, 0)) { buffer.resize(2 * cap); st = buffer.data(); cap = buffer.size(); }
      u64 pk[max_children];
      bool sp[max_children];
      u32 k = 0;
      for (u32 c = cb; c < ce; ++c, ++k) sp[k] = examine(c, y, pk[k], ex);
      nc += k;
      for (u32 j = 0; j < k; ++j) {
        np += sp[j];
        if constexpr (!counts_only_v<Visitor>) if (sp[j]) visit(cb + j, y);
        st[top] = pk[j];
        top += !sp[j];
      }
    }
    pairs += np; examined += nc; exact_tests += ex;
  }
};

template <int D, class Visitor>
u64 run_wspd(const std::vector<WspdNode<D>>& W, double eps, Visitor& visit, std::vector<u32>* split_counts = nullptr) {
  WspdLoop<D, Visitor> loop(W.data(), eps, visit);
  if (split_counts) split_counts->assign(W.size(), 0);
  loop.split_counts = split_counts;
  const u32 N = (u32)W.size();
  for (u32 z = 0; z < N; ++z)
    if (W[z].child_end - W[z].child_begin >= 2) loop.run_from(z);
  return loop.pairs;
}

inline void require_separation(double eps) {
  if (!(eps > 0) || !std::isfinite(eps)) throw std::invalid_argument("lqt: eps must be positive and finite");
}

}

template <int D, class Visitor>
u64 for_each_wspd_pair(const Tree<D>& tree, double eps, Visitor&& visit, SeparationTest test = SeparationTest::BoundingBox,
                       std::vector<u32>* split_counts = nullptr) {
  detail::require_separation(eps);
  auto nodes = detail::make_wspd_nodes<D>(tree, test);
  using V = std::remove_reference_t<Visitor>;
  return detail::run_wspd<D, V>(nodes, eps, visit, split_counts);
}

template <int D> u64 count_wspd_pairs(const Tree<D>& tree, double eps, SeparationTest test = SeparationTest::BoundingBox) {
  CountOnly counter;
  return for_each_wspd_pair<D>(tree, eps, counter, test);
}

}

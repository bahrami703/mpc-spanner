#pragma once

#include "wspd.hpp"

namespace lqt {

struct Edge {
  u32 u, v;
};

enum class Representative { LightChild, FirstLeaf };

enum class SeparationRule { CallahanKosaraju, HarPeled };

struct SpannerParams {
  double t = 2.0;
  SeparationRule rule = SeparationRule::CallahanKosaraju;
  Representative representative = Representative::LightChild;
  double eps = 0;
};

inline double separation_for_stretch(double t, SeparationRule rule) {
  if (!(t > 1) || !std::isfinite(t)) throw std::invalid_argument("lqt: the stretch t must be finite and greater than 1");
  return rule == SeparationRule::HarPeled ? (t - 1) / 16.0 : (t - 1) / (2 * (t + 1));
}

inline double spanner_separation(const SpannerParams& params) {
  if (params.eps > 0) return params.eps;
  return separation_for_stretch(params.t, params.rule);
}

namespace detail {

template <int D, bool Light, class Visitor> struct RepresentativeEdgeVisitor {
  const WspdNode<D>* nodes;
  Visitor& visit;
  void operator()(u32 a, u32 b) {
    if constexpr (Light) visit(nodes[a].light_representative, nodes[b].light_representative);
    else visit(nodes[a].first_representative, nodes[b].first_representative);
  }
};

}

template <int D, class Visitor> u64 for_each_spanner_edge(const Tree<D>& tree, const SpannerParams& params, Visitor&& visit) {
  const double eps = spanner_separation(params);
  detail::require_separation(eps);
  auto nodes = detail::make_wspd_nodes<D>(tree, SeparationTest::BoundingBox);
  using V = std::remove_reference_t<Visitor>;
  if (params.representative == Representative::LightChild) {
    detail::RepresentativeEdgeVisitor<D, true, V> edge{nodes.data(), visit};
    return detail::run_wspd<D>(nodes, eps, edge);
  }
  detail::RepresentativeEdgeVisitor<D, false, V> edge{nodes.data(), visit};
  return detail::run_wspd<D>(nodes, eps, edge);
}

template <int D> std::vector<Edge> build_spanner(const Tree<D>& tree, const SpannerParams& params = {}) {
  std::vector<Edge> edges;
  for_each_spanner_edge<D>(tree, params, [&](u32 p, u32 q) { edges.push_back({p, q}); });
  return edges;
}

inline u32 max_degree(u32 point_count, const std::vector<Edge>& edges) {
  std::vector<u32> degree(point_count, 0);
  for (const Edge& e : edges) { ++degree[e.u]; ++degree[e.v]; }
  u32 m = 0;
  for (u32 d : degree) m = std::max(m, d);
  return m;
}

template <int D> std::vector<Edge> to_input_indices(const PointSet<D>& points, std::vector<Edge> edges) {
  for (Edge& e : edges) { e.u = points.input_index[e.u]; e.v = points.input_index[e.v]; }
  return edges;
}

}

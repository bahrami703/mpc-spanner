#include "lqt/lqt.hpp"

#include <cstdio>
#include <functional>
#include <limits>
#include <queue>
#include <string>

namespace {

using namespace lqt;

int failures = 0;
int checks = 0;

void expect(bool ok, const std::string& what) {
  ++checks;
  if (!ok) {
    ++failures;
    std::fprintf(stderr, "FAIL: %s\n", what.c_str());
  }
}

struct SplitMix {
  u64 state;
  u64 next() {
    u64 z = (state += 0x9e3779b97f4a7c15ull);
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ull;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebull;
    return z ^ (z >> 31);
  }
  double uniform() { return (double)(next() >> 11) * 0x1.0p-53; }
};

template <int D> std::vector<std::array<double, D>> sample(const std::string& kind, u32 n, u64 seed) {
  SplitMix g{seed};
  std::vector<std::array<double, D>> x(n);
  if (kind == "uniform") {
    for (auto& p : x)
      for (int i = 0; i < D; ++i) p[i] = g.uniform();
  } else if (kind == "clusters") {
    u32 clusters = 1;
    while (clusters * clusters < n) ++clusters;
    std::vector<std::array<double, D>> centers(clusters);
    for (auto& c : centers)
      for (int i = 0; i < D; ++i) c[i] = g.uniform() * clusters;
    for (u32 j = 0; j < n; ++j)
      for (int i = 0; i < D; ++i) x[j][i] = centers[j % clusters][i] + g.uniform();
  } else if (kind == "lattice") {
    for (auto& p : x)
      for (int i = 0; i < D; ++i) p[i] = (double)(g.next() % 7);
  } else {
    for (u32 j = 0; j < n; ++j) {
      x[j].fill(0);
      x[j][0] = std::ldexp(1.0, (int)(j % 60)) + j;
    }
  }
  return x;
}

template <int D> double distance(const PointSet<D>& P, u32 a, u32 b) {
  u128 s = 0;
  for (int i = 0; i < D; ++i) {
    u64 g = P.points[a][i] > P.points[b][i] ? P.points[a][i] - P.points[b][i] : P.points[b][i] - P.points[a][i];
    s += (u128)g * g;
  }
  return std::sqrt(detail::to_double(s));
}

template <int D> double max_stretch(const PointSet<D>& P, const std::vector<Edge>& edges) {
  const u32 n = P.size();
  std::vector<std::vector<u32>> adj(n);
  for (const Edge& e : edges) { adj[e.u].push_back(e.v); adj[e.v].push_back(e.u); }
  double worst = 0;
  std::vector<double> dist(n);
  using Item = std::pair<double, u32>;
  for (u32 s = 0; s < n; ++s) {
    std::fill(dist.begin(), dist.end(), std::numeric_limits<double>::infinity());
    std::priority_queue<Item, std::vector<Item>, std::greater<>> heap;
    dist[s] = 0;
    heap.push({0, s});
    while (!heap.empty()) {
      auto [du, u] = heap.top();
      heap.pop();
      if (du > dist[u]) continue;
      for (u32 v : adj[u]) {
        double nd = du + distance<D>(P, u, v);
        if (nd < dist[v]) { dist[v] = nd; heap.push({nd, v}); }
      }
    }
    for (u32 v = s + 1; v < n; ++v) worst = std::max(worst, dist[v] / distance<D>(P, s, v));
  }
  return worst;
}

template <int D> void check_tree(const Tree<D>& T, const PointSet<D>& P, const std::string& name) {
  std::vector<u8> seen(P.size(), 0);
  bool leaves_ok = true, partition_ok = true;
  for (u32 v = 0; v < T.node_count; ++v) {
    if (T.is_leaf(v)) {
      leaves_ok &= T.point[v] >= 0 && T.size[v] == 1 && !seen[T.point[v]];
      if (T.point[v] >= 0) seen[T.point[v]] = 1;
    } else {
      u32 s = 0;
      for (u32 w = T.child_begin[v]; w < T.child_end[v]; ++w) s += T.size[w];
      partition_ok &= s == T.size[v] && T.child_end[v] - T.child_begin[v] <= (2u << D);
    }
  }
  bool all_seen = true;
  for (u8 x : seen) all_seen &= x == 1;
  expect(leaves_ok && all_seen, name + ": one point per leaf");
  expect(partition_ok && T.size[0] == P.size(), name + ": children partition the points");
}

template <int D> void check_wspd(const Tree<D>& T, const PointSet<D>& P, SeparationTest test, double eps, const std::string& name) {
  const u32 n = P.size();
  std::vector<u8> cover((size_t)n * n, 0);
  for_each_wspd_pair<D>(T, eps, [&](u32 a, u32 b) {
    for (u32 p : T.points_of(a))
      for (u32 q : T.points_of(b)) {
        u32 x = std::min(p, q), y = std::max(p, q);
        if (cover[(size_t)x * n + y] < 255) ++cover[(size_t)x * n + y];
      }
  }, test);
  u64 wrong = 0;
  for (u32 p = 0; p < n; ++p)
    for (u32 q = p + 1; q < n; ++q) wrong += cover[(size_t)p * n + q] != 1;
  expect(wrong == 0, name + ": every point pair covered exactly once (eps " + std::to_string(eps) + (test == SeparationTest::Cell ? ", cell)" : ", box)"));
}

template <int D> void check_spanner(const Tree<D>& T, const PointSet<D>& P, double t, SeparationRule rule, Representative rep, const std::string& name) {
  SpannerParams params;
  params.t = t;
  params.rule = rule;
  params.representative = rep;
  auto edges = build_spanner<D>(T, params);
  double s = max_stretch<D>(P, edges);
  expect(s <= t * (1 + 1e-12), name + ": stretch " + std::to_string(s) + " <= t = " + std::to_string(t) +
                                   (rule == SeparationRule::HarPeled ? " (harpeled" : " (ck") + (rep == Representative::LightChild ? ", light)" : ", first)"));
}

template <int D> void run_case(const std::string& kind, u32 n, u64 seed, bool spanner) {
  auto P = PointSet<D>::from_real(sample<D>(kind, n, seed));
  if (P.size() < 2) return;
  const std::string name = kind + " d=" + std::to_string(D) + " n=" + std::to_string(P.size());
  auto level = build_level_quadtree<D>(P, {1.0 / (2 * D), 0.3});
  auto level_default = build_level_quadtree<D>(P);
  auto compressed = build_compressed_quadtree<D>(P);
  for (auto* T : {&level, &level_default, &compressed}) {
    const std::string tn = name + (T == &compressed ? " compressed" : T == &level ? " level(1/2d,0.3)" : " level");
    check_tree<D>(*T, P, tn);
    for (double eps : {1.0, 0.5, 0.25})
      for (SeparationTest test : {SeparationTest::BoundingBox, SeparationTest::Cell}) check_wspd<D>(*T, P, test, eps, tn);
    if (spanner)
      for (double t : {2.0, 1.5, 1.1}) {
        check_spanner<D>(*T, P, t, SeparationRule::CallahanKosaraju, Representative::LightChild, tn);
        check_spanner<D>(*T, P, t, SeparationRule::CallahanKosaraju, Representative::FirstLeaf, tn);
        if (t >= 1.5) check_spanner<D>(*T, P, t, SeparationRule::HarPeled, Representative::LightChild, tn);
      }
  }
}

void check_bounded_degree(const std::string& kind, u32 n, u64 seed) {
  auto P = PointSet<2>::from_real(sample<2>(kind, n, seed));
  const std::string name = "bounded degree " + kind + " n=" + std::to_string(P.size());
  auto trees = std::vector<Tree<2>>{build_level_quadtree<2>(P), build_compressed_quadtree<2>(P)};
  for (const auto& T : trees)
    for (double t : {2.0, 1.5, 1.25, 1.1}) {
      BoundedDegreeParams params;
      params.t = t;
      auto S = build_bounded_degree_spanner(P, T, params);
      const std::string tn = name + " t=" + std::to_string(t);
      double s = max_stretch<2>(P, S.edges);
      expect(s <= t * (1 + 1e-12), tn + ": stretch " + std::to_string(s));
      expect(S.max_degree == max_degree(P.size(), S.edges), tn + ": reported max degree");
      expect(S.max_degree <= S.degree_bound, tn + ": degree within D + (D+1)(kappa+1)");
      expect(S.max_sink_in_degree <= S.kappa + 1, tn + ": sink in-degree <= kappa + 1");
      bool simple = true;
      std::vector<std::pair<u32, u32>> sorted;
      for (const Edge& e : S.edges) { simple &= e.u != e.v; sorted.push_back({std::min(e.u, e.v), std::max(e.u, e.v)}); }
      std::sort(sorted.begin(), sorted.end());
      simple &= std::adjacent_find(sorted.begin(), sorted.end()) == sorted.end();
      expect(simple, tn + ": no loops or repeated edges");
    }
}

void check_point_set() {
  std::vector<std::array<double, 2>> x = {{0, 0}, {1, 1}, {0, 0}, {2, 3}, {1, 1}, {5, 0}};
  auto P = PointSet<2>::from_real(x);
  expect(P.size() == 4 && P.duplicates_removed == 2 && P.input_size == 6, "duplicates removed");
  std::vector<u32> idx(P.input_index.begin(), P.input_index.end());
  std::sort(idx.begin(), idx.end());
  expect(idx == std::vector<u32>{0, 1, 3, 5}, "first occurrences kept as input indices");
  bool keys_increasing = true;
  for (u32 i = 1; i < P.size(); ++i) keys_increasing &= P.keys[i - 1] < P.keys[i];
  expect(keys_increasing, "points in Morton order");
  auto throws = [](auto f) {
    try { f(); } catch (const std::invalid_argument&) { return true; }
    return false;
  };
  expect(throws([] { PointSet<2>::from_grid({{1ull << 20, 0}}, 20); }), "grid coordinate out of range throws");
  expect(throws([] { PointSet<3>::from_real({{0, 0, std::numeric_limits<double>::quiet_NaN()}}); }), "non-finite coordinate throws");
  expect(throws([] { PointSet<2>::from_grid({}, 63); }), "too many bits throws");
  expect(throws([&] { build_level_quadtree<2>(P, {0.75, 0.5}); }), "beta above 1/d throws");
  expect(throws([&] { build_level_quadtree<2>(PointSet<2>{}); }), "empty point set throws");
  auto T = build_level_quadtree<2>(P);
  expect(throws([&] { count_wspd_pairs<2>(T, 0); }), "eps = 0 throws");
  expect(throws([&] { build_spanner<2>(T, {1.0}); }), "t = 1 throws");
  expect(throws([&] { build_bounded_degree_spanner(P, T, {2.0, 0.5, 4}); }), "too few cones throws");
}

}

int main(int argc, char** argv) {
  const bool quick = argc > 1 && std::string(argv[1]) == "--quick";
  check_point_set();
  u64 seed = 1;
  const std::vector<u32> sizes = quick ? std::vector<u32>{50, 200} : std::vector<u32>{50, 300, 1200};
  for (u32 n : sizes)
    for (const char* kind : {"uniform", "clusters", "lattice", "line"}) {
      const bool spanner = n <= 300;
      run_case<2>(kind, n, ++seed, spanner);
      run_case<3>(kind, n, ++seed, spanner);
      run_case<4>(kind, n, ++seed, false);
      run_case<5>(kind, n, ++seed, false);
    }
  for (const char* kind : {"uniform", "clusters", "lattice"}) check_bounded_degree(kind, quick ? 150 : 400, ++seed);
  std::printf("%s: %d checks, %d failures\n", failures ? "FAILED" : "PASSED", checks, failures);
  return failures ? 1 : 0;
}

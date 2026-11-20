# lqt: Level Quadtrees, WSPDs, and geometric spanners

`lqt` is a header-only C++20 library that builds, for a set of points in 2 to 5
dimensions:

- the **Level Quadtree** of the paper below, and the standard **compressed quadtree**;
- the **well-separated pair decomposition (WSPD)** on either tree, for a given ε;
- a **t-spanner** with one edge per WSPD pair;
- a **bounded-degree t-spanner** (2D only).

It is the sequential (RAM-model) implementation of the paper. The paper's
Massively Parallel Computation (MPC) algorithms are not part of this library.
The package also has a command-line tool, `lqt-spanner`, that reads a point file
and writes the spanner edges or the WSPD pairs.

## Citation

M. A. Abam and M. R. Bahrami. *Level Quadtrees, Well-Separated Pair
Decompositions, and Geometric Spanners in Massively Parallel Computation.*
IEEE Access, under submission, 2026.

## Build

The library needs only a C++20 compiler (GCC or Clang; it uses `unsigned __int128`
and `__builtin_clzll`). To use it, add `include/` to the include path and include
`<lqt/lqt.hpp>`. To build the tool and the tests:

```
cmake -S . -B build
cmake --build build -j
cd build && ctest
```

The build type defaults to `Release` (`-O3`). Options:

| option | default | meaning |
|---|---|---|
| `LQT_NATIVE` | `OFF` | compile the tool and the tests with `-march=native` (or `-mcpu=native`) |
| `LQT_BUILD_TOOLS` | `ON` | build `lqt-spanner` |
| `LQT_BUILD_TESTS` | `ON` | build `lqt-tests` |

From another CMake project, add the library with `add_subdirectory(level-quadtree)`
and link `lqt::lqt`. After `cmake --install build`, use `find_package(lqt)` instead.
`lqt::lqt` is an INTERFACE target and adds no compile flags of its own.

`build/lqt-tests` runs the property tests in about 10 seconds; `lqt-tests --quick`
runs a smaller set.

## Library usage

```cpp
#include <lqt/lqt.hpp>

std::vector<std::array<double, 2>> input = ...;
auto points = lqt::PointSet<2>::from_real(input);
auto tree = lqt::build_level_quadtree<2>(points);

std::uint64_t pairs = lqt::count_wspd_pairs<2>(tree, 0.25);
lqt::for_each_wspd_pair<2>(tree, 0.25, [&](std::uint32_t a, std::uint32_t b) {
  for (std::uint32_t p : tree.points_of(a)) { ... }
});

lqt::SpannerParams params;
params.t = 1.5;
std::vector<lqt::Edge> edges = lqt::build_spanner<2>(tree, params);
edges = lqt::to_input_indices(points, std::move(edges));

lqt::for_each_spanner_edge<2>(tree, params, [&](std::uint32_t p, std::uint32_t q) { ... });

lqt::BoundedDegreeParams bounded;
bounded.t = 2;
lqt::BoundedDegreeSpanner sparse = lqt::build_bounded_degree_spanner(points, tree, bounded);
```

### Points

- `PointSet<D>::from_real(points, bits)` takes `std::array<double, D>` points.
- `PointSet<D>::from_grid(points, bits)` takes integer points `lqt::Point<D>`
  (`std::array<std::uint64_t, D>`) with every coordinate below `2^bits`.

Both sort the points in Morton (Z-) order and remove duplicates. Inside the library,
a point is identified by its rank in this order, its **point id**. Every function
reports points by point id. `points.input_index[id]` is the position of the point in
the input vector, and `to_input_indices` maps a whole edge list. Other members are
`points.points[id]` (the integer coordinates), `points.size()`, `points.input_size`,
and `points.duplicates_removed`.

### Trees

- `build_level_quadtree<D>(points, {beta, lambda}, &info)` builds the Level Quadtree.
  - `beta` is the grid exponent, in (0, 1/d]; the default 0 means 1/d.
  - `lambda` is the base exponent, in (0, 1]; the default is 1/2.
  - The optional `LevelQuadtreeInfo` receives the number of recursion levels, the
    non-base calls, the base calls, the forced base calls, and the donuts.
- `build_compressed_quadtree<D>(points)` builds the compressed quadtree.

Both return a `Tree<D>` without empty nodes. Its nodes are numbered from the root
`0` so that the children of node `v` are the consecutive ids
`child_begin[v] .. child_end[v] - 1`, and every child has a larger id than its parent.
Per node, the tree stores:

- the parent, the cell (`cell_corner`, `cell_level`, so the side is `2^cell_level`),
  and the bounding box of its points (`box_lo`, `box_hi`);
- the number of points (`size`) and the points themselves (`points_of(v)`);
- the depth and the `kind` (grid node, region, donut, compressed-quadtree node, or point).

A donut node's hole is given by `hole_corner` and `hole_level`.

### WSPD

`for_each_wspd_pair<D>(tree, eps, visitor, test)` calls `visitor(a, b)` once for
every pair of tree nodes in the WSPD and returns the number of pairs.
`count_wspd_pairs` only counts them. Every pair of distinct points lies in exactly
one pair `(a, b)`, with one point in `points_of(a)` and the other in `points_of(b)`.

A node pair is separated if max(Δ(a), Δ(b)) ≤ ε·δ(a, b). The `test` chooses Δ and δ:

- `SeparationTest::BoundingBox` (default, as in the paper): Δ is the diagonal of the
  bounding box of the node's points, and δ is the distance between the two boxes.
- `SeparationTest::Cell`: Δ is the diameter of the node's cell if it holds at least
  two points, and 0 for a single point. δ is the distance between the cells; a
  one-point node is its point.

The visitor is a template parameter and is inlined, so streaming the pairs costs
nothing extra. An optional last argument `std::vector<std::uint32_t>* split_counts`
receives, for every node, how often the algorithm split it (the paper's |V_x|).

### Spanners

`for_each_spanner_edge<D>(tree, params, visitor)` streams one edge `(p, q)` per
WSPD pair (bounding-box test), and `build_spanner` returns the edges as a vector.
There is one edge per pair, so the edges are distinct. Each edge is oriented from
the representative of the node that the last split introduced; this orientation
bounds the out-degree in the paper. `SpannerParams` has these fields:

- `t`: the stretch factor, greater than 1.
- `rule`: the separation for `t`.
  - `SeparationRule::CallahanKosaraju` (default) uses ε = (t−1)/(2(t+1)).
  - `SeparationRule::HarPeled` uses ε = (t−1)/16, Har-Peled's bound; it gives many
    more edges.

  Both give stretch at most t.
- `eps`: an explicit ε; a positive value overrides `rule`.
- `representative`: the endpoints of the edges.
  - `Representative::LightChild` (default, the paper's rule): a leaf is represented
    by its point, and an internal node by the representative of its child with the
    fewest points.
  - `Representative::FirstLeaf`: the leftmost leaf.

`max_degree(n, edges)` returns the largest degree of an edge list.

### Bounded-degree spanner (2D)

`build_bounded_degree_spanner(points, tree, params)` first builds a WSPD t'-spanner
with t' = t^`split` (default `split` = 1/2) and the representatives "last leaf of the
first child". It then replaces the edges that enter each point q by the q-sink
spanner on κ cones (Narasimhan and Smid, *Geometric Spanner Networks*, Sec. 4.2).
By default, κ is the smallest integer κ ≥ 9 with 1/(cos θ − sin θ) ≤ t/t', where
θ = 2π/κ; `params.kappa` fixes it instead. The result holds:

- the undirected, duplicate-free `edges`;
- `max_degree` and the theoretical `degree_bound` D + (D+1)(κ+1), where D is the
  largest distinct out-degree of the base spanner;
- `kappa`, `base_stretch`, `sink_stretch`, `base_eps`, `base_edges`, and
  `max_sink_in_degree`.

The base edges are stored, about 12 bytes each, while the sink spanners are built.

### Errors

Invalid arguments throw `std::invalid_argument`:

- `bits`, `beta`, `lambda`, `t`, `eps`, `split`, or `kappa` out of range;
- coordinates that are not finite or not below `2^bits`;
- an empty point set;
- `2^32 − 1` or more input points.

A dimension other than 2 to 5 is a compile-time error. The library never prints or exits.

## Command-line tool

```
lqt-spanner [options] [input]
```

The tool reads the points from `input`, or from standard input if `input` is
omitted or `-`. It writes the spanner edges to standard output, one edge `i j` per
line, and statistics to standard error: the points read and removed, the tree size,
the number of edges, and the maximum and average degree.

| option | meaning |
|---|---|
| `--t T` | stretch factor, default 2 |
| `--eps E` | separation; overrides the value derived from `T` |
| `--rule ck\|harpeled` | ε = (t−1)/(2(t+1)) (default) or (t−1)/16 |
| `--representative light\|first` | edge endpoints, default `light` |
| `--tree level\|compressed` | tree, default `level` |
| `--beta B`, `--lambda L` | Level Quadtree parameters, default 1/d and 1/2; fractions such as `1/4` are accepted |
| `--bounded-degree` | bounded-degree t-spanner (2D) |
| `--split A`, `--kappa K` | bounded-degree parameters, see above |
| `--d D` | dimension; default: the number of coordinates on the first line |
| `--bits K` | grid bits per axis |
| `--format auto\|plain\|tsp` | input format, default `auto` |
| `--edges FILE` | write the edges to `FILE` instead of standard output |
| `--no-edges` | write no edges, statistics only |
| `--pairs FILE` | also write the WSPD pairs (for the same ε) to `FILE` |

Examples:

```
lqt-spanner points.txt --t 1.5 > edges.txt
lqt-spanner usa115475.tsp --bounded-degree --t 2 --edges edges.txt
lqt-spanner points.txt --eps 0.25 --no-edges --pairs pairs.txt
```

### Input format

- **Plain text:** one point per line, with 2 to 5 coordinates separated by
  whitespace or commas. Blank lines and lines starting with `#` are skipped.
- **TSPLIB:** a file with a `NODE_COORD_SECTION` (for example `EUC_2D`) is read as
  2D TSPLIB. Each line `id x y` gives a point, up to `EOF`.

### Output format

Points are numbered by their position among the input points, from 0; skipped lines
do not count. A duplicate point keeps the number of its first occurrence.

- **Edges:** each line `i j` is one edge. For a t-spanner, the edge is oriented as in
  the library.
- **WSPD pairs:** each line `i1 i2 ... ; j1 j2 ...` lists the points of the two sets
  of one pair. This output grows with the total size of all sets, so it is meant for
  small inputs.

## Limits

- **Dimensions** 2 to 5. The bounded-degree spanner is 2D only.
- **Integer grid.** Coordinates are quantized to an integer grid with K bits per
  axis, so that every Morton key fits in 128 bits: K = 62, 42, 31, 25 for d = 2, 3, 4, 5
  (the maximum, and the default). A real point x becomes ⌊(x − min_i)/L · 2^K⌋, clamped
  to 2^K − 1, where min_i is the minimum on axis i and L the largest extent over the
  axes. All cells, boxes, and squared distances are then exact integers.
- **Resolution.** Points closer than about L/2^K fall into the same grid point and are
  merged as duplicates. Spreads up to about 2^K are represented exactly.
- **Duplicates** are removed after quantization; the first occurrence in input order
  is kept.
- **Size.** There are fewer than 2^32 − 1 points, with 32-bit ids.
- **Precision of the tests.** Separation tests use exact integer squared lengths.
  Each squared length is rounded once to double, and a double filter decides all
  pairs except those within a relative 10^−9 of the boundary, which take the exact
  path.

## Algorithms and complexity

The bounds below are the paper's, for constant d, β, and λ.

- **Level Quadtree.**
  - *Construction:* every call that is not a base call ranks its points along each
    axis and places a grid at the coordinates of evenly spaced ranks (every
    ⌈m^{1−β}⌉-th of the m points). Every region of the grid tree receives fewer than
    (2d+1)·m^{1−β} of the call's points.
  - *Base calls:* a call with at most n^λ points builds a compressed quadtree. So does
    a call that does not shrink its input, which is possible only for
    n < (2d+1)^{2/(βλ)}.
  - *Bounds:* there are at most ⌈2 ln(1/λ)/β⌉ levels of non-base calls. Construction
    takes O(n log n) time, the tree has O(n) nodes, and its height is
    O(n^β + n^λ).
  - *Implementation:* the per-axis orders are sorted once by radix sort and inherited
    by stable partition; the grid is never materialized.
- **Compressed quadtree.** It is built in linear time from the Morton-sorted points,
  after an O(n) radix sort of the keys.
- **WSPD.** This is Har-Peled's algorithm on the tree. The initial pairs are all pairs
  of distinct children of every node. A pair that is not separated is split at the
  node with the larger Δ. On equal Δ, the node is split that has the lower
  construction level, then the smaller depth inside its call's tree, then the smaller
  smallest input index. The output has O(n/ε^d) pairs and takes O(n log n + n/ε^d)
  time. The pair loop uses an explicit stack, so there is no recursion.
- **Spanner.** It has one edge per pair, so O(n/ε^d) edges, built in the same time as
  the WSPD.
- **Bounded-degree spanner.**
  - *Base spanner:* a WSPD t'-spanner with O(n) edges, whose orientation has bounded
    out-degree D.
  - *Sink spanners:* each in-star is replaced by a sink spanner on κ cones.
  - *Result:* a t-spanner with maximum degree at most D + (D+1)(κ+1) and O(n) edges.
- **Memory.** The trees use O(n) memory. WSPD pairs and t-spanner edges are streamed
  to the visitor and stored only by `build_spanner`.

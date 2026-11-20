#include "lqt/lqt.hpp"

#include <charconv>
#include <cstdlib>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <iterator>
#include <memory>
#include <sstream>
#include <string>

namespace {

const char* const usage_text =
    "usage: lqt-spanner [options] [input]\n"
    "  input                points file (plain or TSPLIB); standard input if omitted or -\n"
    "  --t T                stretch factor, T > 1 (default 2)\n"
    "  --eps E              WSPD separation; overrides the value derived from T\n"
    "  --rule ck|harpeled   separation derived from T: (T-1)/(2(T+1)) or (T-1)/16 (default ck)\n"
    "  --representative light|first   edge endpoints of a pair (default light)\n"
    "  --tree level|compressed        quadtree (default level)\n"
    "  --beta B --lambda L  Level Quadtree parameters (default 1/d and 1/2)\n"
    "  --bounded-degree     bounded-degree t-spanner (d = 2 only)\n"
    "  --split A            bounded degree: base stretch T^A (default 0.5)\n"
    "  --kappa K            bounded degree: number of cones (default: smallest valid)\n"
    "  --d D                dimension 2..5 (default: from the first point)\n"
    "  --bits K             grid bits per axis (default 62, 42, 31, 25 for d = 2..5)\n"
    "  --format auto|plain|tsp        input format (default auto)\n"
    "  --edges FILE         write the spanner edges to FILE (default -, standard output)\n"
    "  --no-edges           do not write edges, print statistics only\n"
    "  --pairs FILE         also write the WSPD pairs to FILE\n"
    "  --help               show this help\n";

struct Options {
  std::string input = "-", format = "auto", edges = "-", pairs, tree = "level", rule = "ck", representative = "light";
  double t = 2, eps = 0, beta = 0, lambda = 0.5, split = 0.5;
  int d = 0, bits = 0;
  unsigned kappa = 0;
  bool bounded_degree = false, write_edges = true;
};

struct UsageError : std::runtime_error {
  using std::runtime_error::runtime_error;
};

double parse_number(const std::string& s) {
  auto slash = s.find('/');
  try {
    if (slash != std::string::npos) return std::stod(s.substr(0, slash)) / std::stod(s.substr(slash + 1));
    return std::stod(s);
  } catch (const std::exception&) {
    throw UsageError("not a number: " + s);
  }
}

Options parse_options(int argc, char** argv) {
  Options o;
  bool have_input = false;
  for (int i = 1; i < argc; ++i) {
    std::string a = argv[i];
    auto value = [&]() -> std::string {
      if (i + 1 >= argc) throw UsageError("missing value for " + a);
      return argv[++i];
    };
    if (a == "--help" || a == "-h") { std::cout << usage_text; std::exit(0); }
    else if (a == "--t") o.t = parse_number(value());
    else if (a == "--eps") o.eps = parse_number(value());
    else if (a == "--rule") o.rule = value();
    else if (a == "--representative") o.representative = value();
    else if (a == "--tree") o.tree = value();
    else if (a == "--beta") o.beta = parse_number(value());
    else if (a == "--lambda") o.lambda = parse_number(value());
    else if (a == "--bounded-degree") o.bounded_degree = true;
    else if (a == "--split") o.split = parse_number(value());
    else if (a == "--kappa") o.kappa = (unsigned)parse_number(value());
    else if (a == "--d") o.d = (int)parse_number(value());
    else if (a == "--bits") o.bits = (int)parse_number(value());
    else if (a == "--format") o.format = value();
    else if (a == "--edges") o.edges = value();
    else if (a == "--no-edges") o.write_edges = false;
    else if (a == "--pairs") o.pairs = value();
    else if (a.size() > 1 && a[0] == '-' && a != "-") throw UsageError("unknown option " + a);
    else if (have_input) throw UsageError("more than one input file");
    else { o.input = a; have_input = true; }
  }
  if (o.rule != "ck" && o.rule != "harpeled") throw UsageError("--rule must be ck or harpeled");
  if (o.representative != "light" && o.representative != "first") throw UsageError("--representative must be light or first");
  if (o.tree != "level" && o.tree != "compressed") throw UsageError("--tree must be level or compressed");
  if (o.format != "auto" && o.format != "plain" && o.format != "tsp") throw UsageError("--format must be auto, plain or tsp");
  if (o.d && (o.d < 2 || o.d > 5)) throw UsageError("--d must be 2, 3, 4 or 5");
  if (o.bounded_degree && o.d && o.d != 2) throw UsageError("--bounded-degree needs d = 2");
  return o;
}

std::string read_all(const std::string& path) {
  if (path == "-") return std::string(std::istreambuf_iterator<char>(std::cin), std::istreambuf_iterator<char>());
  std::ifstream f(path, std::ios::binary);
  if (!f) throw std::runtime_error("cannot open " + path);
  return std::string(std::istreambuf_iterator<char>(f), std::istreambuf_iterator<char>());
}

struct RawPoints {
  int d = 0;
  std::vector<double> coordinates;
  size_t count() const { return d ? coordinates.size() / d : 0; }
};

std::vector<double> numbers_in(const std::string& line) {
  std::vector<double> v;
  const char* p = line.c_str();
  while (true) {
    while (*p == ' ' || *p == '\t' || *p == '\r' || *p == ',') ++p;
    if (!*p) break;
    char* end;
    double x = std::strtod(p, &end);
    if (end == p) throw std::runtime_error("cannot parse line: " + line);
    v.push_back(x);
    p = end;
  }
  return v;
}

RawPoints parse_plain(const std::string& text, int d) {
  RawPoints r;
  r.d = d;
  std::istringstream in(text);
  std::string line;
  size_t line_number = 0;
  while (std::getline(in, line)) {
    ++line_number;
    size_t first = line.find_first_not_of(" \t\r");
    if (first == std::string::npos || line[first] == '#') continue;
    auto v = numbers_in(line);
    if (!r.d) r.d = (int)v.size();
    if ((int)v.size() != r.d)
      throw std::runtime_error("line " + std::to_string(line_number) + ": expected " + std::to_string(r.d) + " coordinates");
    r.coordinates.insert(r.coordinates.end(), v.begin(), v.end());
  }
  return r;
}

RawPoints parse_tsp(const std::string& text, int d) {
  if (d && d != 2) throw UsageError("TSPLIB input is two-dimensional");
  RawPoints r;
  r.d = 2;
  std::istringstream in(text);
  std::string line;
  bool in_section = false;
  while (std::getline(in, line)) {
    if (!in_section) {
      if (line.find("NODE_COORD_SECTION") != std::string::npos) in_section = true;
      continue;
    }
    if (line.find("EOF") != std::string::npos) break;
    std::istringstream ss(line);
    double id, x, y;
    if (!(ss >> id >> x >> y)) continue;
    r.coordinates.push_back(x);
    r.coordinates.push_back(y);
  }
  return r;
}

class Writer {
 public:
  explicit Writer(const std::string& path) {
    if (path == "-") file_ = stdout;
    else {
      file_ = std::fopen(path.c_str(), "w");
      if (!file_) throw std::runtime_error("cannot write " + path);
      owned_ = true;
    }
    buffer_.reserve(capacity + 64);
  }
  ~Writer() {
    if (!buffer_.empty()) std::fwrite(buffer_.data(), 1, buffer_.size(), file_);
    if (owned_) std::fclose(file_);
    else std::fflush(file_);
  }
  Writer(const Writer&) = delete;
  Writer& operator=(const Writer&) = delete;
  void number(lqt::u64 x) {
    char tmp[24];
    auto res = std::to_chars(tmp, tmp + sizeof tmp, x);
    buffer_.append(tmp, res.ptr);
    if (buffer_.size() >= capacity) flush();
  }
  void put(char c) { buffer_.push_back(c); }
  void flush() {
    if (!buffer_.empty() && std::fwrite(buffer_.data(), 1, buffer_.size(), file_) != buffer_.size())
      throw std::runtime_error("write failed");
    buffer_.clear();
  }

 private:
  static constexpr size_t capacity = 1 << 20;
  std::FILE* file_ = nullptr;
  bool owned_ = false;
  std::string buffer_;
};

template <int D> int run(const Options& o, const RawPoints& raw) {
  std::vector<std::array<double, D>> real(raw.count());
  for (size_t j = 0; j < real.size(); ++j)
    for (int i = 0; i < D; ++i) real[j][i] = raw.coordinates[j * D + i];
  const int bits = o.bits ? o.bits : lqt::max_bits_per_axis<D>();
  auto points = lqt::PointSet<D>::from_real(real, bits);
  real.clear();
  real.shrink_to_fit();
  std::fprintf(stderr, "points: %llu read, %u distinct, %llu duplicates removed, d = %d, %d bits per axis\n",
               (unsigned long long)points.input_size, points.size(), (unsigned long long)points.duplicates_removed, D, bits);
  if (points.size() < 2) throw std::runtime_error("need at least two distinct points");
  lqt::Tree<D> tree;
  if (o.tree == "level") {
    lqt::LevelQuadtreeInfo info;
    tree = lqt::build_level_quadtree<D>(points, {o.beta, o.lambda}, &info);
    std::fprintf(stderr, "tree: Level Quadtree, %u nodes, height %u, %d recursion levels, %llu non-base calls\n", tree.node_count,
                 tree.height, info.levels, (unsigned long long)info.calls);
  } else {
    tree = lqt::build_compressed_quadtree<D>(points);
    std::fprintf(stderr, "tree: compressed quadtree, %u nodes, height %u\n", tree.node_count, tree.height);
  }
  lqt::SpannerParams params;
  params.t = o.t;
  params.eps = o.eps;
  params.rule = o.rule == "harpeled" ? lqt::SeparationRule::HarPeled : lqt::SeparationRule::CallahanKosaraju;
  params.representative = o.representative == "first" ? lqt::Representative::FirstLeaf : lqt::Representative::LightChild;
  const double eps = lqt::spanner_separation(params);
  const auto& index = points.input_index;
  if (!o.pairs.empty()) {
    Writer w(o.pairs);
    auto write_set = [&](lqt::u32 v) {
      bool first = true;
      for (lqt::u32 p : tree.points_of(v)) {
        if (!first) w.put(' ');
        first = false;
        w.number(index[p]);
      }
    };
    lqt::u64 pairs = lqt::for_each_wspd_pair<D>(tree, eps, [&](lqt::u32 a, lqt::u32 b) {
      write_set(a);
      w.put(' '); w.put(';'); w.put(' ');
      write_set(b);
      w.put('\n');
    });
    w.flush();
    std::fprintf(stderr, "wspd: %llu pairs, eps = %.6g\n", (unsigned long long)pairs, eps);
  }
  std::vector<lqt::u32> degree(points.size(), 0);
  lqt::u64 edges = 0;
  std::unique_ptr<Writer> out;
  if (o.write_edges) out = std::make_unique<Writer>(o.edges);
  auto emit = [&](lqt::u32 p, lqt::u32 q) {
    ++degree[p];
    ++degree[q];
    if (out) { out->number(index[p]); out->put(' '); out->number(index[q]); out->put('\n'); }
  };
  if (o.bounded_degree) {
    if constexpr (D == 2) {
      lqt::BoundedDegreeParams bp;
      bp.t = o.t;
      bp.split = o.split;
      bp.kappa = o.kappa;
      auto result = lqt::build_bounded_degree_spanner(points, tree, bp);
      for (const lqt::Edge& e : result.edges) emit(e.u, e.v);
      edges = result.edges.size();
      std::fprintf(stderr, "bounded degree: base stretch %.6g (eps %.6g, %llu base edges), %u cones, sink stretch %.6g\n",
                   result.base_stretch, result.base_eps, (unsigned long long)result.base_edges, result.kappa, result.sink_stretch);
    } else {
      throw UsageError("--bounded-degree needs d = 2");
    }
  } else {
    edges = lqt::for_each_spanner_edge<D>(tree, params, emit);
    std::fprintf(stderr, "spanner: t = %.6g, eps = %.6g\n", o.t, eps);
  }
  if (out) out->flush();
  lqt::u32 max_degree = 0;
  for (lqt::u32 x : degree) max_degree = std::max(max_degree, x);
  std::fprintf(stderr, "edges: %llu, max degree %u, average degree %.3f\n", (unsigned long long)edges, max_degree,
               2.0 * (double)edges / points.size());
  return 0;
}

}

int main(int argc, char** argv) {
  try {
    Options o = parse_options(argc, argv);
    std::string text = read_all(o.input);
    bool tsp = o.format == "tsp" || (o.format == "auto" && text.find("NODE_COORD_SECTION") != std::string::npos);
    RawPoints raw = tsp ? parse_tsp(text, o.d) : parse_plain(text, o.d);
    text.clear();
    text.shrink_to_fit();
    if (raw.count() == 0) throw std::runtime_error("no points in the input");
    if (raw.d < 2 || raw.d > 5) throw std::runtime_error("points must have 2 to 5 coordinates");
    if (o.bounded_degree && raw.d != 2) throw UsageError("--bounded-degree needs d = 2");
    switch (raw.d) {
      case 2: return run<2>(o, raw);
      case 3: return run<3>(o, raw);
      case 4: return run<4>(o, raw);
      default: return run<5>(o, raw);
    }
  } catch (const UsageError& e) {
    std::fprintf(stderr, "lqt-spanner: %s\n%s", e.what(), usage_text);
    return 2;
  } catch (const std::exception& e) {
    std::fprintf(stderr, "lqt-spanner: %s\n", e.what());
    return 1;
  }
}

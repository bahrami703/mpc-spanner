#pragma once

#include "core.hpp"

#include <limits>
#include <stdexcept>
#include <string>

namespace lqt {

template <int D> struct PointSet {
  static_assert(D >= 2 && D <= 5, "lqt supports dimensions 2 to 5");

  int bits = 0;
  std::vector<Point<D>> points;
  std::vector<u128> keys;
  std::vector<u32> input_index;
  u64 input_size = 0;
  u64 duplicates_removed = 0;

  u32 size() const { return (u32)points.size(); }

  static PointSet from_grid(const std::vector<Point<D>>& grid_points, int bits = max_bits_per_axis<D>()) {
    if (bits < 1 || bits > max_bits_per_axis<D>())
      throw std::invalid_argument("lqt: bits per axis must be in [1, " + std::to_string(max_bits_per_axis<D>()) + "]");
    if (grid_points.size() >= std::numeric_limits<u32>::max())
      throw std::invalid_argument("lqt: too many points for 32-bit point ids");
    const u64 limit = 1ull << bits;
    for (const Point<D>& p : grid_points)
      for (int i = 0; i < D; ++i)
        if (p[i] >= limit) throw std::invalid_argument("lqt: grid coordinate not below 2^bits");
    PointSet s;
    s.bits = bits;
    s.sort_and_deduplicate(grid_points);
    return s;
  }

  static PointSet from_real(const std::vector<std::array<double, D>>& real_points, int bits = max_bits_per_axis<D>()) {
    return from_grid(quantize(real_points, bits), bits);
  }

  static std::vector<Point<D>> quantize(const std::vector<std::array<double, D>>& x, int bits = max_bits_per_axis<D>()) {
    if (bits < 1 || bits > max_bits_per_axis<D>())
      throw std::invalid_argument("lqt: bits per axis must be in [1, " + std::to_string(max_bits_per_axis<D>()) + "]");
    std::vector<Point<D>> q(x.size());
    if (x.empty()) return q;
    for (const auto& p : x)
      for (int i = 0; i < D; ++i)
        if (!std::isfinite(p[i])) throw std::invalid_argument("lqt: coordinate is not finite");
    std::array<double, D> lo, hi;
    for (int i = 0; i < D; ++i) lo[i] = hi[i] = x[0][i];
    for (const auto& p : x)
      for (int i = 0; i < D; ++i) { lo[i] = std::min(lo[i], p[i]); hi[i] = std::max(hi[i], p[i]); }
    double extent = 0;
    for (int i = 0; i < D; ++i) extent = std::max(extent, hi[i] - lo[i]);
    const u64 top = bits >= 64 ? ~0ull : ((1ull << bits) - 1);
    const double scale = extent > 0 ? std::ldexp(1.0, bits) / extent : 0.0;
    for (size_t j = 0; j < x.size(); ++j)
      for (int i = 0; i < D; ++i) {
        double t = (x[j][i] - lo[i]) * scale;
        q[j][i] = t <= 0 ? 0 : (t >= (double)top ? top : (u64)t);
      }
    return q;
  }

 private:
  struct KeyRecord {
    u64 lo, hi;
    u32 index;
  };

  void sort_and_deduplicate(const std::vector<Point<D>>& in) {
    const size_t n0 = in.size();
    input_size = n0;
    std::vector<KeyRecord> r(n0);
    for (size_t i = 0; i < n0; ++i) {
      u128 k = detail::morton_key<D>(in[i]);
      r[i] = {(u64)k, (u64)(k >> 64), (u32)i};
    }
    if (n0 < 4096) {
      std::stable_sort(r.begin(), r.end(), [](const KeyRecord& a, const KeyRecord& b) {
        return a.hi != b.hi ? a.hi < b.hi : a.lo < b.lo;
      });
    } else {
      detail::radix_sort(r, D * bits, [](const KeyRecord& x, int p) -> u32 {
        int s = 16 * p;
        return s < 64 ? (u32)((x.lo >> s) & 0xffff) : (u32)((x.hi >> (s - 64)) & 0xffff);
      });
    }
    std::vector<u8> duplicate(n0, 0);
    u64 dups = 0;
    for (size_t i = 1; i < n0; ++i)
      if (r[i].lo == r[i - 1].lo && r[i].hi == r[i - 1].hi) { duplicate[r[i].index] = 1; ++dups; }
    const size_t n = n0 - dups;
    duplicates_removed = dups;
    points.resize(n);
    keys.resize(n);
    input_index.resize(n);
    size_t j = 0;
    for (size_t i = 0; i < n0; ++i) {
      if (duplicate[r[i].index]) continue;
      points[j] = in[r[i].index];
      keys[j] = ((u128)r[i].hi << 64) | r[i].lo;
      input_index[j] = r[i].index;
      ++j;
    }
  }
};

}

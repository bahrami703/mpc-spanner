#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <vector>

namespace lqt {

using u8 = std::uint8_t;
using u16 = std::uint16_t;
using u32 = std::uint32_t;
using i32 = std::int32_t;
using u64 = std::uint64_t;
using i64 = std::int64_t;
using u128 = unsigned __int128;

template <int D> using Point = std::array<u64, D>;

template <int D> constexpr int max_bits_per_axis() { return D == 2 ? 62 : D == 3 ? 42 : D == 4 ? 31 : 25; }

namespace detail {

inline int bit_length(u64 x) { return x ? 64 - __builtin_clzll(x) : 0; }

inline int highest_bit(u128 x) {
  u64 hi = (u64)(x >> 64);
  return hi ? 127 - __builtin_clzll(hi) : 63 - __builtin_clzll((u64)x);
}

inline double to_double(u128 x) {
  u64 hi = (u64)(x >> 64);
  if (!hi) return (double)(u64)x;
  int shift = 64 - __builtin_clzll(hi);
  u64 lo = (u64)x;
  u64 top = (u64)(x >> shift);
  u64 mask = shift == 64 ? ~0ull : ((1ull << shift) - 1);
  top |= (lo & mask) != 0;
  u64 bits = (u64)(1023 + shift) << 52;
  double scale;
  std::memcpy(&scale, &bits, 8);
  return (double)top * scale;
}

inline u64 spread_bits_2d(u64 x) {
  x &= 0xffffffffull;
  x = (x | (x << 16)) & 0x0000FFFF0000FFFFull;
  x = (x | (x << 8)) & 0x00FF00FF00FF00FFull;
  x = (x | (x << 4)) & 0x0F0F0F0F0F0F0F0Full;
  x = (x | (x << 2)) & 0x3333333333333333ull;
  x = (x | (x << 1)) & 0x5555555555555555ull;
  return x;
}

inline u64 spread_bits_3d(u64 x) {
  x &= 0x1fffffull;
  x = (x | (x << 32)) & 0x1f00000000ffffull;
  x = (x | (x << 16)) & 0x1f0000ff0000ffull;
  x = (x | (x << 8)) & 0x100f00f00f00f00full;
  x = (x | (x << 4)) & 0x10c30c30c30c30c3ull;
  x = (x | (x << 2)) & 0x1249249249249249ull;
  return x;
}

template <int D> inline u128 morton_key(const Point<D>& p) {
  if constexpr (D == 2) {
    u64 lo = spread_bits_2d(p[0]) | (spread_bits_2d(p[1]) << 1);
    u64 hi = spread_bits_2d(p[0] >> 32) | (spread_bits_2d(p[1] >> 32) << 1);
    return ((u128)hi << 64) | lo;
  } else if constexpr (D == 3) {
    u64 lo = spread_bits_3d(p[0]) | (spread_bits_3d(p[1]) << 1) | (spread_bits_3d(p[2]) << 2);
    u64 hi = spread_bits_3d(p[0] >> 21) | (spread_bits_3d(p[1] >> 21) << 1) | (spread_bits_3d(p[2] >> 21) << 2);
    return ((u128)hi << 63) | lo;
  } else {
    u128 k = 0;
    for (int b = 0; b < 128 / D; ++b)
      for (int i = 0; i < D; ++i) k |= (u128)((p[i] >> b) & 1) << (D * b + i);
    return k;
  }
}

template <int D> inline int common_cell_level(u128 a, u128 b) { return highest_bit(a ^ b) / D + 1; }

template <int D> inline Point<D> cell_corner_of(const Point<D>& p, int level) {
  Point<D> c;
  for (int i = 0; i < D; ++i) c[i] = level >= 64 ? 0 : (p[i] >> level) << level;
  return c;
}

inline u64 ceil_power(u64 m, double a) {
  double ma = std::pow((double)m, a);
  long long k = std::max(1LL, (long long)std::ceil(ma - 1e-9));
  while (k - 1 >= 1 && (double)(k - 1) >= ma - 1e-12) --k;
  while ((double)k < ma - 1e-9) ++k;
  return (u64)k;
}

template <class Record, class DigitOf>
void radix_sort(std::vector<Record>& a, int bits, DigitOf digit_of) {
  const size_t n = a.size();
  if (n < 2) return;
  const int passes = (bits + 15) / 16;
  std::vector<Record> b(n);
  std::vector<u32> histogram((size_t)passes << 16, 0);
  for (size_t i = 0; i < n; ++i)
    for (int p = 0; p < passes; ++p) ++histogram[((size_t)p << 16) + digit_of(a[i], p)];
  Record* src = a.data();
  Record* dst = b.data();
  for (int p = 0; p < passes; ++p) {
    u32* h = &histogram[(size_t)p << 16];
    bool single_bucket = false;
    for (u32 d = 0; d < 65536; ++d)
      if (h[d]) { single_bucket = (h[d] == n); break; }
    if (single_bucket) continue;
    u32 sum = 0;
    for (u32 d = 0; d < 65536; ++d) { u32 c = h[d]; h[d] = sum; sum += c; }
    for (size_t i = 0; i < n; ++i) dst[h[digit_of(src[i], p)]++] = src[i];
    std::swap(src, dst);
  }
  if (src != a.data()) std::memcpy(a.data(), src, n * sizeof(Record));
}

}

}

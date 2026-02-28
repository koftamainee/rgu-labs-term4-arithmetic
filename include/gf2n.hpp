#ifndef RGU_LABS_TERM4_ARITHMETIC_GFN_HPP
#define RGU_LABS_TERM4_ARITHMETIC_GFN_HPP
#include <cstdint>
#include <cstddef>
#include <stdexcept>
#include <string>
#include "polynomial.hpp"
class GF2n final {
public:
  using elem = uint64_t;
  GF2n(int n, elem mod) :m_n(n), m_mod(mod) {
    if (n <= 1 || n > 64) {
      throw std::invalid_argument("invalid n");
    }
  }

  elem add(elem a, elem b) const {
    return a ^ b;
  }

  elem mul(elem a, elem b) const {
    elem res = 0;
    const elem top_bit = static_cast<elem>(1) << (m_n - 1);
    while (b != 0) {
      if (b & 1) {
        res ^= a;
      }
      b >>= 1;
      const bool carry = (a & top_bit) != 0;
      a <<= 1;
      a &= (top_bit << 1) - 1;
      if (carry) {
        a ^= (m_mod & ((top_bit << 1) - 1));
      }
    }
    return res;
  }
  elem divmod(elem a, elem b, elem &r) const {
    r = a;
    elem q = 0;
    const int deg_b = degree(b);
    while (degree(r) >= deg_b) {
      const int shift = degree(r) - deg_b;
      q ^= static_cast<elem>(1) << shift;
      r ^= b << shift;
    }
    return q;
  }
  elem egcd(elem a, elem b, elem &x, elem &y) const {
    if (b == 0) {
      x = 1;
      y = 0;
      return a;
    }
    elem x1;
    elem y1;
    elem r;
    elem q = divmod(a, b, r);
    const elem g = egcd(b, r, x1, y1);
    x = y1;
    y = x1 ^ mul(q, y1);
    return g;
  }
  elem inv(elem a) const {
    if (a == 0) {
      throw std::invalid_argument("zero has no inverse");
    }
    elem x, y;
    const elem g = egcd(a, m_mod, x, y);
    if (g != 1) {
      throw std::runtime_error("no inverse");
    }
    return x;
  }
  int degree(elem a) const {
    if (a == 0) {
      return -1;
    }
    int deg = 0;
    if (a >> 32) { a >>= 32; deg += 32; }
    if (a >> 16) { a >>= 16; deg += 16; }
    if (a >> 8)  { a >>= 8;  deg += 8;  }
    if (a >> 4)  { a >>= 4;  deg += 4;  }
    if (a >> 2)  { a >>= 2;  deg += 2;  }
    if (a >> 1)  {           deg += 1;  }
    return deg;
  }
  elem from_polynomial(const Polynomial &p) {
    const Vector &v = p.coefficients();
    elem a = 0;
    const size_t n = v.dimension();
    for (size_t i = 0; i < n && i < 64; i++) {
      if (v[i] != 0) {
        a |= (static_cast<elem>(1) << i);
      }
    }
    return a;
  }
  Polynomial to_polynomial(elem a) const {
    std::vector<bigfloat> coeffs;
    // Include up to degree m_n (the leading term of the modulus)
    for (int i = 0; i <= m_n; i++) {
      coeffs.emplace_back((a >> i) & 1 ? 1 : 0);
    }
    return {Vector(coeffs)};
  }
  std::string to_string(elem a) const {
    return to_polynomial(a).to_string();
  }
private:
  int m_n;
  elem m_mod;
};
#endif //RGU_LABS_TERM4_ARITHMETIC_GFN_HPP
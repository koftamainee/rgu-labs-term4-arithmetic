#pragma once

#include <cstdint>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <memory>
#include <vector>

#include "poly_ring.hpp"
#include "monomial.hpp"
#include "polynomial.hpp"

class GF2n final {
public:
  using elem = uint64_t;
  using size_type = std::size_t;

  GF2n(int n, elem mod) : m_n(n), m_mod(mod) {
    if (n <= 1 || n > 64) {
      throw std::invalid_argument("GF2n::GF2n: n must be in range (1, 64]");
    }
  }

  int n() const noexcept { return m_n; }
  elem mod() const noexcept { return m_mod; }

  elem add(elem a, elem b) const noexcept {
    return a ^ b;
  }

  elem mul(elem a, elem b) const noexcept {
    elem res = 0;
    const elem top_bit = static_cast<elem>(1) << (m_n - 1);
    const elem mask = (top_bit << 1) - 1;
    while (b != 0) {
      if (b & 1) {
        res ^= a;
      }
      b >>= 1;
      const bool carry = (a & top_bit) != 0;
      a <<= 1;
      a &= mask;
      if (carry) {
        a ^= (m_mod & mask);
      }
    }
    return res;
  }

  elem divmod(elem a, elem b, elem& r) const {
    if (b == 0) {
      throw std::invalid_argument("GF2n::divmod: division by zero");
    }
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

  elem inv(elem a) const {
    if (a == 0) {
      throw std::invalid_argument("GF2n::inv: zero has no inverse");
    }
    elem x = 0;
    elem y = 0;
    const elem g = egcd(a, m_mod, x, y);
    if (g != 1) {
      throw std::runtime_error("GF2n::inv: element has no inverse");
    }
    return x;
  }

  int degree(elem a) const noexcept {
    if (a == 0) {
      return -1;
    }
    int deg = 0;
    if (a >> 32) {
      a >>= 32;
      deg += 32;
    }
    if (a >> 16) {
      a >>= 16;
      deg += 16;
    }
    if (a >> 8) {
      a >>= 8;
      deg += 8;
    }
    if (a >> 4) {
      a >>= 4;
      deg += 4;
    }
    if (a >> 2) {
      a >>= 2;
      deg += 2;
    }
    if (a >> 1) { deg += 1; }
    return deg;
  }

  elem from_polynomial(const Polynomial<int>& p) const {
    elem a = 0;
    const auto ring = p.ring();
    const size_type vars = p.n_vars();
    if (vars != 1) {
      throw std::invalid_argument("GF2n::from_polynomial: polynomial must be univariate");
    }
    for (const auto& mono : p.supp()) {
      const int exp = mono[0];
      if (exp < 0 || exp >= 64) {
        continue;
      }
      const int* coeff = p.get(mono);
      if (coeff != nullptr && (*coeff % 2) != 0) {
        a |= (static_cast<elem>(1) << exp);
      }
    }
    return a;
  }

  Polynomial<int> to_polynomial(elem a, const std::shared_ptr<PolyRing<int>> &ring) const {
    if (ring->n_vars() != 1) {
      throw std::invalid_argument("Invalid ring for gf2n polynomial");
    }
    Polynomial<int> p(ring);
    for (int i = 0; i <= m_n; ++i) {
      if ((a >> i) & 1) {
        Monomial mono(1);
        mono[0] = i;
        p.set(mono, 1);
      }
    }
    return p;
  }

  std::string to_string(elem a) const {
    return to_polynomial(a).to_string();
  }

  std::string to_latex(elem a) const {
    return to_polynomial(a).to_latex();
  }

private:
  elem egcd(elem a, elem b, elem& x, elem& y) const {
    if (b == 0) {
      x = 1;
      y = 0;
      return a;
    }
    elem x1 = 0;
    elem y1 = 0;
    elem r = 0;
    const elem q = divmod(a, b, r);
    const elem g = egcd(b, r, x1, y1);
    x = y1;
    y = x1 ^ mul(q, y1);
    return g;
  }

  int m_n;
  elem m_mod;
};

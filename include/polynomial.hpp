#ifndef POLYNOMIAL_HPP
#define POLYNOMIAL_HPP

#include "bigfloat.h"
#include "vector.h"
#include <string>

#include "poly_tostring.hpp"

class Polynomial {
private:
  Vector coeffs_;
  bigfloat a_;

public:
  Polynomial(const Vector &coeffs, const bigfloat &a = bigfloat(0))
      : coeffs_(coeffs), a_(a) {}

  const Vector &coefficients() const { return coeffs_; }

  const bigfloat &expansion_point() const { return a_; }

  size_t degree() const {
    size_t deg = coeffs_.dimension();
    if (deg == 0) {
      return 0;
    }

    deg--;
    while (deg > 0 && coeffs_[deg] == bigfloat(0)) {
      deg--;
    }
    return deg;
  }

  bigfloat evaluate(const bigfloat &x) const {
    size_t n = coeffs_.dimension();
    if (n == 0) {
      return bigfloat(0);
    }

    bigfloat result = coeffs_[n - 1];
    bigfloat dx = x - a_;

    for (int i = static_cast<int>(n) - 2; i >= 0; i--) {
      result = result * dx + coeffs_[i];
    }
    return result;
  }

  Polynomial change_expansion_point(const bigfloat &B) const {
    if (B == a_) {
      return *this;
    }

    Vector new_coeffs = coeffs_;
    bigfloat c = B - a_;
    size_t n = new_coeffs.dimension();

    for (size_t i = 0; i < n; i++) {
      for (size_t j = n - 1; j > i; j--) {
        new_coeffs[j - 1] = new_coeffs[j - 1] + c * new_coeffs[j];
      }
    }

    return Polynomial(new_coeffs, B);
  }

  size_t zero_order() const {
    size_t n = coeffs_.dimension();
    for (size_t i = 0; i < n; i++) {
      if (coeffs_[i] != 0) {
        return i;
      }
    }
    return n;
  }

  bool is_zero() const { return zero_order() == coeffs_.dimension(); }

  std::string to_string() const { return poly_tostring(coeffs_, a_); }
};

#endif

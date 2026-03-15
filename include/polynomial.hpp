#ifndef POLYNOMIAL_HPP
#define POLYNOMIAL_HPP

#include "bigmath/bigfloat.hpp"
#include "VectorBF.h"
#include <string>

#include "poly_tostring.hpp"

class Polynomial {
private:
  VectorBF coeffs_;
  bigfloat a_;

public:
  Polynomial(const VectorBF &coeffs, const bigfloat &a = bigfloat(0))
      : coeffs_(coeffs), a_(a) {}

  const VectorBF &coefficients() const { return coeffs_; }

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

    VectorBF new_coeffs = coeffs_;
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

  Polynomial operator+(const Polynomial& other) const {
    const size_t n = std::max(coeffs_.dimension(), other.coeffs_.dimension());
    std::vector<bigfloat> result(n, bigfloat(0));
    for (size_t i = 0; i < coeffs_.dimension(); ++i) result[i] += coeffs_[i];
    for (size_t i = 0; i < other.coeffs_.dimension(); ++i) result[i] += other.coeffs_[i];
    return Polynomial(VectorBF(result));
  }

  Polynomial operator-(const Polynomial& other) const {
    const size_t n = std::max(coeffs_.dimension(), other.coeffs_.dimension());
    std::vector<bigfloat> result(n, bigfloat(0));
    for (size_t i = 0; i < coeffs_.dimension(); ++i) result[i] += coeffs_[i];
    for (size_t i = 0; i < other.coeffs_.dimension(); ++i) result[i] -= other.coeffs_[i];
    return Polynomial(VectorBF(result));
  }

  Polynomial operator*(const bigfloat& scalar) const {
    std::vector<bigfloat> result(coeffs_.dimension());
    for (size_t i = 0; i < coeffs_.dimension(); ++i) result[i] = coeffs_[i] * scalar;
    return Polynomial(VectorBF(result));
  }

  Polynomial operator*(const Polynomial& other) const {
    const size_t fn = coeffs_.dimension();
    const size_t gn = other.coeffs_.dimension();
    if (fn == 0 || gn == 0)
      return Polynomial(VectorBF(std::vector<bigfloat>{bigfloat(0)}));
    std::vector<bigfloat> result(fn + gn - 1, bigfloat(0));
    for (size_t i = 0; i < fn; ++i)
      for (size_t j = 0; j < gn; ++j)
        result[i + j] += coeffs_[i] * other.coeffs_[j];
    return Polynomial(VectorBF(result));
  }

  Polynomial rem_xn_minus_1(const size_t n) const {
    std::vector<bigfloat> result(n, bigfloat(0));
    for (size_t i = 0; i < coeffs_.dimension(); ++i) result[i % n] += coeffs_[i];
    return Polynomial(VectorBF(result));
  }

  Polynomial rem_xn_plus_1(const size_t n) const {
    std::vector<bigfloat> result(n, bigfloat(0));
    for (size_t i = 0; i < coeffs_.dimension(); ++i) {
      if ((i / n) % 2 == 0) result[i % n] += coeffs_[i];
      else                   result[i % n] -= coeffs_[i];
    }
    return Polynomial(VectorBF(result));
  }

  Polynomial substitute_wx(const std::vector<bigfloat>& omega_powers) const {
    const size_t n = coeffs_.dimension();
    std::vector<bigfloat> result(n);
    for (size_t i = 0; i < n; ++i)
      result[i] = coeffs_[i] * omega_powers[i % omega_powers.size()];
    return Polynomial(VectorBF(result));
  }

  Polynomial first_n_coeffs(const size_t n) const {
    std::vector<bigfloat> result(n, bigfloat(0));
    for (size_t i = 0; i < n && i < coeffs_.dimension(); ++i)
      result[i] = coeffs_[i];
    return Polynomial(VectorBF(result));
  }

  std::string to_string() const { return poly_tostring(coeffs_, a_); }
};

#endif
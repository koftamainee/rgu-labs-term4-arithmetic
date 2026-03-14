#pragma once

#include <stdexcept>
#include <vector>

#include "bigmath/bigfloat.hpp"

class PowerSeries {
private:
  std::vector<bigfloat> coeffs_;

public:
  PowerSeries() = default;

  explicit PowerSeries(std::vector<bigfloat> coeffs)
    : coeffs_(std::move(coeffs)) {}

  explicit PowerSeries(size_t n)
    : coeffs_(n, bigfloat(0)) {}

  size_t size() const noexcept { return coeffs_.size(); }

  bigfloat& operator[](size_t i) { return coeffs_[i]; }
  const bigfloat& operator[](size_t i) const { return coeffs_[i]; }

  PowerSeries first_n(size_t n) const {
    PowerSeries result(n);
    for (size_t i = 0; i < n && i < coeffs_.size(); i++)
      result[i] = coeffs_[i];
    return result;
  }

  PowerSeries operator+(const PowerSeries& other) const {
    const size_t n = std::max(coeffs_.size(), other.coeffs_.size());
    PowerSeries result(n);
    for (size_t i = 0; i < coeffs_.size(); i++) result[i] += coeffs_[i];
    for (size_t i = 0; i < other.coeffs_.size(); i++) result[i] += other.coeffs_[i];
    return result;
  }

  PowerSeries operator-(const PowerSeries& other) const {
    const size_t n = std::max(coeffs_.size(), other.coeffs_.size());
    PowerSeries result(n);
    for (size_t i = 0; i < coeffs_.size(); ++i) result[i] += coeffs_[i];
    for (size_t i = 0; i < other.coeffs_.size(); ++i) result[i] -= other.coeffs_[i];
    return result;
  }

  PowerSeries operator*(const PowerSeries& other) const {
    if (coeffs_.empty() || other.coeffs_.empty())
      return PowerSeries(std::vector<bigfloat>{bigfloat(0)});
    const size_t n = coeffs_.size() + other.coeffs_.size() - 1;
    PowerSeries result(n);
    for (size_t i = 0; i < coeffs_.size(); ++i)
      for (size_t j = 0; j < other.coeffs_.size(); ++j)
        result[i + j] += coeffs_[i] * other.coeffs_[j];
    return result;
  }

  PowerSeries operator*(const bigfloat& scalar) const {
    PowerSeries result(coeffs_.size());
    for (size_t i = 0; i < coeffs_.size(); ++i)
      result[i] = coeffs_[i] * scalar;
    return result;
  }

  PowerSeries inverse(size_t n) const {
    if (coeffs_.empty() || coeffs_[0] == bigfloat(0))
      throw std::domain_error("PowerSeries::inverse: constant term is zero");

    PowerSeries g(n);
    const bigfloat inv_f0 = coeffs_[0].reciprocal();
    g[0] = inv_f0;

    for (size_t i = 1; i < n; ++i) {
      bigfloat sum(0);
      for (size_t k = 1; k <= i; ++k) {
        bigfloat fk = k < coeffs_.size() ? coeffs_[k] : bigfloat(0);
        sum += fk * g[i - k];
      }
      g[i] = -inv_f0 * sum;
    }

    return g;
  }

  const std::vector<bigfloat>& coefficients() const { return coeffs_; }
};

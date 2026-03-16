#pragma once

#include <algorithm>
#include <cstddef>
#include <stdexcept>
#include <vector>

#include "fft.hpp"

class PowerSeries {
private:
  std::vector<double> coeffs_;

public:
  PowerSeries() = default;

  explicit PowerSeries(std::vector<double> coeffs)
    : coeffs_(std::move(coeffs)) {}

  explicit PowerSeries(size_t n)
    : coeffs_(n, 0.0) {}

  size_t size() const noexcept { return coeffs_.size(); }

  double& operator[](size_t i) { return coeffs_[i]; }
  const double& operator[](size_t i) const { return coeffs_[i]; }

  const std::vector<double>& coefficients() const { return coeffs_; }

  PowerSeries first_n(size_t n) const {
    PowerSeries result(n);
    for (size_t i = 0; i < n && i < coeffs_.size(); ++i)
      result[i] = coeffs_[i];
    return result;
  }

  PowerSeries operator+(const PowerSeries& other) const {
    const size_t n = std::max(coeffs_.size(), other.coeffs_.size());
    PowerSeries result(n);
    for (size_t i = 0; i < coeffs_.size(); ++i) result[i] += coeffs_[i];
    for (size_t i = 0; i < other.coeffs_.size(); ++i) result[i] += other.coeffs_[i];
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
      return PowerSeries(std::vector<double>{0.0});

    const size_t result_size = coeffs_.size() + other.coeffs_.size() - 1;

    size_t n = 1;
    while (n < result_size)
      n <<= 1;

    std::vector<CD> f(n), g(n);

    for (size_t i = 0; i < coeffs_.size(); ++i)
      f[i] = coeffs_[i];

    for (size_t i = 0; i < other.coeffs_.size(); ++i)
      g[i] = other.coeffs_[i];

    fft(f, false);
    fft(g, false);

    for (size_t i = 0; i < n; ++i)
      f[i] *= g[i];

    fft(f, true);

    std::vector<double> result(result_size);

    for (size_t i = 0; i < result_size; ++i)
      result[i] = f[i].real() / static_cast<double>(n);

    return PowerSeries(result);
  }

  PowerSeries operator*(double scalar) const {
    PowerSeries result(coeffs_.size());
    for (size_t i = 0; i < coeffs_.size(); ++i)
      result[i] = coeffs_[i] * scalar;
    return result;
  }

  friend PowerSeries operator*(double scalar, const PowerSeries& ps) {
    return ps * scalar;
  }

  PowerSeries inverse(size_t n) const {
    if (coeffs_.empty() || coeffs_[0] == 0.0)
      throw std::domain_error("PowerSeries::inverse: constant term is zero");

    PowerSeries g(1);
    g[0] = 1.0 / coeffs_[0];
    size_t m = 1;
    while (m < n) {
      size_t m2 = std::min(2 * m, n);
      PowerSeries f_cut = this->first_n(m2);
      PowerSeries fg = (f_cut * g).first_n(m2);
      PowerSeries two(m2);
      two[0] = 2.0;
      PowerSeries correction = two - fg;
      g = (g * correction).first_n(m2);
      m = m2;
    }
    return g.first_n(n);
  }
};

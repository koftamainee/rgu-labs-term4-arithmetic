#pragma once

#include <algorithm>
#include <cstddef>
#include <stdexcept>
#include <vector>

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

    double &operator[](size_t i) { return coeffs_[i]; }
    const double &operator[](size_t i) const { return coeffs_[i]; }

    const std::vector<double> &coefficients() const { return coeffs_; }

    PowerSeries first_n(size_t n) const {
        PowerSeries result(n);
        for (size_t i = 0; i < n && i < coeffs_.size(); ++i)
            result[i] = coeffs_[i];
        return result;
    }

    PowerSeries operator+(const PowerSeries &other) const {
        const size_t n = std::max(coeffs_.size(), other.coeffs_.size());
        PowerSeries result(n);
        for (size_t i = 0; i < coeffs_.size(); ++i)       result[i] += coeffs_[i];
        for (size_t i = 0; i < other.coeffs_.size(); ++i) result[i] += other.coeffs_[i];
        return result;
    }

    PowerSeries operator-(const PowerSeries &other) const {
        const size_t n = std::max(coeffs_.size(), other.coeffs_.size());
        PowerSeries result(n);
        for (size_t i = 0; i < coeffs_.size(); ++i)       result[i] += coeffs_[i];
        for (size_t i = 0; i < other.coeffs_.size(); ++i) result[i] -= other.coeffs_[i];
        return result;
    }

    PowerSeries operator*(const PowerSeries &other) const {
        if (coeffs_.empty() || other.coeffs_.empty())
            return PowerSeries(std::vector<double>{0.0});
        const size_t n = coeffs_.size() + other.coeffs_.size() - 1;
        PowerSeries result(n);
        for (size_t i = 0; i < coeffs_.size(); ++i)
            for (size_t j = 0; j < other.coeffs_.size(); ++j)
                result[i + j] += coeffs_[i] * other.coeffs_[j];
        return result;
    }

    PowerSeries operator*(double scalar) const {
        PowerSeries result(coeffs_.size());
        for (size_t i = 0; i < coeffs_.size(); ++i)
            result[i] = coeffs_[i] * scalar;
        return result;
    }

    friend PowerSeries operator*(double scalar, const PowerSeries &ps) {
        return ps * scalar;
    }

    PowerSeries inverse(size_t n) const {
        if (coeffs_.empty() || coeffs_[0] == 0.0)
            throw std::domain_error("PowerSeries::inverse: constant term is zero");
        PowerSeries g(n);
        const double inv_f0 = 1.0 / coeffs_[0];
        g[0] = inv_f0;
        for (size_t i = 1; i < n; ++i) {
            double sum = 0.0;
            for (size_t k = 1; k <= i; ++k) {
                double fk = k < coeffs_.size() ? coeffs_[k] : 0.0;
                sum += fk * g[i - k];
            }
            g[i] = -inv_f0 * sum;
        }
        return g;
    }
};
#pragma once

#include <cmath>
#include <cstddef>
#include <functional>
#include <vector>

#include "../primitives/vector.hpp"
#include "../matrix.hpp"

template <typename T>
struct SystemResult {
  Vector<T> solution;
  std::size_t iterations;
  bool converged;
};

template <typename T>
struct SystemStep {
  std::size_t iteration;
  Vector<T> approximation;
};

template <typename T>
SystemResult<T> newton_system(
  const std::function<Vector<T>(const Vector<T>&)>& f,
  const std::function<Matrix<T>(const Vector<T>&)>& j,
  const Vector<T>& x0,
  T eps = T{1e-10},
  std::size_t max_iter = 100000,
  std::vector<SystemStep<T>>* steps = nullptr
) {
  Vector<T> x = x0;
  for (std::size_t i = 0; i < max_iter; ++i) {
    Vector<T> fx = f(x);
    Vector<T> delta(j(x).solve_gauss(fx.components()));
    Vector<T> x_new = x - delta;
    if (steps != nullptr) {
      steps->push_back({i + 1, x_new});
    }
    if (f(x_new).norm() < eps && delta.norm() < eps) {
      return {x_new, i + 1, true};
    }
    x = x_new;
  }
  return {x, max_iter, false};
}

template <typename T>
SystemResult<T> newton_system_inf(
  const std::function<Vector<T>(const Vector<T>&)>& f,
  const std::function<Matrix<T>(const Vector<T>&)>& j,
  const Vector<T>& x0,
  T eps = T{1e-10},
  std::size_t max_iter = 100000,
  std::vector<SystemStep<T>>* steps = nullptr
) {
  Vector<T> x = x0;
  for (std::size_t i = 0; i < max_iter; ++i) {
    Vector<T> fx = f(x);
    Matrix<T> j_inv = j(x).inverse();
    const std::size_t n = x.dimension();
    std::vector<T> delta_raw(n, T{});
    for (std::size_t row = 0; row < n; ++row) {
      for (std::size_t col = 0; col < n; ++col) {
        delta_raw[row] += j_inv.at(row, col) * fx.components()[col];
      }
    }
    Vector<T> delta(delta_raw);
    Vector<T> x_new = x - delta;
    if (steps != nullptr) {
      steps->push_back({i + 1, x_new});
    }
    T inf_norm{};
    for (std::size_t k = 0; k < delta.dimension(); ++k) {
      inf_norm = std::fmax(inf_norm, std::fabs(delta[k]));
    }
    if (inf_norm < eps) {
      return {x_new, i + 1, true};
    }
    x = x_new;
  }
  return {x, max_iter, false};
}

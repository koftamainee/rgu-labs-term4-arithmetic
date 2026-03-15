#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <vector>

#include "vector.hpp"
#include "matrix.hpp"

struct SystemResult {
  Vector solution;
  size_t iterations;
  bool converged;
};

struct SystemStep {
  size_t iteration;
  Vector approximation;
};

inline SystemResult newton_system(
  const std::function<Vector(const Vector&)>& F,
  const std::function<Matrix(const Vector&)>& J,
  const Vector& x0,
  double eps = 1e-10,
  size_t max_iter = 100000,
  std::vector<SystemStep>* steps = nullptr
) {
  Vector x = x0;
  for (size_t i = 0; i < max_iter; ++i) {
    Vector fx = F(x);
    std::vector<double> fx_raw(fx.components().begin(), fx.components().end());
    std::vector<double> delta_raw = J(x).solve_gauss(fx_raw);
    Vector delta(delta_raw);
    Vector x_new = x - delta;
    if (steps)
      steps->push_back({i + 1, x_new});
    if (F(x_new).norm() < eps && delta.norm() < eps)
      return {x_new, i + 1, true};
    x = x_new;
  }
  return {x, max_iter, false};
}

inline SystemResult newton_system_inf(

  const std::function<Vector(const Vector&)>& F,

  const std::function<Matrix(const Vector&)>& J,

  const Vector& x0,

  double eps = 1e-10,
  size_t max_iter = 100000,
  std::vector<SystemStep>* steps = nullptr
) {
  Vector x = x0;
  for (size_t i = 0; i < max_iter; ++i) {
    Vector fx = F(x);
    std::vector<double> fx_raw(fx.components().begin(), fx.components().end());
    Matrix J_inv = J(x).inverse();
    const size_t n = x.dimension();
    std::vector<double> delta_raw(n, 0.0);
    for (size_t row = 0; row < n; ++row)
      for (size_t col = 0; col < n; ++col)
        delta_raw[row] += J_inv.at(row, col) * fx_raw[col];
    Vector delta(delta_raw);
    Vector x_new = x - delta;
    if (steps)
      steps->push_back({i + 1, x_new});
    double inf_norm = 0.0;
    for (size_t k = 0; k < delta.dimension(); ++k)
      inf_norm = std::fmax(inf_norm, std::fabs(delta[k]));
    if (inf_norm < eps)
      return {x_new, i + 1, true};
    x = x_new;
  }
  return {x, max_iter, false};
}

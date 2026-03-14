#pragma once

#include <functional>
#include <stdexcept>
#include <vector>

#include "bigmath/bigfloat.hpp"
#include "matrix.h"
#include "vector.h"

struct SystemResult {
  Vector solution;
  size_t iterations;
  bool   converged;
};

struct SystemStep {
  size_t iteration;
  Vector approximation;
};

inline SystemResult newton_system(
    const std::function<Vector(const Vector&)>& F,
    const std::function<Matrix(const Vector&)>& J,
    const Vector&   x0,
    const bigfloat& eps      = bigfloat::DEFAULT_EPS,
    size_t          max_iter = 100000,
    std::vector<SystemStep>* steps = nullptr)
{
  Vector x = x0;

  for (size_t i = 0; i < max_iter; ++i) {
    Vector fx = F(x);

    std::vector<bigfloat> fx_raw(fx.components().begin(), fx.components().end());
    std::vector<bigfloat> delta_raw = J(x).solve_gauss(fx_raw);

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
    std::function<Vector(const Vector&)>  F,
    std::function<Matrix(const Vector&)>  J,
    const Vector&   x0,
    const bigfloat& eps      = bigfloat::DEFAULT_EPS,
    size_t          max_iter = 100000,
    std::vector<SystemStep>* steps = nullptr)
{
  Vector x = x0;

  for (size_t i = 0; i < max_iter; ++i) {
    Vector fx = F(x);

    std::vector<bigfloat> fx_raw(fx.components().begin(), fx.components().end());
    Matrix J_inv = J(x).inverse();

    std::vector<bigfloat> delta_raw(x.dimension());
    for (size_t row = 0; row < x.dimension(); ++row) {
      bigfloat sum(0);
      for (size_t col = 0; col < x.dimension(); ++col)
        sum += J_inv.at(row, col) * fx_raw[col];
      delta_raw[row] = sum;
    }

    Vector delta(delta_raw);
    Vector x_new = x - delta;

    if (steps)
      steps->push_back({i + 1, x_new});

    bigfloat inf_norm(0);
    for (size_t k = 0; k < delta.dimension(); ++k) {
      bigfloat d = delta[k].abs();
      if (d > inf_norm) inf_norm = d;
    }

    if (inf_norm < eps)
      return {x_new, i + 1, true};

    x = x_new;
  }

  return {x, max_iter, false};
}
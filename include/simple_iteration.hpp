#pragma once

#include <functional>

#include "bigmath/bigfloat.hpp"

struct IterationResult {
  bigfloat root;
  size_t   iterations;
  bool     converged;
};

inline IterationResult simple_iteration(
    const std::function<bigfloat(const bigfloat&)>& phi,
    const bigfloat& x0,
    const bigfloat& eps      = bigfloat::DEFAULT_EPS,
    size_t   max_iter = 100000)
{
  bigfloat x = x0;
  for (size_t i = 0; i < max_iter; i++) {
    bigfloat x_new = phi(x);
    if ((x_new - x).abs() < eps)
      return {x_new, i + 1, true};
    x = x_new;
  }
  return {x, max_iter, false};
}
#pragma once

#include <functional>
#include <cmath>

struct IterationResult {
  double root;
  size_t iterations;
  bool converged;
};


inline IterationResult simple_iteration(
  const std::function<double(double)>& phi,
  double x0,
  double eps = 1e-10,
  size_t max_iter = 100000) {
  double x = x0;
  for (size_t i = 0; i < max_iter; ++i) {
    const double x_new = phi(x);
    if (std::fabs(x_new - x) < eps) {
      return {x_new, i + 1, true};
    }
    x = x_new;
  }
  return {x, max_iter, false};
}

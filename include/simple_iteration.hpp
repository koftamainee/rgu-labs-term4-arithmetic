#pragma once

#include <cmath>
#include <cstddef>
#include <functional>

template <typename T>
struct IterationResult {
  T root;
  std::size_t iterations;
  bool converged;
};

template <typename T>
IterationResult<T> simple_iteration(
    const std::function<T(T)>& phi,
    T x0,
    T eps = T{1e-10},
    std::size_t max_iter = 100000
) {
  T x = x0;
  for (std::size_t i = 0; i < max_iter; ++i) {
    T x_new = phi(x);
    if (std::fabs(x_new - x) < eps) {
      return {x_new, i + 1, true};
    }
    x = x_new;
  }
  return {x, max_iter, false};
}
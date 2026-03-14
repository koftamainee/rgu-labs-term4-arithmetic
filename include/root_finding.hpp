#pragma once

#include <functional>
#include <stdexcept>
#include <vector>

#include "bigmath/bigfloat.hpp"

struct RootResult {
  bigfloat root;
  size_t   iterations;
  bool     converged;
};

struct RootStep {
  size_t   iteration;
  bigfloat approximation;
};

inline RootResult bisection(
    std::function<bigfloat(const bigfloat&)> f,
    const bigfloat& a,
    const bigfloat& b,
    const bigfloat& eps      = bigfloat::DEFAULT_EPS,
    size_t          max_iter = 100000,
    std::vector<RootStep>* steps = nullptr)
{
  bigfloat lo = a;
  bigfloat hi = b;

  if ((f(lo) > bigfloat(0)) == (f(hi) > bigfloat(0)))
    throw std::invalid_argument("bisection: f(a) and f(b) must have opposite signs");

  for (size_t i = 0; i < max_iter; ++i) {
    bigfloat mid = (lo + hi) / bigfloat(2);

    if (steps)
      steps->push_back({i + 1, mid});

    if ((hi - lo) / bigfloat(2) < eps)
      return {mid, i + 1, true};

    if ((f(mid) > bigfloat(0)) == (f(lo) > bigfloat(0)))
      lo = mid;
    else
      hi = mid;
  }

  return {(lo + hi) / bigfloat(2), max_iter, false};
}

inline RootResult newton(
    std::function<bigfloat(const bigfloat&)> f,
    std::function<bigfloat(const bigfloat&)> df,
    const bigfloat& x0,
    const bigfloat& eps      = bigfloat::DEFAULT_EPS,
    size_t          max_iter = 100000,
    std::vector<RootStep>* steps = nullptr)
{
  bigfloat x = x0;

  for (size_t i = 0; i < max_iter; ++i) {
    bigfloat fx  = f(x);
    bigfloat dfx = df(x);

    if (dfx == bigfloat(0))
      throw std::runtime_error("newton: derivative is zero");

    bigfloat x_new = x - fx / dfx;

    if (steps)
      steps->push_back({i + 1, x_new});

    if ((x_new - x).abs() < eps)
      return {x_new, i + 1, true};

    x = x_new;
  }

  return {x, max_iter, false};
}

inline RootResult newton_modified(
    std::function<bigfloat(const bigfloat&)> f,
    std::function<bigfloat(const bigfloat&)> df,
    const bigfloat& x0,
    const bigfloat& sigma,
    const bigfloat& eps      = bigfloat::DEFAULT_EPS,
    size_t          max_iter = 100000,
    std::vector<RootStep>* steps = nullptr)
{
  bigfloat x = x0;

  for (size_t i = 0; i < max_iter; ++i) {
    bigfloat fx  = f(x);
    bigfloat dfx = df(x);

    if (dfx == bigfloat(0))
      throw std::runtime_error("newton_modified: derivative is zero");

    bigfloat x_new = x - sigma * fx / dfx;

    if (steps)
      steps->push_back({i + 1, x_new});

    if ((x_new - x).abs() < eps)
      return {x_new, i + 1, true};

    x = x_new;
  }

  return {x, max_iter, false};
}
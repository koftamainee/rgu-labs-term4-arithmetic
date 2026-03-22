#pragma once

#include <cmath>
#include <cstddef>
#include <functional>
#include <stdexcept>
#include <vector>

template <typename T>
struct RootResult {
    T root;
    std::size_t iterations;
    bool converged;
};

template <typename T>
struct RootStep {
    std::size_t iteration;
    T approximation;
};

template <typename T>
RootResult<T> bisection(
    const std::function<T(T)>& f,
    T a,
    T b,
    T eps = T{1e-10},
    std::size_t max_iter = 100000,
    std::vector<RootStep<T>>* steps = nullptr
) {
    T lo = a;
    T hi = b;
    if ((f(lo) > T{}) == (f(hi) > T{})) {
        throw std::invalid_argument("bisection: f(a) and f(b) must have opposite signs");
    }
    for (std::size_t i = 0; i < max_iter; ++i) {
        T mid = (lo + hi) / T{2};
        if (steps != nullptr) {
            steps->push_back({i + 1, mid});
        }
        if ((hi - lo) / T{2} < eps) {
            return {mid, i + 1, true};
        }
        if ((f(mid) > T{}) == (f(lo) > T{})) {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    return {(lo + hi) / T{2}, max_iter, false};
}

template <typename T>
RootResult<T> newton_modified(
    const std::function<T(T)>& f,
    const std::function<T(T)>& df,
    T x0,
    T sigma,
    T eps = T{1e-10},
    std::size_t max_iter = 100000,
    std::vector<RootStep<T>>* steps = nullptr
) {
    T x = x0;
    for (std::size_t i = 0; i < max_iter; ++i) {
        T fx = f(x);
        T dfx = df(x);
        if (dfx == T{}) {
            throw std::runtime_error("newton_modified: derivative is zero");
        }
        T x_new = x - sigma * fx / dfx;
        if (steps != nullptr) {
            steps->push_back({i + 1, x_new});
        }
        if (std::fabs(x_new - x) < eps) {
            return {x_new, i + 1, true};
        }
        x = x_new;
    }
    return {x, max_iter, false};
}

template <typename T>
RootResult<T> newton(
    const std::function<T(T)>& f,
    const std::function<T(T)>& df,
    T x0,
    T eps = T{1e-10},
    std::size_t max_iter = 100000,
    std::vector<RootStep<T>>* steps = nullptr
) {
    return newton_modified<T>(f, df, x0, T{1}, eps, max_iter, steps);
}
#pragma once

#include <cmath>
#include <cstddef>
#include <functional>
#include <stdexcept>
#include <vector>

struct RootResult {
    double root;
    size_t iterations;
    bool   converged;
};

struct RootStep {
    size_t iteration;
    double approximation;
};

inline RootResult bisection(
    std::function<double(double)> f,
    double a,
    double b,
    double eps      = 1e-10,
    size_t max_iter = 100000,
    std::vector<RootStep> *steps = nullptr)
{
    double lo = a;
    double hi = b;
    if ((f(lo) > 0.0) == (f(hi) > 0.0))
        throw std::invalid_argument("bisection: f(a) and f(b) must have opposite signs");
    for (size_t i = 0; i < max_iter; ++i) {
        double mid = (lo + hi) / 2.0;
        if (steps)
            steps->push_back({i + 1, mid});
        if ((hi - lo) / 2.0 < eps)
            return {mid, i + 1, true};
        if ((f(mid) > 0.0) == (f(lo) > 0.0))
            lo = mid;
        else
            hi = mid;
    }
    return {(lo + hi) / 2.0, max_iter, false};
}

inline RootResult newton(
    std::function<double(double)> f,
    std::function<double(double)> df,
    double x0,
    double eps      = 1e-10,
    size_t max_iter = 100000,
    std::vector<RootStep> *steps = nullptr)
{
    double x = x0;
    for (size_t i = 0; i < max_iter; ++i) {
        double fx  = f(x);
        double dfx = df(x);
        if (dfx == 0.0)
            throw std::runtime_error("newton: derivative is zero");
        double x_new = x - fx / dfx;
        if (steps)
            steps->push_back({i + 1, x_new});
        if (std::fabs(x_new - x) < eps)
            return {x_new, i + 1, true};
        x = x_new;
    }
    return {x, max_iter, false};
}

inline RootResult newton_modified(
    std::function<double(double)> f,
    std::function<double(double)> df,
    double x0,
    double sigma,
    double eps      = 1e-10,
    size_t max_iter = 100000,
    std::vector<RootStep> *steps = nullptr)
{
    double x = x0;
    for (size_t i = 0; i < max_iter; ++i) {
        double fx  = f(x);
        double dfx = df(x);
        if (dfx == 0.0)
            throw std::runtime_error("newton_modified: derivative is zero");
        double x_new = x - sigma * fx / dfx;
        if (steps)
            steps->push_back({i + 1, x_new});
        if (std::fabs(x_new - x) < eps)
            return {x_new, i + 1, true};
        x = x_new;
    }
    return {x, max_iter, false};
}
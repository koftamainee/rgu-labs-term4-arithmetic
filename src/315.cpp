#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

#include "simple_iteration.hpp"

void print_result(double x0, const IterationResult &r) {
    std::cout << "    x0=" << std::fixed << std::setprecision(2) << std::setw(7) << x0 << "  =>  ";
    if (r.converged)
        std::cout << "converged in " << std::setw(6) << r.iterations
                  << " iter,  root = " << std::setprecision(10) << r.root << "\n";
    else
        std::cout << "did NOT converge after " << r.iterations << " iter\n";
}

int main() {
    constexpr double EPS = 1e-10;

    std::cout << "=== a) x_{n+1} = 2^(x_n - 1) ===\n";
    std::cout << "    phi'(x) = ln2 * 2^(x-1),  |phi'| < 1 iff x < 1 + log2(1/ln2) ~ 1.528\n";
    std::cout << "    Fixed point: x ~ 1.0  (converges for x0 < ~1.528)\n\n";

    auto phi_a = [](double x) { return std::pow(2.0, x - 1.0); };
    for (double x0 : {-2.0, 0.0, 0.5, 1.0, 1.3, 1.5, 1.52, 2.0, 5.0})
        print_result(x0, simple_iteration(phi_a, x0, EPS));
    std::cout << "\n";

    std::cout << "=== b) x_{n+1} = e^(2x_n) - 1 ===\n";
    std::cout << "    phi'(x) = 2*e^(2x) >= 2 > 1 for all x  =>  diverges for ANY x0\n\n";

    auto phi_b = [](double x) { return std::exp(2.0 * x) - 1.0; };
    for (double x0 : {-5.0, -1.0, -0.1, 0.0, 0.1, 1.0})
        print_result(x0, simple_iteration(phi_b, x0, EPS, 200));
    std::cout << "\n";

    std::cout << "=== c) x_{n+1} = A - ln(x_n) ===\n";
    std::cout << "    phi'(x) = -1/x,  |phi'| < 1 iff x > 1\n";
    std::cout << "    Converges for x0 > 1 (fixed point depends on A)\n\n";

    struct CaseC { double A; std::vector<double> x0s; };
    std::vector<CaseC> cases_c = {
        { 2.0, {0.1, 0.5, 1.5, 2.0, 5.0, 10.0} },
        { 3.0, {0.1, 0.5, 1.5, 3.0, 8.0       } },
    };
    for (const auto &c : cases_c) {
        std::cout << "  A=" << std::fixed << std::setprecision(1) << c.A << "\n";
        auto phi_c = [&](double x) { return c.A - std::log(x); };
        for (double x0 : c.x0s)
            print_result(x0, simple_iteration(phi_c, x0, EPS));
        std::cout << "\n";
    }

    std::cout << "=== d) x_{n+1} = alpha*e^(-x_n) + beta*x_n ===\n";
    std::cout << "    phi'(x) = -alpha*e^(-x) + beta\n";
    std::cout << "    |phi'| < 1 requires |beta| < 1 and |alpha| < (1-|beta|)*e^x near root\n\n";

    struct CaseD { double alpha, beta; std::vector<double> x0s; };
    std::vector<CaseD> cases_d = {
        { 0.5,  0.3, {-2.0, 0.0, 0.5, 1.0, 3.0} },
        { 1.0,  0.5, {-1.0, 0.0, 0.5, 2.0      } },
        { 1.0, -0.5, { 0.0, 0.5, 1.0, 3.0      } },
        { 2.0,  0.9, { 0.0, 0.5, 1.0            } },
        { 1.0,  1.5, { 0.0, 1.0                 } },
    };
    for (const auto &c : cases_d) {
        std::cout << std::fixed << std::setprecision(4)
                  << "  alpha=" << c.alpha << "  beta=" << c.beta
                  << "  |beta|=" << std::fabs(c.beta)
                  << "  =>  " << (std::fabs(c.beta) < 1.0 ? "potentially CONVERGES" : "DIVERGES") << "\n";
        auto phi_d = [&](double x) { return c.alpha * std::exp(-x) + c.beta * x; };
        for (double x0 : c.x0s)
            print_result(x0, simple_iteration(phi_d, x0, EPS));
        std::cout << "\n";
    }

    return 0;
}
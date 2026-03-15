#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

#include "root_finding.hpp"

static const double PI = std::acos(-1.0);

auto f = [](double x) {
    return std::pow(x - 1.0, 3) * std::sin(PI * x) * (std::cos(2.0 * PI * x) - 1.0);
};

auto df = [](double x) {
    double a  = std::pow(x - 1.0, 3),  da = 3.0 * (x - 1.0) * (x - 1.0);
    double b  = std::sin(PI * x),       db = PI * std::cos(PI * x);
    double c  = std::cos(2.0*PI*x)-1.0, dc = -2.0 * PI * std::sin(2.0 * PI * x);
    return da*b*c + a*db*c + a*b*dc;
};

struct Root { double x0; double exact; int p; };

int main() {
    const int    n   = 6;
    const double eps = std::pow(10.0, -static_cast<double>(n));

    std::cout << "f(x) = (x-1)^3 * sin(pi*x) * (cos(2*pi*x) - 1)\n\n";
    std::cout << "Roots and their multiplicities:\n";
    std::cout << "  x=1: (x-1)^3 ~ (x-1)^3,  sin(pi*x) ~ -pi(x-1),  cos(2pi*x)-1 ~ -2pi^2(x-1)^2\n";
    std::cout << "       => p = 3+1+2 = 6,  q = 1 - 1/6 = 0.833333\n";
    std::cout << "  x=k (k=2,3,...): sin(pi*x) ~ +-pi(x-k),  cos(2pi*x)-1 ~ -2pi^2(x-k)^2\n";
    std::cout << "       => p = 1+2 = 3,  q = 1 - 1/3 = 0.666667\n\n";
    std::cout << "eps = 1e-" << n << "\n\n";

    std::vector<Root> roots = {
        {0.90, 1.0, 6},
        {1.90, 2.0, 3},
        {2.90, 3.0, 3},
    };

    for (const auto &r : roots) {
        std::cout << std::string(60, '-') << "\n";
        std::cout << "root x* = " << r.exact << "  (p=" << r.p
                  << ", q=1-1/p=" << std::fixed << std::setprecision(6)
                  << (1.0 - 1.0/r.p) << ")  x0=" << r.x0 << "\n\n";

        std::vector<RootStep> steps;
        auto res = newton(f, df, r.x0, eps, 100000, &steps);

        std::cout << "  iter   approximation       error           q=err/prev_err\n";
        double prev_err = -1.0;
        for (const auto &s : steps) {
            double err = std::fabs(s.approximation - r.exact);
            std::cout << "  [" << std::setw(3) << s.iteration << "]  "
                      << std::fixed << std::setprecision(12) << s.approximation
                      << "  " << std::scientific << std::setprecision(3) << err;
            if (prev_err > 0.0 && err > 0.0)
                std::cout << "  " << std::fixed << std::setprecision(6) << (err / prev_err);
            std::cout << "\n";
            prev_err = err;
        }
        std::cout << "\n  converged: " << (res.converged ? "yes" : "no")
                  << "  in " << res.iterations << " iterations\n\n";
    }

    return 0;
}
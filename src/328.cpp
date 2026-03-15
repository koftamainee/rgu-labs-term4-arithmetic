#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

#include "root_finding.hpp"

static const double PI = std::acos(-1.0);

auto f = [](double x) {
    return std::pow(x-1.0,3) * std::sin(PI*x) * (std::cos(2.0*PI*x) - 1.0);
};
auto df = [](double x) {
    double a=(x-1.0)*(x-1.0)*(x-1.0), da=3.0*(x-1.0)*(x-1.0);
    double b=std::sin(PI*x),            db=PI*std::cos(PI*x);
    double c=std::cos(2.0*PI*x)-1.0,    dc=-2.0*PI*std::sin(2.0*PI*x);
    return da*b*c + a*db*c + a*b*dc;
};

struct Root { double x0; double exact; int p; };

int main() {
    std::cout << "f(x) = (x-1)^3 * sin(pi*x) * (cos(2*pi*x) - 1)\n\n";
    std::cout << "Modified Newton: x_{n+1} = x_n - p * f(x_n)/f'(x_n)\n";
    std::cout << "sigma = -p eliminates the linear error term => quadratic convergence\n\n";

    std::vector<Root> roots = {{0.90, 1.0, 6}, {1.90, 2.0, 3}, {2.90, 3.0, 3}};

    for (int n : {3, 4, 5, 6}) {
        double eps = std::pow(10.0, -static_cast<double>(n));
        std::cout << std::string(70, '=') << "\n";
        std::cout << "eps = 1e-" << n << "\n\n";

        for (const auto &r : roots) {
            std::cout << "  root x*=" << r.exact << "  p=" << r.p
                      << "  x0=" << r.x0 << "\n";

            std::vector<RootStep> steps_n, steps_m;

            newton(f, df, r.x0, eps, 100000, &steps_n);
            newton_modified(f, df, r.x0, static_cast<double>(r.p), eps, 100000, &steps_m);

            std::cout << "    standard Newton    : " << std::setw(5) << steps_n.size() << " iter\n";
            std::cout << "    modified Newton    : " << std::setw(5) << steps_m.size() << " iter";
            std::cout << "  (speedup x"
                      << std::fixed << std::setprecision(1)
                      << static_cast<double>(steps_n.size()) / steps_m.size() << ")\n\n";
        }
    }

    std::cout << std::string(70, '=') << "\n";
    std::cout << "Modified Newton iterations (detailed, eps=1e-6):\n\n";

    double eps = 1e-6;
    for (const auto &r : roots) {
        std::cout << "  root x*=" << r.exact << "  p=" << r.p << "  sigma=-p=" << -r.p << "\n";
        std::cout << "  iter   approximation         error          q=err/prev^2\n";
        std::vector<RootStep> steps;
        newton_modified(f, df, r.x0, static_cast<double>(r.p), eps, 100000, &steps);
        double prev = -1.0, prev2 = -1.0;
        for (const auto &s : steps) {
            double err = std::fabs(s.approximation - r.exact);
            std::cout << "  [" << std::setw(2) << s.iteration << "]  "
                      << std::fixed << std::setprecision(12) << s.approximation
                      << "  " << std::scientific << std::setprecision(3) << err;
            if (prev2 > 1e-15 && err > 1e-15)
                std::cout << "  " << std::fixed << std::setprecision(4) << err/(prev*prev);
            std::cout << "\n";
            prev2 = prev; prev = err;
        }
        std::cout << "\n";
    }

    return 0;
}
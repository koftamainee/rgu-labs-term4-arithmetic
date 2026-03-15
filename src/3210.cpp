#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "newton_system.hpp"

static const double PI = std::acos(-1.0);

void print_system_result(const std::string &label,
                         const SystemResult &r,
                         const std::vector<SystemStep> &steps) {
    std::cout << label << "\n";
    for (const auto &s : steps) {
        std::cout << "  [" << std::setw(3) << s.iteration << "] ";
        for (size_t i = 0; i < s.approximation.dimension(); ++i)
            std::cout << " x" << i+1 << "=" << std::fixed << std::setprecision(10)
                      << s.approximation[i];
        std::cout << "\n";
    }
    if (r.converged) {
        std::cout << "  => converged in " << r.iterations << " iter: (";
        for (size_t i = 0; i < r.solution.dimension(); ++i) {
            if (i) std::cout << ", ";
            std::cout << std::fixed << std::setprecision(10) << r.solution[i];
        }
        std::cout << ")\n\n";
    } else {
        std::cout << "  => did NOT converge\n\n";
    }
}

int main() {
    const double eps = 1e-10;
    std::cout << "=== a) ===\n\n";
    {
        auto F = [](const Vector &v) -> Vector {
            double x1=v[0], x2=v[1], x3=v[2];
            return Vector{
                x1*x1*x1 + x1*x1*x2 - x1*x3 + 6.0,
                std::exp(x1) + std::exp(x2) - x3,
                x2*x2 - 2.0*x1*x3 - 4.0
            };
        };
        auto J = [](const Vector &v) -> Matrix {
            double x1=v[0], x2=v[1], x3=v[2];
            return Matrix({{3*x1*x1 + 2*x1*x2 - x3,  x1*x1,       -x1      },
                           {std::exp(x1),              std::exp(x2), -1.0    },
                           {-2.0*x3,                   2.0*x2,       -2.0*x1 }});
        };

        for (auto x0 : std::vector<Vector>{ {-1.4,-1.7,0.5}, {-2.0,0.2,1.0} }) {
            std::vector<SystemStep> steps;
            auto r = newton_system_inf(F, J, x0, eps, 100000, &steps);
            std::string lbl = "  x0=(" + std::to_string(x0[0]) + "," +
                              std::to_string(x0[1]) + "," + std::to_string(x0[2]) + ")";
            print_system_result(lbl, r, steps);
        }
    }

    std::cout << "=== b) ===\n\n";
    {
        auto F = [](const Vector &v) -> Vector {
            double x1=v[0], x2=v[1], x3=v[2];
            return Vector{
                6.0*x1 - 2.0*std::cos(x2*x3) - 1.0,
                9.0*x2 + std::sqrt(x1*x1 + std::sin(x3) + 1.06) + 0.9,
                60.0*x3 + 3.0*std::exp(-x1*x2) + 10.0*PI - 3.0
            };
        };
        auto J = [](const Vector &v) -> Matrix {
            double x1=v[0], x2=v[1], x3=v[2];
            double sq = std::sqrt(x1*x1 + std::sin(x3) + 1.06);
            return Matrix({
                { 6.0,                      2.0*x3*std::sin(x2*x3),  2.0*x2*std::sin(x2*x3)       },
                { x1/sq,                    9.0,                      std::cos(x3)/(2.0*sq)         },
                {-3.0*x2*std::exp(-x1*x2), -3.0*x1*std::exp(-x1*x2), 60.0                          }
            });
        };

        std::vector<SystemStep> steps;
        auto r = newton_system_inf(F, J, Vector{0.5, -0.1, -0.5}, eps, 100000, &steps);
        print_system_result("  x0=(0.5,-0.1,-0.5)", r, steps);
    }

    std::cout << "=== c) ===\n";
    std::cout << "  (eigenvalue problem: A*x=x4*x, |x|=1, A=[[4,-1,1],[-1,3,-2],[1,-2,3]])\n";
    std::cout << "  analytical eigenvalues: 1, 3, 6\n\n";
    {
        auto F = [](const Vector &v) -> Vector {
            double x1=v[0], x2=v[1], x3=v[2], x4=v[3];
            return Vector{
                4*x1 - x2 + x3 - x1*x4,
               -x1 + 3*x2 - 2*x3 - x2*x4,
                x1 - 2*x2 + 3*x3 - x3*x4,
                x1*x1 + x2*x2 + x3*x3 - 1.0
            };
        };
        auto J = [](const Vector &v) -> Matrix {
            double x1=v[0], x2=v[1], x3=v[2], x4=v[3];
            return Matrix({
                { 4-x4,  -1,    1,   -x1 },
                {-1,   3-x4,   -2,   -x2 },
                { 1,    -2,  3-x4,   -x3 },
                {2*x1, 2*x2, 2*x3,   0.0 }
            });
        };

        std::vector<std::pair<Vector,std::string>> x0s = {
            { { 0.0,   0.707, 0.707, 1.0}, "lambda~1  x0=(0, 0.707, 0.707, 1)" },
            { {-0.816,-0.408, 0.408, 3.0}, "lambda~3  x0=(-0.816,-0.408,0.408,3)" },
            { {-0.577, 0.577,-0.577, 6.0}, "lambda~6  x0=(-0.577,0.577,-0.577,6)" },
        };
        for (auto &[x0, label] : x0s) {
            std::vector<SystemStep> steps;
            auto r = newton_system_inf(F, J, x0, eps, 100000, &steps);
            print_system_result("  " + label, r, steps);
        }
    }

    return 0;
}
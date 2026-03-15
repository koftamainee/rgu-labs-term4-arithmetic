#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "matrix.hpp"
#include "newton_system.hpp"
#include "vector.hpp"

static const double PI = std::acos(-1.0);


void print_steps(const std::vector<SystemStep> &steps) {
    for (const auto &s : steps)
        std::cout << "    [" << std::setw(3) << s.iteration << "]  x="
                  << std::fixed << std::setprecision(10) << s.approximation[0]
                  << "  y=" << s.approximation[1] << "\n";
}

int main() {
    const double eps = 1e-10;

    // -----------------------------------------------------------------------
    // a) tan(xy + A) = x^2,  x^2/a2 + y^2/b2 = 1
    //    F1 = tan(xy+A) - x^2
    //    F2 = x^2/a2 + y^2/b2 - 1
    //    J = [ y/cos^2(xy+A) - 2x,   x/cos^2(xy+A)  ]
    //        [ 2x/a2,                 2y/b2           ]
    // -----------------------------------------------------------------------
    std::cout << "=== a) tan(xy+A) = x^2,  x^2/a^2 + y^2/b^2 = 1 ===\n\n";

    struct CaseA { double A, a2, b2; std::vector<std::pair<double,double>> x0s; std::string name; };
    std::vector<CaseA> cases_a = {
        { 0.2, 1.0/0.6, 0.5, { {-0.2, 0.7}, {0.2,-0.7}, {0.9, 0.5} }, "i.   A=0.2 a2=1/0.6 b2=1/2" },
        { 0.4, 1.0/0.8, 0.5, { {-0.4, 0.7}, {0.4,-0.7}, {1.0, 0.4} }, "ii.  A=0.4 a2=1/0.8 b2=1/2" },
        { 0.3, 1.0/0.2, 1.0/3, { {-1.0,-0.5}, {1.0, 0.5}, {0.3,-0.6} }, "iii. A=0.3 a2=1/0.2 b2=1/3" },
        { 0.0, 1.0/0.6, 0.5, { {0.0,-0.7}, {0.0, 0.7}, {0.7, 0.6} }, "iv.  A=0.0 a2=1/0.6 b2=1/2" },
    };

    for (const auto &c : cases_a) {
        std::cout << c.name << "\n";
        auto F = [&](const Vector &v) -> Vector {
            double x = v[0], y = v[1];
            return Vector{ std::tan(x*y + c.A) - x*x,
                           x*x/c.a2 + y*y/c.b2 - 1.0 };
        };
        auto J = [&](const Vector &v) -> Matrix {
            double x = v[0], y = v[1];
            double cs2 = std::cos(x*y + c.A); cs2 *= cs2;
            return Matrix({{ y/cs2 - 2*x,   x/cs2      },
                           { 2*x/c.a2,       2*y/c.b2   }});
        };
        for (auto [x0, y0] : c.x0s) {
            std::vector<SystemStep> steps;
            auto r = newton_system(F, J, Vector{x0, y0}, eps, 100000, &steps);
            std::cout << "  x0=(" << std::fixed << std::setprecision(2) << x0 << ","
                      << std::setw(5) << y0 << ")  ";
            if (r.converged)
                std::cout << "converged in " << std::setw(3) << r.iterations
                          << " iter  =>  (" << std::setprecision(10) << r.solution[0]
                          << ", " << r.solution[1] << ")\n";
            else
                std::cout << "did NOT converge\n";
        }
        std::cout << "\n";
    }

    // -----------------------------------------------------------------------
    // b) x1^2 + x2^2 - 2 = 0,  e^(x1-1) + x2^3 - 2 = 0
    //    J = [ 2x1,          2x2        ]
    //        [ e^(x1-1),     3*x2^2     ]
    // -----------------------------------------------------------------------
    std::cout << "=== b) x1^2 + x2^2 = 2,  e^(x1-1) + x2^3 = 2 ===\n\n";

    auto Fb = [](const Vector &v) -> Vector {
        double x = v[0], y = v[1];
        return Vector{ x*x + y*y - 2.0,
                       std::exp(x - 1.0) + y*y*y - 2.0 };
    };
    auto Jb = [](const Vector &v) -> Matrix {
        double x = v[0], y = v[1];
        return Matrix({{ 2*x,             2*y       },
                       { std::exp(x-1.0), 3*y*y     }});
    };

    for (auto [x0, y0] : std::vector<std::pair<double,double>>{{1.0,1.0},{0.5,1.2},{-1.0,1.0}}) {
        std::vector<SystemStep> steps;
        auto r = newton_system(Fb, Jb, Vector{x0, y0}, eps, 100000, &steps);
        std::cout << "  x0=(" << std::fixed << std::setprecision(2) << x0 << ","
                  << std::setw(5) << y0 << ")  ";
        if (r.converged) {
            std::cout << "converged in " << std::setw(3) << r.iterations
                      << " iter  =>  (" << std::setprecision(10) << r.solution[0]
                      << ", " << r.solution[1] << ")\n";
            print_steps(steps);
        } else {
            std::cout << "did NOT converge\n";
        }
    }

    return 0;
}
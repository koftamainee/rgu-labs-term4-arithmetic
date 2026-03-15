#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "root_finding.hpp"

static const double PI = std::acos(-1.0);

struct Root {
  double a, b;
  double x0_newton;
};

struct Equation {
  std::string label;
  std::function<double(double)> f;
  std::function<double(double)> df;
  std::vector<Root> roots;
  bool bisection_applicable;
};

void solve(const Equation& eq, int n) {
  double eps = std::pow(10.0, -static_cast<double>(n));
  std::cout << eq.label << "  (eps = 1e-" << n << ")\n";

  for (const auto& r : eq.roots) {
    if (eq.bisection_applicable) {
      std::cout << "  [bisection]  interval [" << std::fixed << std::setprecision(4)
        << r.a << ", " << r.b << "]\n";
      std::vector<RootStep> steps;
      auto res = bisection(eq.f, r.a, r.b, eps, 100000, &steps);
      for (const auto& s : steps)
        std::cout << "    iter " << std::setw(4) << s.iteration
          << "  x = " << std::setprecision(n + 3) << s.approximation << "\n";
      std::cout << "    => root = " << std::setprecision(n + 3) << res.root
        << "  (" << res.iterations << " iter)\n\n";
    }

    std::cout << "  [newton]     x0 = " << std::fixed << std::setprecision(4)
      << r.x0_newton << "\n";
    std::vector<RootStep> steps;
    try {
      auto res = newton(eq.f, eq.df, r.x0_newton, eps, 100000, &steps);
      for (const auto& s : steps)
        std::cout << "    iter " << std::setw(4) << s.iteration
          << "  x = " << std::setprecision(n + 3) << s.approximation << "\n";
      std::cout << "    => root = " << std::setprecision(n + 3) << res.root
        << "  (" << res.iterations << " iter)\n\n";
    }
    catch (const std::exception& e) {
      std::cout << "    failed: " << e.what() << "\n\n";
    }
  }
}

int main() {
  int n = 6;

  std::vector<Equation> equations = {
    {
      "a) sin(x) - 2x^2 + 0.5 = 0",
      [](double x) { return std::sin(x) - 2 * x * x + 0.5; },
      [](double x) { return std::cos(x) - 4 * x; },
      {
        {-0.35, -0.28, -0.31},
        {0.73, 0.80, 0.76}
      },
      true
    },
    {
      "b) x^n = a  (n=3, a=2)  <=>  x^3 - 2 = 0",
      [](double x) { return x * x * x - 2.0; },
      [](double x) { return 3.0 * x * x; },
      {{1.0, 2.0, 1.5}},
      true
    },
    {
      "c) sqrt(1-x^2) - e^x + 0.1 = 0",
      [](double x) { return std::sqrt(1.0 - x * x) - std::exp(x) + 0.1; },
      [](double x) { return -x / std::sqrt(1.0 - x * x) - std::exp(x); },
      {
        {-0.97, -0.95, -0.96},
        {0.08, 0.10, 0.09}
      },
      true
    },
    {
      "d) x^6 - 5x^3 - 2 = 0",
      [](double x) { return std::pow(x, 6) - 5 * x * x * x - 2; },
      [](double x) { return 6 * std::pow(x, 5) - 15 * x * x; },
      {
        {-0.75, -0.70, -0.72},
        {1.74, 1.76, 1.75}
      },
      true
    },
    {
      "e) log2(x) - 1/(1+x^2) = 0",
      [](double x) { return std::log2(x) - 1.0 / (1.0 + x * x); },
      [](double x) { return 1.0 / (x * std::log(2.0)) + 2 * x / std::pow(1 + x * x, 2); },
      {{1.28, 1.32, 1.30}},
      true
    },
    {
      "f) sin(x/2) - 1 = 0  [root has multiplicity 2 => bisection inapplicable]",
      [](double x) { return std::sin(x / 2.0) - 1.0; },
      [](double x) { return 0.5 * std::cos(x / 2.0); },
      {
        {0.0, 0.0, PI},
        {0.0, 0.0, -3 * PI},
        {0.0, 0.0, 5 * PI}
      },
      false
    },
    {
      "g) ln(x) - 1 = 0",
      [](double x) { return std::log(x) - 1.0; },
      [](double x) { return 1.0 / x; },
      {{2.0, 3.0, 2.5}},
      true
    },
  };

  for (const auto& eq : equations) {
    std::cout << std::string(70, '-') << "\n";
    solve(eq, n);
  }

  return 0;
}

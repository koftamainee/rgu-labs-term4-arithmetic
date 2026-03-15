#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <vector>

#include "simple_iteration.hpp"

double phi_a(double x, double a, double b, double c) {
  return c + a * std::sin(x) * std::sin(x) + b * std::cos(x) * std::cos(x);
}

double phi_b(double x, double a, double b, double c) {
  return c + a * std::exp(-b * x * x);
}

bool converges_a(double a, double b) { return std::fabs(a - b) < 1.0; }
bool converges_b(double a, double b) { return b > 0.0 && 2.0 * a * a * b < std::exp(1.0); }
double bound_a(double a, double b) { return std::fabs(a - b); }

double bound_b(double a, double b) {
  if (b <= 0.0) return std::numeric_limits<double>::infinity();
  return std::fabs(a) * std::sqrt(2.0 * b / std::exp(1.0));
}

void print_result(double x0, const IterationResult& r) {
  std::cout << "    x0=" << std::fixed << std::setprecision(1) << std::setw(6) << x0 << "  =>  ";
  if (r.converged)
    std::cout << "converged in " << std::setw(5) << r.iterations
      << " iter,  root = " << std::setprecision(10) << r.root << "\n";
  else
    std::cout << "did NOT converge after " << r.iterations << " iter\n";
}

int main() {
  std::cout << "=== a) phi(x) = c + a*sin^2(x) + b*cos^2(x) ===\n";
  std::cout << "    Convergence for ANY x0 if |a - b| < 1  (c is free)\n\n";

  struct CaseA {
    double a, b, c;
    std::vector<double> x0s;
  };
  const std::vector<CaseA> cases_a = {
    {0.3, 0.5, 1.0, {0.0, 2.0, -3.0, 10.0}},
    {0.8, -0.1, 0.5, {0.0, 1.0, -5.0}},
    {1.5, 0.2, 0.0, {0.0, 1.0, -5.0, 10.0}},
  };

  for (const auto& c : cases_a) {
    std::cout << std::fixed << std::setprecision(4)
      << "  a=" << c.a << "  b=" << c.b << "  c=" << c.c
      << "  |a-b|=" << bound_a(c.a, c.b)
      << "  =>  " << (converges_a(c.a, c.b) ? "CONVERGES" : "DIVERGES") << "\n";
    for (double x0 : c.x0s)
      print_result(x0, simple_iteration([&](double x) { return phi_a(x, c.a, c.b, c.c); }, x0));
    std::cout << "\n";
  }

  std::cout << "=== b) phi(x) = c + a*exp(-b*x^2) ===\n";
  std::cout << "    Convergence for ANY x0 if b > 0  AND  2*a^2*b < e = "
    << std::fixed << std::setprecision(10) << std::exp(1.0) << "\n\n";

  struct CaseB {
    double a, b, c;
    std::vector<double> x0s;
  };
  std::vector<CaseB> cases_b = {
    {0.5, 0.3, 0.0, {0.0, 1.0, -2.0, 5.0}},
    {1.0, 1.0, 0.0, {0.0, 0.5, -1.0}},
    {2.0, 1.0, 0.0, {0.0, 1.0, -3.0}},
    {0.5, -0.5, 0.0, {0.0, 1.0, 5.0}},
  };

  for (const auto& c : cases_b) {
    double q = bound_b(c.a, c.b);
    std::cout << std::fixed << std::setprecision(4)
      << "  a=" << c.a << "  b=" << c.b << "  c=" << c.c
      << "  |a|*sqrt(2b/e)=";
    if (std::isfinite(q)) std::cout << std::setprecision(4) << q;
    else std::cout << "inf";
    std::cout << "  =>  " << (converges_b(c.a, c.b) ? "CONVERGES" : "DIVERGES") << "\n";
    for (double x0 : c.x0s)
      print_result(x0, simple_iteration([&](double x) { return phi_b(x, c.a, c.b, c.c); }, x0));
    std::cout << "\n";
  }

  return 0;
}

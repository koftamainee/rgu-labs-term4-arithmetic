#include "simple_iteration.hpp"

#include <iomanip>
#include <iostream>
#include <string>



inline IterationResult eq_a_phi1(
  double x0,
  double eps = 1e-10,
  size_t max_iter = 100000) {
  auto phi = [](double x) -> double {
    return std::cbrt(1.0 - 3.0 * x * x);
  };
  return simple_iteration(phi, x0, eps, max_iter);
}

IterationResult eq_a_phi2(
  double x0 = -0.7,
  double eps = 1e-10,
  size_t max_iter = 100000) {
  auto phi = [](double x) -> double {
    return -std::sqrt(1.0 / (x + 3.0));
  };
  return simple_iteration(phi, x0, eps, max_iter);
}


IterationResult eq_b_phi1(
  double x0 = 0.1,
  double eps = 1e-10,
  size_t max_iter = 100000) {
  auto phi = [](double x) -> double {
    return x * x * x;
  };
  return simple_iteration(phi, x0, eps, max_iter);
}

inline IterationResult eq_b_phi2(
  double x0 = 0.9,
  double eps = 1e-10,
  size_t max_iter = 100000) {
  auto phi = [](double x) -> double {
    return std::pow(x, 0.75);
  };
  return simple_iteration(phi, x0, eps, max_iter);
}

IterationResult eq_c_phi1(
  double x0 = 0.5,
  double eps = 1e-10,
  size_t max_iter = 100000) {
  auto phi = [](double x) -> double {
    return (x * x + 2.0) / 3.0;
  };
  return simple_iteration(phi, x0, eps, max_iter);
}

inline IterationResult eq_c_phi2(
  double x0 = 2.5,
  double eps = 1e-10,
  size_t max_iter = 100000) {
  auto phi = [](double x) -> double {
    return 3.0 - 2.0 / x;
  };
  return simple_iteration(phi, x0, eps, max_iter);
}


inline IterationResult eq_a_phi3(
  double x0 = 0.5,
  double eps = 1e-10,
  size_t max_iter = 100000) {
  auto phi = [](double x) -> double {
    return std::sqrt((1.0 - x * x * x) / 3.0);
  };
  return simple_iteration(phi, x0, eps, max_iter);
}

static void print_result(const std::string &label, const IterationResult &r)
{
    std::cout << label << "\n";
    if (r.converged) {
        std::cout << "  Status : converged in " << r.iterations << " iteration(s)\n";
        std::cout << "  Root   : " << std::fixed << std::setprecision(10) << r.root << "\n";
    } else {
        std::cout << "  Status : did NOT converge after " << r.iterations << " iteration(s)\n";
        std::cout << "  Last x : " << std::fixed << std::setprecision(10) << r.root << "\n";
    }
    std::cout << "\n";
}

int main()
{
    constexpr double eps = 1e-10;

    std::cout << "=== Simple Iteration Method ===\n\n";

    std::cout << "--- a) x^3 + 3x^2 - 1 = 0 ----------------------------\n\n";

    print_result("phi3(x) = sqrt((1-x^3)/3),  x0 =  0.5  =>  root ~  0.532",
                 eq_a_phi3(0.5, eps));

    print_result("phi1(x) = cbrt(1 - 3x^2),   x0 = -3.0  =>  root ~ -2.879",
                 eq_a_phi1(-3.0, eps));

    print_result("phi2(x) = -sqrt(1/(x+3)),   x0 = -0.7  =>  root ~ -0.653",
                 eq_a_phi2(-0.7, eps));

    std::cout << "--- b) x^4 - x^3 = 0 ---------------------------------\n\n";

    print_result("phi1(x) = x^3,      x0 = 0.1  =>  root = 0",
                 eq_b_phi1(0.1, eps));

    print_result("phi2(x) = x^(3/4),  x0 = 0.9  =>  root = 1",
                 eq_b_phi2(0.9, eps));

    std::cout << "--- c) x^2 - 3x + 2 = 0 ------------------------------\n\n";

    print_result("phi1(x) = (x^2 + 2)/3,  x0 = 0.5  =>  root = 1",
                 eq_c_phi1(0.5, eps));

    print_result("phi2(x) = 3 - 2/x,      x0 = 2.5  =>  root = 2",
                 eq_c_phi2(2.5, eps));

    return 0;
}
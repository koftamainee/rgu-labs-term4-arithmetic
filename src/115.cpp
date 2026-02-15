#include "bigfloat.h"
#include "composite_rational_function.hpp"
#include "limit.hpp"
#include "polynomial.hpp"
#include "vector.h"
#include <iostream>

int main() {

  std::cout << "Task: Compute limits of the composite rational function\n";
  std::cout << "T(x) = f_1 (s_1(x))^k  / f_2 (s_2(x))^l\n";
  std::cout << "             \n";
  std::cout << "where f_1, f_2, s_1, s_2 are polynomials and k, l are positive "
               "integers.\n";
  std::cout << "The program calculates:\n";
  std::cout << "1) lim(x -> c) T(x) at finite points c,\n";
  std::cout << "2) lim(x -> +infinity) T(x),\n";
  std::cout << "3) lim(x -> -infinity) T(x).\n";

  std::cout << "Example 1:" << std::endl;

  Vector f1_coeffs = {-1, 1};
  Vector s1_coeffs = {0, 0, 1};
  size_t k = 2;

  Vector f2_coeffs = {1, 1};
  Vector s2_coeffs = {0, 1};
  size_t l = 1;

  Polynomial f1(f1_coeffs, 0);
  Polynomial s1(s1_coeffs, 0);
  Polynomial f2(f2_coeffs, 0);
  Polynomial s2(s2_coeffs, 0);

  CompositeRationalFunction T1(f1, s1, k, f2, s2, l);

  std::cout << "T(x) = " << T1.to_string() << std::endl;

  Limit lim1_at_1 = T1.limit_at_point(1);
  std::cout << "lim(x->1) T(x) = " << lim1_at_1.to_string() << std::endl;

  Limit lim1_at_neg1 = T1.limit_at_point(-1);
  std::cout << "lim(x->-1) T(x) = " << lim1_at_neg1.to_string() << std::endl;

  Limit lim1_plus_inf = T1.limit_at_plus_infinity();
  std::cout << "lim(x->+inf) T(x) = " << lim1_plus_inf.to_string() << std::endl;

  Limit lim1_minus_inf = T1.limit_at_minus_infinity();
  std::cout << "lim(x->-inf) T(x) = " << lim1_minus_inf.to_string()
            << std::endl;

  std::cout << std::endl;

  std::cout << "Example 2: " << std::endl;

  Vector f1_coeffs2 = {0, 1};
  Vector s1_coeffs2 = {-1, 1};
  size_t k2 = 1;

  Vector f2_coeffs2 = {0, 1};
  Vector s2_coeffs2 = {0, 0, 1};
  size_t l2 = 2;

  Polynomial f1_2(f1_coeffs2, 0);
  Polynomial s1_2(s1_coeffs2, 0);
  Polynomial f2_2(f2_coeffs2, 0);
  Polynomial s2_2(s2_coeffs2, 0);

  CompositeRationalFunction T2(f1_2, s1_2, k2, f2_2, s2_2, l2);

  std::cout << "T(x) = " << T2.to_string() << std::endl;

  Limit lim2_at_0 = T2.limit_at_point(0);
  std::cout << "lim(x->0) T(x) = " << lim2_at_0.to_string() << std::endl;

  Limit lim2_at_1 = T2.limit_at_point(1);
  std::cout << "lim(x->1) T(x) = " << lim2_at_1.to_string() << std::endl;

  Limit lim2_plus_inf = T2.limit_at_plus_infinity();
  std::cout << "lim(x->+inf) T(x) = " << lim2_plus_inf.to_string() << std::endl;

  Limit lim2_minus_inf = T2.limit_at_minus_infinity();
  std::cout << "lim(x->-inf) T(x) = " << lim2_minus_inf.to_string()
            << std::endl;

  std::cout << std::endl;

  std::cout << "Example 3" << std::endl;

  Vector f1_coeffs3 = {1};
  Vector s1_coeffs3 = {0, 1};
  size_t k3 = 1;

  Vector f2_coeffs3 = {0, 1};
  Vector s2_coeffs3 = {0, 1};
  size_t l3 = 1;

  Polynomial f1_3(f1_coeffs3, 0);
  Polynomial s1_3(s1_coeffs3, 0);
  Polynomial f2_3(f2_coeffs3, 0);
  Polynomial s2_3(s2_coeffs3, 0);

  CompositeRationalFunction T3(f1_3, s1_3, k3, f2_3, s2_3, l3);

  std::cout << "T(x) = " << T3.to_string() << std::endl;

  Limit lim3_at_0 = T3.limit_at_point(0);
  std::cout << "lim(x->0) T(x) = " << lim3_at_0.to_string() << std::endl;

  Limit lim3_plus_inf = T3.limit_at_plus_infinity();
  std::cout << "lim(x->+inf) T(x) = " << lim3_plus_inf.to_string() << std::endl;

  Limit lim3_minus_inf = T3.limit_at_minus_infinity();
  std::cout << "lim(x->-inf) T(x) = " << lim3_minus_inf.to_string()
            << std::endl;

  return 0;
}

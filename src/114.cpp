#include "bigfloat.h"
#include "limit.hpp"
#include "rational_function.hpp"
#include "vector.h"
#include <iostream>

int main(void) {

  std::cout << "Task: Compute limits of rational functions.\n";
  std::cout << "Given a rational function R(x) = P(x)/Q(x), where P and Q are "
               "polynomials,\n";
  std::cout << "the program calculates:\n";
  std::cout << "1) lim(x -> c) R(x) at a finite point c,\n";
  std::cout << "2) lim(x -> +infinity) R(x),\n";
  std::cout << "3) lim(x -> -infinity) R(x).\n";

  std::cout << "Example 1: " << std::endl;
  Vector f1 = {-1, 0, 1};
  Vector g1 = {-1, 1};
  RationalFunction R1(f1, g1);

  std::cout << R1.to_string() << std::endl;

  Limit lim1_at_1 = R1.limit_at_point(1);
  std::cout << "lim(x->1) R(x) = " << lim1_at_1.to_string() << std::endl;

  Limit lim1_plus_inf = R1.limit_at_plus_infinity();
  std::cout << "lim(x->+inf) R(x) = " << lim1_plus_inf.to_string() << std::endl;

  Limit lim1_minus_inf = R1.limit_at_minus_infinity();
  std::cout << "lim(x->-inf) R(x) = " << lim1_minus_inf.to_string()
            << std::endl;

  std::cout << std::endl;

  std::cout << "Example 2: " << std::endl;
  Vector f2 = {1};
  Vector g2 = {0, 1};
  RationalFunction R2(f2, g2);

  std::cout << R2.to_string() << std::endl;

  Limit lim2_at_0 = R2.limit_at_point(0);
  std::cout << "lim(x->0) R(x) = " << lim2_at_0.to_string() << std::endl;

  Limit lim2_plus_inf = R2.limit_at_plus_infinity();
  std::cout << "lim(x->+inf) R(x) = " << lim2_plus_inf.to_string() << std::endl;

  Limit lim2_minus_inf = R2.limit_at_minus_infinity();
  std::cout << "lim(x->-inf) R(x) = " << lim2_minus_inf.to_string()
            << std::endl;

  std::cout << std::endl;

  std::cout << "Example 3: " << std::endl;
  Vector f3 = {-4, 0, 1};
  Vector g3 = {0, -2, 1};
  RationalFunction R3(f3, g3);

  std::cout << R3.to_string() << std::endl;

  Limit lim3_at_0 = R3.limit_at_point(0);
  std::cout << "lim(x->0) R(x) = " << lim3_at_0.to_string() << std::endl;

  Limit lim3_at_2 = R3.limit_at_point(2);
  std::cout << "lim(x->2) R(x) = " << lim3_at_2.to_string() << std::endl;

  Limit lim3_plus_inf = R3.limit_at_plus_infinity();
  std::cout << "lim(x->+inf) R(x) = " << lim3_plus_inf.to_string() << std::endl;

  std::cout << std::endl;

  return 0;
}

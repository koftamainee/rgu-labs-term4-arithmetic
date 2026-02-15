#include "bigfloat.h"
#include "poly_tostring.hpp"
#include <iostream>

int main() {

  std::cout << "Task: Convert a polynomial from standard coefficient form to "
               "(x - a) form.\n";
  std::cout << "Given coefficients represent the polynomial:\n";
  std::cout << "  P(x) = c0 + c1*x + c2*x^2 + ...\n";
  std::cout << "The program converts it to the form based on (x - a) using "
               "Horner's scheme.\n\n";

  std::vector<bigfloat> coeffs = {1, 2, 3};
  bigfloat a = 5;

  size_t n = coeffs.size();
  std::vector<bigfloat> new_coeffs = coeffs;

  for (int j = n - 2; j >= 0; --j) {
    new_coeffs[j] += a * new_coeffs[j + 1];
  }

  std::cout << "Initial poly: " << poly_tostring(coeffs, 0) << std::endl;
  std::cout << "Shift value a = " << a << std::endl;
  std::cout << "Shifted: " << poly_tostring(new_coeffs, a) << std::endl;

  return 0;
}

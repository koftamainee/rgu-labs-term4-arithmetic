#include "bigfloat.h"
#include "poly_tostring.hpp"
#include <iostream>

int main() {
  std::vector<bigfloat> coeffs = {1, 2, 3};
  bigfloat a = 5;

  size_t n = coeffs.size();
  std::vector<bigfloat> new_coeffs = coeffs;

  for (int j = n - 2; j >= 0; --j) {
    new_coeffs[j] += a * new_coeffs[j + 1];
  }

  std::cout << "Initial poly: " << poly_tostring(coeffs, 0) << std::endl;
  std::cout << "In (x - a) form: " << poly_tostring(new_coeffs, a) << std::endl;

  return 0;
}

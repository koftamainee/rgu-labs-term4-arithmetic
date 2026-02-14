#include "bigfloat.h"
#include "vector.h"
#include <iostream>

int main() {
  Vector coeffs = {1, 2, 3};
  bigfloat a = 5;

  size_t n = coeffs.dimension();
  Vector new_coeffs = coeffs;

  for (int j = n - 2; j >= 0; --j) {
    new_coeffs[j] += a * new_coeffs[j + 1];
  }

  std::cout << "Coefficients in (x-a)^k form: " << new_coeffs.to_string()
            << std::endl;

  return 0;
}

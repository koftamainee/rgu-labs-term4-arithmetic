#include "bigfloat.h"
#include "vector.h"
#include <iostream>

int main() {
  Vector coeffs = {1, 2, 3};
  bigfloat a = 1;
  bigfloat B = 2;
  size_t n = coeffs.dimension();

  Vector new_coeffs = coeffs;

  bigfloat h = B - a;

  for (size_t i = 0; i < n; i++) {
    for (size_t j = n - 1; j > i; j--) {
      new_coeffs[j - 1] += h * new_coeffs[j];
    }
  }

  std::cout << "Coefficients in (x-B)^k form: " << new_coeffs.to_string()
            << std::endl;
  return 0;
}

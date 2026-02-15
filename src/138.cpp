#include "bigfloat.h"
#include "vector.h"
#include <iostream>

bigfloat evaluate_polynomial_factorial(const Vector &coeffs,
                                       const bigfloat &x) {
  size_t n = coeffs.dimension();
  if (n == 0) {
    return bigfloat(0);
  }

  bigfloat res = coeffs[0];
  bigfloat xp = 1;

  for (size_t k = 1; k < n; k++) {
    xp *= (x - (k - 1));
    res += coeffs[k] * xp;
  }

  return res;
}

int main(void) {
  Vector coeffs({1, 2, 3, 4});

  bigfloat x = bigfloat(5);

  bigfloat val_factorial = evaluate_polynomial_factorial(coeffs, x);

  std::cout << "Polynomial coefficients: ";
  for (size_t i = 0; i < coeffs.dimension(); ++i) {
    std::cout << coeffs[i] << " ";
  }
  std::cout << "\n";

  std::cout << "x = " << x << "\n";
  std::cout << "Factorial powers evaluation: " << val_factorial << "\n";

  return 0;
}

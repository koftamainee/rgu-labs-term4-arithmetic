
#include "bigfloat.h"
#include <complex>
#include <iostream>
#include <vector>

std::complex<bigfloat> horner(const std::vector<std::complex<bigfloat>> &coeffs,
                              std::complex<bigfloat> z) {
  std::complex<bigfloat> result = coeffs.back();
  for (int i = coeffs.size() - 2; i >= 0; --i) {
    result = result * z + coeffs[i];
  }
  return result;
}

std::complex<bigfloat>
scheme3(const std::vector<std::complex<bigfloat>> &coeffs,
        std::complex<bigfloat> z) {

  int n = coeffs.size() - 1;
  bigfloat x = z.real();
  bigfloat y = z.imag();

  std::complex<bigfloat> a = coeffs[n];
  std::complex<bigfloat> b = coeffs[n - 1];

  bigfloat r = 2 * x;
  bigfloat s = x * x + y * y;

  for (int j = 1; j < n; ++j) {
    std::complex<bigfloat> a_new = b + r * a;
    std::complex<bigfloat> b_new = coeffs[n - j - 1] - s * a;
    a = a_new;
    b = b_new;
  }
  return a * z + b;
}

int main(void) {
  std::vector<std::complex<bigfloat>> coeffs = {
      {1, 0}, {2, 0}, {3, 1}, {4, -2}};
  std::complex<bigfloat> z(1.0, 2.0);

  std::complex<bigfloat> val_horner = horner(coeffs, z);
  std::complex<bigfloat> val_scheme3 = scheme3(coeffs, z);

  std::cout << "Horner: " << val_horner << "\n";
  std::cout << "Scheme3: " << val_scheme3 << "\n";

  return 0;
}

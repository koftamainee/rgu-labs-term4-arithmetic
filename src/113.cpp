#include "bigfloat.h"
#include "poly_tostring.hpp"

int main() {
  std::vector<bigfloat> coeffs = {1, 2, 3};
  bigfloat a = 1;
  bigfloat B = 2;
  size_t n = coeffs.size();

  std::vector<bigfloat> new_coeffs = coeffs;

  bigfloat h = B - a;

  for (size_t i = 0; i < n; i++) {
    for (size_t j = n - 1; j > i; j--) {
      new_coeffs[j - 1] += h * new_coeffs[j];
    }
  }

  std::cout << "Poly in (x - a) form: " << poly_tostring(coeffs, a)
            << std::endl;
  std::cout << "In (x - B) form: " << poly_tostring(new_coeffs, B) << std::endl;

  return 0;
}

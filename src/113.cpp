#include "bigfloat.h"
#include "poly_tostring.hpp"

int main() {

  std::cout << "Task: Convert a polynomial P(x) from (x - a) form to (x - B) "
               "form.\n";
  std::cout
      << "If P(x) = c_0 + c_1*(x - a) + c_2*(x - a)^2 + ... + c_n*(x - a)^n,\n";
  std::cout << "compute coefficients d0, d1, ..., dn such that:\n";
  std::cout
      << "P(x) = d_0 + d_1*(x - B) + d_2*(x - B)^2 + ... + d_n*(x - B)^n.\n\n";

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

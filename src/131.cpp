#include "bigfloat.h"
#include "poly_tostring.hpp"
#include <vector>

bigfloat odd_poly_eval(bigfloat x, const std::vector<bigfloat> &coeffs) {

  std::cout << "Task: Evaluate an odd polynomial at a given point.\n";
  std::cout << "An odd polynomial has the form:\n";
  std::cout << "  P(x) = c0*x + c1*x^3 + c2*x^5 + ... + cn*x^(2n+1)\n\n";

  bigfloat t = x * x;
  int n = coeffs.size();

  bigfloat p = coeffs[n - 1];

  for (int i = n - 2; i >= 0; --i) {
    p = p * t + coeffs[i];
  }

  return x * p;
}

int main(void) {
  std::vector<bigfloat> u = {2, -3, 1};
  bigfloat x = 2.0;

  auto val = odd_poly_eval(x, u);

  std::cout << "P(x) = " << odd_poly_tostring(u, 0) << std::endl;
  std::cout << "P(" << x << ") = " << val << std::endl;
}

#include "bigfloat.h"
#include <vector>

bigfloat odd_poly_eval(bigfloat x, const std::vector<bigfloat> &coeffs) {
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

  auto val = odd_poly_eval(2.0, u);

  std::cout << "x^5 - 3x^3 + 2x at 2.0 is " << val << std::endl;
}

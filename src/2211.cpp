#include <iostream>
#include "bigmath/bigfloat.hpp"


std::pair<bigfloat, bigfloat> divide_complex_7(
  bigfloat const &a0, const bigfloat& a1,
  bigfloat const &b0, const bigfloat& b1) {

  const bigfloat denom = b0 * b0 + b1 * b1;
  bigfloat re = (a0 * b0 + a1 * b1) / denom;
  bigfloat im = (a1 * b0 - a0 * b1) / denom;

  return {re, im};
}

std::pair<bigfloat, bigfloat> divide_complex_6(
  const bigfloat &a0, const bigfloat &a1,
  const bigfloat &b0, const bigfloat &b1) {
  bigfloat re, im;


  if (b0.abs() >= b1.abs()) {
    const bigfloat r = b1 / b0;
    const bigfloat d = b0 + b1 * r;
    re = (a0 + a1 * r) / d;
    im = (a1 - a0 * r) / d;
  }
  else {
    const bigfloat r = b0 / b1;
    const bigfloat d = b1 + b0 * r;
    re = (a0 * r + a1) / d;
    im = (a1 * r - a0) / d;
  }

  return {re, im};
}

int main() {
  std::cout << "=== Complex Division ===\n\n";

  {
    const bigfloat a0(3), a1(4), b0(1), b1(2);

    std::cout << "z1 = 3 + 4i,  z2 = 1 + 2i\n";
    std::cout << "Expected: 2.2 - 0.4i\n";

    auto [re7, im7] = divide_complex_7(a0, a1, b0, b1);
    std::cout << "Algorithm A (7-op): "
      << re7.to_decimal(6) << " + " << im7.to_decimal(6) << "i\n";

    auto [re6, im6] = divide_complex_6(a0, a1, b0, b1);
    std::cout << "Algorithm B (6-op): "
      << re6.to_decimal(6) << " + " << im6.to_decimal(6) << "i\n\n";
  }

  return 0;
}

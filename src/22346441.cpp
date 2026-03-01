#include <iostream>
#include "bigfloat.h"

std::pair<bigfloat, bigfloat> complex_mul_3(
    const bigfloat& a, const bigfloat& b,
    const bigfloat& c, const bigfloat& d)
{
  const bigfloat p1 = a * c;
  const bigfloat p2 = b * d;
  const bigfloat p3 = (a + b) * (c + d);
  const bigfloat re = p1 - p2;
  const bigfloat im = p3 - p1 - p2;
  return {re, im};
}

std::pair<bigfloat, bigfloat> complex_mul_naive(
    const bigfloat& a, const bigfloat& b,
    const bigfloat& c, const bigfloat& d)
{
  return {a * c - b * d, a * d + b * c};
}

int main() {
  const bigfloat a(3), b(4), c(1), d(2);

  std::cout << "z1 = " << a << " + " << b << "i\n";
  std::cout << "z2 = " << c << " + " << d << "i\n\n";

  auto [re3, im3] = complex_mul_3(a, b, c, d);
  std::cout << "3-mul:  " << re3 << " + " << im3 << "i\n";

  auto [re4, im4] = complex_mul_naive(a, b, c, d);
  std::cout << "naive:  " << re4 << " + " << im4 << "i\n\n";

  const bigfloat a2(-2), b2(5), c2(3), d2(-1);

  std::cout << "z1 = " << a2 << " + " << b2 << "i\n";
  std::cout << "z2 = " << c2 << " + " << d2 << "i\n\n";

  auto [re3b, im3b] = complex_mul_3(a2, b2, c2, d2);
  std::cout << "3-mul:  " << re3b << " + " << im3b << "i\n";

  auto [re4b, im4b] = complex_mul_naive(a2, b2, c2, d2);
  std::cout << "naive:  " << re4b << " + " << im4b << "i\n";

  return 0;
}
#include <iostream>
#include "gf2n.hpp"

int main() {
  std::cout << "Task: Extended Euclidean algorithm in GF(2^n).\n\n";

  struct Test {
    int n;
    GF2n::elem mod;
    GF2n::elem a;
    GF2n::elem b;
  };

  Test tests[] = {
    {4, 0b10011, 0b1011, 0b0111},    // GF(2^4), x^4 + x + 1
    {4, 0b10011, 0b1101, 0b1110},    // GF(2^4), x^4 + x + 1
    {5, 0b100101, 0b10111, 0b11010}, // GF(2^5), x^5 + x^2 + 1
    {5, 0b100101, 0b11001, 0b10011}  // GF(2^5), x^5 + x^2 + 1
  };

  for (size_t i = 0; i < std::size(tests); i++) {
    GF2n gf(tests[i].n, tests[i].mod);
    std::cout << "Example " << i + 1
              << ": GF(2^" << tests[i].n
              << ") with irreducible polynomial "
              << gf.to_string(tests[i].mod) << "\n";
    std::cout << "a = " << gf.to_string(tests[i].a) << "\n";
    std::cout << "b = " << gf.to_string(tests[i].b) << "\n";

    GF2n::elem x, y;
    GF2n::elem g = gf.egcd(tests[i].a, tests[i].b, x, y);

    std::cout << "gcd(a, b)   = " << gf.to_string(g) << "\n";
    std::cout << "x           = " << gf.to_string(x) << "\n";
    std::cout << "y           = " << gf.to_string(y) << "\n";

    const GF2n::elem lhs = gf.add(gf.mul(tests[i].a, x), gf.mul(tests[i].b, y));
    std::cout << "a*x + b*y   = " << gf.to_string(lhs)
              << (lhs == g ? "  [OK]" : "  [FAIL]") << "\n\n";
  }

  return 0;
}
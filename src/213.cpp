#include <iostream>
#include "gf2n.hpp"

int main() {
  std::cout << "Task: Multiply elements in GF(2^n) with different irreducible polynomials.\n\n";

  struct Test {
    int n;
    GF2n::elem mod;
    GF2n::elem a;
    GF2n::elem b;
  };

  Test tests[] = {
    {4, 0b10011, 0b1011, 0b110},     // GF(2^4), x^4 + x + 1
    {4, 0b10011, 0b111, 0b101},      // GF(2^4), x^4 + x + 1
    {5, 0b100101, 0b10101, 0b11011}, // GF(2^5), x^5 + x^2 + 1
    {5, 0b100101, 0b11100, 0b10110}  // GF(2^5), x^5 + x^2 + 1
  };

  for (size_t i = 0; i < std::size(tests); i++) {
    GF2n gf(tests[i].n, tests[i].mod);
    std::cout << "Example " << i+1
              << ": GF(2^" << tests[i].n
              << ") with irreducible polynomial "
              << gf.to_string(tests[i].mod) << "\n";
    std::cout << "a = " << gf.to_string(tests[i].a) << "\n";
    std::cout << "b = " << gf.to_string(tests[i].b) << "\n";
    GF2n::elem c = gf.mul(tests[i].a, tests[i].b);
    std::cout << "a * b = " << gf.to_string(c) << "\n\n";
  }

  return 0;
}
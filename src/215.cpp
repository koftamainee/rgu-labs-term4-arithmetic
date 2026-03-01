#include <iostream>
#include "gf2n.hpp"

int main() {
  std::cout << "Task: Find multiplicative inverse of an element in GF(2^n).\n\n";

  struct Test {
    int n;
    GF2n::elem mod;
    GF2n::elem a;
  };

  constexpr Test tests[] = {
    {4, 0b10011, 0b1011},   // GF(2^4), x^4 + x + 1
    {4, 0b10011, 0b0111},   // GF(2^4), x^4 + x + 1
    {5, 0b100101, 0b10111}, // GF(2^5), x^5 + x^2 + 1
    {5, 0b100101, 0b11001}  // GF(2^5), x^5 + x^2 + 1
  };

  for (size_t i = 0; i < std::size(tests); i++) {
    GF2n gf(tests[i].n, tests[i].mod);
    std::cout << "Example " << i + 1
              << ": GF(2^" << tests[i].n
              << ") with irreducible polynomial "
              << gf.to_string(tests[i].mod) << "\n";
    std::cout << "a     = " << gf.to_string(tests[i].a) << "\n";

    const GF2n::elem a_inv = gf.inv(tests[i].a);
    std::cout << "a^-1  = " << gf.to_string(a_inv) << "\n";


    const GF2n::elem check = gf.mul(tests[i].a, a_inv);
    std::cout << "a * a^-1 = " << gf.to_string(check)
              << (check == 1 ? "  [OK]" : "  [FAIL]") << "\n\n";
  }

  return 0;
}
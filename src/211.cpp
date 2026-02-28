#include <bitset>

#include "gf2n.hpp"
#include "polynomial.hpp"
#include <iostream>

int main() {
  std::cout << "Task: GF(2^n) element representation in polynomial form.\n";

  GF2n gf(4, 0b10011);

  std::cout << "Field: GF(2^4) with irreducible polynomial x^4 + x + 1\n";


  GF2n::elem elems[] = {0b1011, 0b0110, 0b1111};

  for (const auto a : elems) {
    std::cout << "Element a = 0b" << std::bitset<4>(a) << "\n";

    const Polynomial p = gf.to_polynomial(a);
    std::cout << "Polynomial form: " << p.to_string() << "\n";

    const GF2n::elem b = gf.from_polynomial(p);
    std::cout << "Back to GF(2^n): 0b" << std::bitset<4>(b) << "\n";

    std::cout << std::endl;
  }

  return 0;
}
#include <iostream>
#include <cstdint>
#include <bitset>

#include "gf2n.hpp"

int main() {
  std::cout << "Task: Multiply two binary polynomials of degree <= 32.\n\n";

  const GF2n gf(4, 0b10011);

  constexpr GF2n::elem polys[][2] = {
    {0b1011, 0b110},
    {0b111, 0b101},
    {0b1001, 0b11}
  };

  for (int i = 0; i < 3; i++) {
    const GF2n::elem a = polys[i][0];
    const GF2n::elem b = polys[i][1];
    std::cout << "Example " << i+1 << ":\n";
    std::cout << "a = " << gf.to_string(a) << "\n";
    std::cout << "b = " << gf.to_string(b) << "\n";
    const GF2n::elem c = gf.mul(a, b);
    std::cout << "a * b = " << gf.to_string(c) << "\n\n";
  }

  return 0;
}
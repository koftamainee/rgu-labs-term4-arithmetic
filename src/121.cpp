#include <iostream>
#include <vector>

struct DivResult {
  std::vector<int> quotient;
  int remainder;
};

DivResult short_division(const std::vector<int> &u, int base, int divisor) {

  size_t n = u.size();

  std::vector<int> quotient(n);
  int remainder = 0;

  for (size_t i = 0; i < n; ++i) {
    int temp = remainder * base + u[i];
    quotient[i] = temp / divisor;
    remainder = temp % divisor;
  }

  return {quotient, remainder};
}

int main(void) {
  std::vector<int> u = {1, 2, 3, 4};
  int b = 10;
  int v = 7;

  auto res = short_division(u, b, v);

  for (auto d : u) {
    std::cout << d;
  }
  std::cout << " / " << v << " (in base " << b << ")" << ":\n";

  std::cout << "quotient: ";
  for (auto d : res.quotient) {
    std::cout << d;
  }
  std::cout << "\nremainder: " << res.remainder << "\n";
}

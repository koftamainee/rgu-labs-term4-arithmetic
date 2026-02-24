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

  int shift = 0;

  for (size_t i = 0; i < n; ++i) {
    int temp = remainder * base + u[i];
    auto div_res = temp / divisor;
    if (div_res != 0) {
      quotient[i - shift] = temp / divisor;
    } else {
      shift++;
    }
    remainder = temp % divisor;
  }

  return {quotient, remainder};
}

int main(void) {

  std::cout
      << "Task: Perform short division of a number represented in base b.\n";
  std::cout << "Given a number u = [u0, u1, ..., un-1] in base b and an "
               "integer divisor v,\n";
  std::cout << "compute the quotient q and remainder r such that:\n";
  std::cout << "  u = q * v + r \n";
  std::cout << "The algorithm processes digits from most significant to least "
               "significant,\n";
  std::cout << "updating the remainder at each step.\n\n";

  std::vector<int> u(250);
  for (int i = 0; i < 250; i++) {
    u[i] = i + 1;
  }

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

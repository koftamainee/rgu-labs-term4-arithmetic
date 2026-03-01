#include <iostream>
#include <vector>
#include <complex>
#include <cmath>
#include "fft.hpp"

class Polynomial {
public:
  std::vector<double> coeffs;

  explicit Polynomial(std::vector<double> c) : coeffs(std::move(c)) {}

  Polynomial operator*(const Polynomial& other) const {
    const size_t result_size = coeffs.size() + other.coeffs.size() - 1;
    size_t n = 1;
    while (n < result_size) { n <<= 1; }

    std::vector<CD> f(n), g(n);
    for (size_t i = 0; i < coeffs.size(); ++i) f[i] = coeffs[i];
    for (size_t i = 0; i < other.coeffs.size(); ++i) g[i] = other.coeffs[i];

    fft(f, false);
    fft(g, false);
    for (size_t i = 0; i < n; ++i) {f[i] *= g[i];}
    fft(f, true);

    std::vector<double> result(result_size);
    for (size_t i = 0; i < result_size; ++i)
      result[i] = std::round(f[i].real() / static_cast<double>(n));

    return Polynomial(result);
  }

  size_t dimension() const { return coeffs.size(); }
  double operator[](size_t i) const { return coeffs[i]; }
};

Polynomial set_to_poly(const std::vector<int>& S, int max_val) {
  std::vector<double> coeffs(max_val + 1, 0.0);
  for (const int x : S) coeffs[x] = 1.0;
  return Polynomial(coeffs);
}

int main() {
  constexpr int n = 5;
  const std::vector<int> A = {0, 1, 3, 7, 10};
  const std::vector<int> B = {0, 2, 4, 6, 9};
  constexpr int max_val = 10 * n;

  const Polynomial pA = set_to_poly(A, max_val);
  const Polynomial pB = set_to_poly(B, max_val);
  const Polynomial pC = pA * pB;

  std::cout << "Cartesian sum C = {x + y : x in A, y in B}\n\n";
  std::cout << "A = { ";
  for (const int x : A) std::cout << x << " ";
  std::cout << "}\n";
  std::cout << "B = { ";
  for (const int x : B) std::cout << x << " ";
  std::cout << "}\n\n";
  std::cout << "Elements of C and their multiplicities:\n";
  for (size_t i = 0; i < pC.dimension(); ++i) {
    if (pC[i] != 0.0) {
      std::cout << "  c = " << i << "  count = " << pC[i] << "\n";
    }
  }
  return 0;
}

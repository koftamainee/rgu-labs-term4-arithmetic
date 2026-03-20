#include <iostream>
#include <gmpxx.h>

#include "polynomial.hpp"
#include "poly_ring.hpp"
#include "monomial.hpp"

int main() {
  auto ring = make_ring<mpz_class>({"x_1", "x_2", "x_3"});

  TPolynomial p(ring);
  p.set(Monomial({2, 1, 0}), mpz_class(3));
  p.set(Monomial({0, 0, 1}), mpz_class(5));
  p.set(Monomial({1, 0, 0}), mpz_class(2));

  TPolynomial<mpz_class> q(ring);
  q.set(Monomial({2, 1, 0}), mpz_class(1));
  q.set(Monomial({0, 1, 0}), mpz_class(4));

  std::cout << "=== polynomials ===" << std::endl;
  std::cout << "p = " << p.to_string() << std::endl;
  std::cout << "q = " << q.to_string() << std::endl;

  std::cout << "\n=== arithmetic ===" << std::endl;
  std::cout << "p + q = " << (p + q).to_string() << std::endl;
  std::cout << "p - q = " << (p - q).to_string() << std::endl;
  std::cout << "p * q = " << (p * q).to_string() << std::endl;

  std::cout << "\n=== supp(p) ===" << std::endl;
  for (const auto& m : p.supp()) {
    std::cout << "  [";
    for (size_t i = 0; i < m.n_vars(); ++i) {
      std::cout << m[i];
      if (i + 1 < m.n_vars()) {
        std::cout << ", ";
      }
    }
    std::cout << "]" << std::endl;
  }

  std::cout << "\n=== equality ===" << std::endl;
  std::cout << "p == p: " << (p == p) << std::endl;
  std::cout << "p == q: " << (p == q) << std::endl;

  std::cout << "\n=== evaluate ===" << std::endl;
  std::vector<mpz_class> point = {mpz_class(1), mpz_class(2), mpz_class(3)};
  std::cout << "p(1, 2, 3) = " << p.evaluate(point) << std::endl;

  std::cout << "\n=== homogeneous degree ===" << std::endl;
  TPolynomial<mpz_class> h(ring);
  h.set(Monomial({2, 1, 0}), mpz_class(3));
  h.set(Monomial({0, 2, 1}), mpz_class(5));
  h.set(Monomial({1, 1, 1}), mpz_class(7));

  std::cout << "h = " << h.to_string() << std::endl;

  auto deg = h.homogeneous_degree();
  if (deg.has_value()) {
    std::cout << "h is homogeneous of degree " << deg.value() << std::endl;
  }
  else {
    std::cout << "h is not homogeneous" << std::endl;
  }

  auto deg_p = p.homogeneous_degree();
  if (deg_p.has_value()) {
    std::cout << "p is homogeneous of degree " << deg_p.value() << std::endl;
  }
  else {
    std::cout << "p is not homogeneous" << std::endl;
  }

  std::cout << "\n=== homogeneous component ===" << std::endl;
  std::cout << "p homogeneous component of degree 1: "
    << p.homogeneous_component(1).to_string() << std::endl;
  std::cout << "p homogeneous component of degree 3: "
    << p.homogeneous_component(3).to_string() << std::endl;

  std::cout << "\n=== zero polynomial ===" << std::endl;
  TPolynomial<mpz_class> zero(ring);
  std::cout << "zero = " << zero.to_string() << std::endl;
  std::cout << "zero.is_zero(): " << zero.is_zero() << std::endl;

  return 0;
}

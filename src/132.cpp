#include "bigfloat.h"
#include "polynomial.hpp"
#include "vector.h"
#include <iostream>
#include <vector>

class PolynomialOfPolynomials {
private:
  std::vector<Polynomial> coeffs_;

public:
  PolynomialOfPolynomials(const std::vector<Polynomial> &coeffs)
      : coeffs_(coeffs) {}

  std::vector<Polynomial> initialize() const { return coeffs_; }

  std::vector<Polynomial> horner_step(const bigfloat &x0) const {
    size_t n = coeffs_.size();
    if (n == 0)
      return {};

    std::vector<Polynomial> v = coeffs_;

    for (size_t k = 0; k < n; k++) {
      for (int j = static_cast<int>(n) - 1; j >= static_cast<int>(k); j--) {
        if (j < static_cast<int>(n) - 1) {
          Polynomial shifted =
              v[j + 1].change_expansion_point(v[j + 1].expansion_point() + x0);

          Vector coeffs_j = v[j].coefficients();
          Vector coeffs_jp1 = shifted.coefficients();

          size_t max_dim =
              std::max(coeffs_j.dimension(), coeffs_jp1.dimension());
          Vector new_coeffs(max_dim);

          for (size_t i = 0; i < max_dim; i++) {
            bigfloat c1 =
                (i < coeffs_j.dimension()) ? coeffs_j[i] : bigfloat(0);
            bigfloat c2 =
                (i < coeffs_jp1.dimension()) ? coeffs_jp1[i] : bigfloat(0);
            new_coeffs[i] = c1 + c2;
          }

          v[j] = Polynomial(new_coeffs, v[j].expansion_point());
        }
      }
    }

    return v;
  }

  std::vector<Polynomial> compose_direct(const bigfloat &x0) const {
    std::vector<Polynomial> result;
    result.reserve(coeffs_.size());

    for (const auto &poly : coeffs_) {
      result.push_back(
          poly.change_expansion_point(poly.expansion_point() + x0));
    }

    return result;
  }

  void print(const std::string &label) const {
    std::cout << label << ":\n";
    for (size_t i = 0; i < coeffs_.size(); i++) {
      std::cout << "  v[" << i << "] = " << coeffs_[i].to_string() << "\n";
    }
    std::cout << "\n";
  }
};

void demonstrate_horner_with_polynomials(void) {
  Polynomial v0(Vector({1, 1}), 0);
  Polynomial v1(Vector({2, 3}), 0);
  Polynomial v2(Vector({1, -1}), 0);

  std::vector<Polynomial> coeffs = {v0, v1, v2};
  PolynomialOfPolynomials poly_of_poly(coeffs);

  std::cout << "Original polynomial u(x) = v₀ + v₁·x + v₂·x²\n";
  std::cout << "where coefficients are polynomials in t:\n";
  poly_of_poly.print("Coefficients");

  bigfloat x0(1);

  auto result = poly_of_poly.horner_step(x0);

  std::cout << "Result u(x + 1) = w₀ + w₁·x + w₂·x²:\n";
  for (size_t i = 0; i < result.size(); i++) {
    std::cout << "  w[" << i << "] = " << result[i].to_string() << "\n";
  }
  std::cout << "\n";
}

int main(void) {
  try {
    demonstrate_horner_with_polynomials();
  } catch (const std::exception &e) {
    std::cerr << "Error: " << e.what() << std::endl;
    return 1;
  }

  return 0;
}

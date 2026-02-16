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
            bigfloat c1 = (i < coeffs_j.dimension()) ? coeffs_j[i] : 0;
            bigfloat c2 = (i < coeffs_jp1.dimension()) ? coeffs_jp1[i] : 0;
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

int main(void) {
  std::cout << "Task: Perform polynomial composition using Horner's method for "
               "polynomials with polynomial coefficients.\n";
  std::cout
      << "Given a polynomial u(x) = v_0 + v_1*x + v_2*x^2 + ... + v_n*x^n,\n";
  std::cout << "where each v_i is itself a polynomial in another variable\n";
  std::cout << "the program computes u(x + x0) using nested Horner "
               "steps.\n\n";

  Polynomial v0(Vector({1, 1}), 0);
  Polynomial v1(Vector({2, 3}), 0);
  Polynomial v2(Vector({1, -1}), 0);

  std::vector<Polynomial> coeffs = {v0, v1, v2};
  PolynomialOfPolynomials poly_of_poly(coeffs);

  std::cout << "Original polynomial u(x) = v_0 + v_1*x + v_2·x^2\n";
  std::cout << "where coefficients are polynomials in t:\n";
  poly_of_poly.print("Coefficients");

  bigfloat x0 = 1;

  auto result = poly_of_poly.horner_step(x0);

  std::cout << "Result u(x + " << x0 << ") = w_0 + w_1·x + w_2*x^2:\n";
  for (size_t i = 0; i < result.size(); i++) {
    std::cout << "  w[" << i << "] = " << result[i].to_string() << "\n";
  }
  std::cout << "\n";

  return 0;
}

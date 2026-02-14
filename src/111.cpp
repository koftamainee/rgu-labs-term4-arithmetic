#include <exception>
#include <iostream>
#include <vector>

#include "matrix.h"
#include "vector.h"

int main(void) {

  std::vector<std::vector<bigfloat>> basis = {
      {15, 0, 0}, {0, 42, 0}, {0, 0, 727}};

  std::vector<bigfloat> f_coeffs = {5, 7, 3};

  Matrix M(basis);
  Vector f(f_coeffs);

  std::vector<bigfloat> solution;
  try {
    solution = M.solve_gauss(f_coeffs);
  } catch (const std::exception &e) {
    std::cout << "f(x) is NOT in the span" << std::endl;
    std::cout << "reason: " << e.what() << std::endl;
    return 0;
  }

  Vector A(solution);
  std::cout << "f(x) belongs to the span. Coefficients A_l:\n";
  std::cout << A.to_string() << "\n";

  return 0;
}

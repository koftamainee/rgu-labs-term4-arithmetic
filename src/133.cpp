#include "bigfloat.h"
#include <iostream>
#include <sstream>
#include <vector>

std::string poly2d_tostring(const std::vector<std::vector<bigfloat>> &coeffs,
                            int n) {
  std::stringstream out;
  bool first = true;

  for (int i = n - 1; i >= 0; --i) {
    for (int j = n - 1; j >= 0; --j) {
      if (i + j < n && coeffs[i][j] != 0) {
        if (!first) {
          out << " + ";
        }
        out << coeffs[i][j];
        if (i > 0) {
          out << "*x";
        }
        if (i > 1) {
          out << "^" << i;
        }
        if (j > 0) {
          out << "*y";
        }
        if (j > 1) {
          out << "^" << j;
        }
        first = false;
      }
    }
  }

  if (first)
    out << "0";
  return out.str();
}

bigfloat poly2d_eval(const std::vector<std::vector<bigfloat>> &u, bigfloat x,
                     bigfloat y, int n) {
  std::vector<bigfloat> horner_x(n, 0);

  for (int j = 0; j < n; ++j) {
    bigfloat p = 0;
    for (int i = n - 1 - j; i >= 0; --i) {
      p = p * x + u[i][j];
    }
    horner_x[j] = p;
  }

  bigfloat result = 0;
  for (int j = n - 1; j >= 0; --j) {
    result = result * y + horner_x[j];
  }

  return result;
}

int main() {
  std::cout << "Task: Evaluate a two-variable polynomial P(x, y) with terms i "
               "+ j < n.\n";
  std::cout << "P(x, y) = sum_{i+j < n} u[i][j] * x^i * y^j\n\n";

  int n = 4;
  std::vector<std::vector<bigfloat>> coeffs = {
      {1, 3, 6, 10}, {2, 5, 9, 0}, {4, 8, 0, 0}, {7, 0, 0, 0}};

  bigfloat x = {3, 2}, y = 2;

  bigfloat val = poly2d_eval(coeffs, x, y, n);

  std::cout << "P(x, y) = " << poly2d_tostring(coeffs, n) << std::endl;
  std::cout << "P(" << x << ", " << y << ") = " << val << std::endl;

  std::cout << "Polynomial of two variables with terms i + j < n, n = " << n
            << "\n\n";

  return 0;
}

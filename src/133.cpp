#include "bigfloat.h"
#include <iostream>
#include <sstream>
#include <vector>

std::string poly2d_tostring(const std::vector<std::vector<bigfloat>> &coeffs) {
  std::stringstream out;
  int n = coeffs.size() - 1;
  bool first = true;

  for (int i = n; i >= 0; i--) {
    for (int j = n; j >= 0; j--) {
      if (coeffs[i][j] != 0) {
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
  return out.str();
}

bigfloat poly2d_eval(const std::vector<std::vector<bigfloat>> &coeffs,
                     bigfloat x, bigfloat y) {
  int n = coeffs.size() - 1;
  bigfloat result = 0;

  for (int j = n; j >= 0; --j) {
    bigfloat row_sum = 0.0;

    for (int i = n - j; i >= 0; --i) {
      row_sum = row_sum * x + coeffs[i][j];
    }

    result = result * y + row_sum;
  }

  return result;
}

int main(void) {

  std::vector<std::vector<bigfloat>> coeffs = {
      {1, 3, 6, 10}, {2, 5, 9, 0}, {4, 8, 0, 0}, {7, 0, 0, 0}};

  bigfloat x = {3, 2}, y = 2;

  bigfloat val = poly2d_eval(coeffs, x, y);

  std::cout << "P(x, y) = " << poly2d_tostring(coeffs) << std::endl;

  std::cout << "P(" << x << ", " << y << ") = " << val << std::endl;

  return 0;
}

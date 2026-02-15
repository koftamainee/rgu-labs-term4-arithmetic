#ifndef POLY_TOSTRING_HPP
#define POLY_TOSTRING_HPP

#include "bigfloat.h"
#include "vector.h"
#include <sstream>
#include <string>
#include <vector>

inline std::string poly_tostring(const std::vector<bigfloat> &coeffs,
                                 bigfloat a) {
  std::stringstream out;
  int n = coeffs.size() - 1;
  bool first = true;
  for (int i = n; i >= 0; --i) {
    if (coeffs[i] != 0) {
      if (!first) {
        if (coeffs[i] < 0)
          out << " - ";
        else
          out << " + ";
      } else if (coeffs[i] < 0) {
        out << "-";
      }

      if (!(coeffs[i].abs() == 1 && i > 0)) {
        out << coeffs[i].abs();
        if (i > 0) {
          out << "*";
        }
      }

      if (i > 0) {
        if (a != 0) {
          out << "(x-" << a << ")";
        } else {
          out << "x";
        }
        if (i > 1)
          out << "^" << i;
      }
      first = false;
    }
  }
  return out.str();
}

inline std::string poly_tostring(const Vector &coeffs, bigfloat a) {
  return poly_tostring(coeffs.components(), a);
}

#endif

#ifndef POLY_TOSTRING_HPP
#define POLY_TOSTRING_HPP

#include "bigfloat.h"
#include "vector.h"
#include <complex>
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

inline std::string odd_poly_tostring(const std::vector<bigfloat> &coeffs,
                                     bigfloat a) {
  std::vector<bigfloat> full_poly;
  full_poly.reserve(coeffs.size() * 2);

  for (auto const &odd_coeff : coeffs) {
    full_poly.push_back(0);
    full_poly.push_back(odd_coeff);
  }

  return poly_tostring(full_poly, a);
}

inline std::string
complex_poly_tostring(const std::vector<std::complex<bigfloat>> &coeffs,
                      std::complex<bigfloat> a) {
  std::stringstream out;
  int n = coeffs.size() - 1;
  bool first = true;

  for (int i = n; i >= 0; --i) {
    const auto &c = coeffs[i];
    if (c.real() != 0 || c.imag() != 0) {
      if (!first) {
        out << (c.real() < 0 || c.imag() < 0 ? " - " : " + ");
      } else if (c.real() < 0 || c.imag() < 0) {
        out << "-";
      }

      if (!(c.real() == 1 && c.imag() == 0 && i > 0)) {
        out << "(" << c.real();
        if (c.imag() != 0) {
          if (c.imag() > 0)
            out << "+";
          out << c.imag() << "i";
        }
        out << ")";
        if (i > 0)
          out << "*";
      }

      if (i > 0) {
        if (a != std::complex<bigfloat>(0, 0)) {
          out << "(x-(" << a.real();
          if (a.imag() != 0) {
            if (a.imag() > 0)
              out << "+";
            out << a.imag() << "i";
          }
          out << "))";
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

#endif

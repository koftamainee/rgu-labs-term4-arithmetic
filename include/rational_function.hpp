#ifndef RATIONAL_FUNCTION_HPP
#define RATIONAL_FUNCTION_HPP

#include "bigfloat.h"
#include "limit.hpp"
#include "polynomial.hpp"
#include "vector.h"
#include <stdexcept>

class RationalFunction {
private:
  Polynomial numerator_;
  Polynomial denominator_;

  static int sign(const bigfloat &x) {
    if (x > 0) {
      return 1;
    }
    if (x < 0) {
      return -1;
    }
    return 0;
  }

public:
  RationalFunction(const Polynomial &num, const Polynomial &den)
      : numerator_(num), denominator_(den) {
    if (denominator_.is_zero()) {
      throw std::invalid_argument("Denominator cannot be zero polynomial");
    }
  }

  RationalFunction(const Vector &num_coeffs, const Vector &den_coeffs,
                   const bigfloat &a = 0)
      : numerator_(num_coeffs, a), denominator_(den_coeffs, a) {
    if (denominator_.is_zero()) {
      throw std::invalid_argument("Denominator cannot be zero polynomial");
    }
  }

  bigfloat evaluate(const bigfloat &x) const {
    bigfloat den_val = denominator_.evaluate(x);

    if (den_val == bigfloat(0)) {
      throw std::domain_error("Division by zero at x = " + x.to_decimal());
    }

    bigfloat num_val = numerator_.evaluate(x);

    return num_val / den_val;
  }

  Limit limit_at_point(const bigfloat &A) const {
    Polynomial f = numerator_.change_expansion_point(A);
    Polynomial g = denominator_.change_expansion_point(A);

    size_t k_f = f.zero_order();
    size_t k_g = g.zero_order();

    size_t n_f = f.coefficients().dimension();
    size_t n_g = g.coefficients().dimension();

    if (k_f == n_f) {
      return {LimitResult::DOES_NOT_EXIST, 0};
    }
    if (k_g == n_g) {
      return {LimitResult::DOES_NOT_EXIST, 0};
    }

    if (k_f > k_g) {
      return {LimitResult::FINITE, 0};
    } else if (k_f < k_g) {
      bigfloat f_lead = f.coefficients()[k_f];
      bigfloat g_lead = g.coefficients()[k_g];

      int coeff_sign = sign(f_lead / g_lead);
      int power_diff = static_cast<int>(k_g - k_f);

      if (power_diff % 2 == 0) {
        if (coeff_sign > 0) {
          return {LimitResult::PLUS_INFINITY, 0};
        } else if (coeff_sign < 0) {
          return {LimitResult::MINUS_INFINITY, 0};
        } else {
          return {LimitResult::DOES_NOT_EXIST, 0};
        }
      } else {
        return {LimitResult::DOES_NOT_EXIST, 0};
      }
    } else {
      bigfloat f_lead = f.coefficients()[k_f];
      bigfloat g_lead = g.coefficients()[k_g];
      return {LimitResult::FINITE, f_lead / g_lead};
    }
  }

  Limit limit_at_plus_infinity() const {
    size_t deg_f = numerator_.degree();
    size_t deg_g = denominator_.degree();

    bigfloat f_lead = numerator_.coefficients()[deg_f];
    bigfloat g_lead = denominator_.coefficients()[deg_g];

    if (deg_f < deg_g) {
      return {LimitResult::FINITE, 0};
    } else if (deg_f == deg_g) {
      return {LimitResult::FINITE, f_lead / g_lead};
    } else {
      int result_sign = sign(f_lead / g_lead);
      if (result_sign > 0) {
        return {LimitResult::PLUS_INFINITY, 0};
      } else if (result_sign < 0) {
        return {LimitResult::MINUS_INFINITY, 0};
      } else {
        return {LimitResult::DOES_NOT_EXIST, 0};
      }
    }
  }

  Limit limit_at_minus_infinity() const {
    size_t deg_f = numerator_.degree();
    size_t deg_g = denominator_.degree();

    bigfloat f_lead = numerator_.coefficients()[deg_f];
    bigfloat g_lead = denominator_.coefficients()[deg_g];

    if (deg_f < deg_g) {
      return {LimitResult::FINITE, 0};
    } else if (deg_f == deg_g) {
      return {LimitResult::FINITE, f_lead / g_lead};
    } else {
      int degree_diff = static_cast<int>(deg_f - deg_g);
      int result_sign = sign(f_lead / g_lead);

      if (degree_diff % 2 == 0) {
        if (result_sign > 0) {
          return {LimitResult::PLUS_INFINITY, 0};
        } else if (result_sign < 0) {
          return {LimitResult::MINUS_INFINITY, 0};
        } else {
          return {LimitResult::DOES_NOT_EXIST, 0};
        }
      } else {
        if (result_sign > 0) {
          return {LimitResult::MINUS_INFINITY, 0};
        } else if (result_sign < 0) {
          return {LimitResult::PLUS_INFINITY, 0};
        } else {
          return {LimitResult::DOES_NOT_EXIST, 0};
        }
      }
    }
  }

  std::string to_string() const {
    return "R(x) = [" + numerator_.to_string() + "] / [" +
           denominator_.to_string() + "]";
  }
};

#endif

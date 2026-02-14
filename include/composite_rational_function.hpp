#ifndef COMPOSITE_RATIONAL_FUNCTION_HPP
#define COMPOSITE_RATIONAL_FUNCTION_HPP

#include "bigfloat.h"
#include "limit.hpp"
#include "polynomial.hpp"
#include "vector.h"
#include <stdexcept>

class CompositeRationalFunction {
private:
  Polynomial f1_;
  Polynomial s1_;
  size_t k_;
  Polynomial f2_;
  Polynomial s2_;
  size_t l_;

  static int sign(const bigfloat &x) {
    if (x > bigfloat(0)) {
      return 1;
    }
    if (x < bigfloat(0)) {
      return -1;
    }
    return 0;
  }

public:
  CompositeRationalFunction(const Polynomial &f1, const Polynomial &s1,
                            size_t k, const Polynomial &f2,
                            const Polynomial &s2, size_t l)
      : f1_(f1), s1_(s1), k_(k), f2_(f2), s2_(s2), l_(l) {
    if (f2_.is_zero() || s2_.is_zero()) {
      throw std::invalid_argument("Denominator polynomials cannot be zero");
    }
  }

  bigfloat evaluate(const bigfloat &x) const {
    bigfloat s1_val = s1_.evaluate(x);
    bigfloat s2_val = s2_.evaluate(x);

    bigfloat f1_val = f1_.evaluate(s1_val);
    bigfloat f2_val = f2_.evaluate(s2_val);

    if (f2_val == bigfloat(0)) {
      throw std::domain_error("Division by zero in f2(s2(x))");
    }

    bigfloat num = pow(f1_val, k_);
    bigfloat den = pow(f2_val, l_);

    if (den == bigfloat(0)) {
      throw std::domain_error("Division by zero");
    }

    return num / den;
  }

  Limit limit_at_point(const bigfloat &A) const {
    bigfloat s1_limit = s1_.evaluate(A);
    bigfloat s2_limit = s2_.evaluate(A);

    Polynomial f1_at_s1 = f1_.change_expansion_point(s1_limit);
    Polynomial f2_at_s2 = f2_.change_expansion_point(s2_limit);

    size_t k_f1 = f1_at_s1.zero_order();
    size_t k_f2 = f2_at_s2.zero_order();

    size_t n_f1 = f1_at_s1.coefficients().dimension();
    size_t n_f2 = f2_at_s2.coefficients().dimension();

    bool f1_zero = (k_f1 == n_f1);
    bool f2_zero = (k_f2 == n_f2);

    if (f1_zero && f2_zero && k_ > 0 && l_ > 0) {
      return {LimitResult::DOES_NOT_EXIST, 0};
    }
    if (f1_zero && k_ > 0) {
      return {LimitResult::FINITE, 0};
    }
    if (f2_zero && l_ > 0) {
      return {LimitResult::DOES_NOT_EXIST, 0};
    }

    bigfloat f1_limit_val =
        (k_f1 < n_f1) ? f1_at_s1.coefficients()[k_f1] : bigfloat(0);
    bigfloat f2_limit_val =
        (k_f2 < n_f2) ? f2_at_s2.coefficients()[k_f2] : bigfloat(1);

    if (k_f1 == 0 && k_f2 == 0) {
      if (f2_limit_val == bigfloat(0)) {
        return {LimitResult::DOES_NOT_EXIST, 0};
      }

      bigfloat num = pow(f1_limit_val, k_);
      bigfloat den = pow(f2_limit_val, l_);

      return {LimitResult::FINITE, num / den};
    }

    int effective_zero_order = static_cast<int>(k_f1) * static_cast<int>(k_) -
                               static_cast<int>(k_f2) * static_cast<int>(l_);

    if (effective_zero_order > 0) {
      return {LimitResult::FINITE, 0};
    } else if (effective_zero_order < 0) {
      if (f1_limit_val == 0 || f2_limit_val == 0) {
        return {LimitResult::DOES_NOT_EXIST, 0};
      }

      bigfloat ratio = pow(f1_limit_val, k_) / pow(f2_limit_val, l_);
      int ratio_sign = sign(ratio);

      if ((-effective_zero_order) % 2 == 0) {
        if (ratio_sign > 0) {
          return {LimitResult::PLUS_INFINITY, 0};
        } else if (ratio_sign < 0) {
          return {LimitResult::MINUS_INFINITY, 0};
        } else {
          return {LimitResult::DOES_NOT_EXIST, 0};
        }
      } else {
        return {LimitResult::DOES_NOT_EXIST, 0};
      }
    } else {
      if (f1_limit_val == bigfloat(0) || f2_limit_val == 0) {
        return {LimitResult::DOES_NOT_EXIST, 0};
      }

      bigfloat num = pow(f1_limit_val, k_);
      bigfloat den = pow(f2_limit_val, l_);

      return {LimitResult::FINITE, num / den};
    }
  }

  Limit limit_at_plus_infinity() const {
    size_t deg_s1 = s1_.degree();
    size_t deg_f1 = f1_.degree();
    size_t deg_s2 = s2_.degree();
    size_t deg_f2 = f2_.degree();

    bigfloat s1_lead = s1_.coefficients()[deg_s1];
    bigfloat f1_lead = f1_.coefficients()[deg_f1];
    bigfloat s2_lead = s2_.coefficients()[deg_s2];
    bigfloat f2_lead = f2_.coefficients()[deg_f2];

    int deg_num = static_cast<int>(deg_f1 * deg_s1 * k_);
    int deg_den = static_cast<int>(deg_f2 * deg_s2 * l_);

    bigfloat lead_num = pow(f1_lead * pow(s1_lead, deg_f1), k_);
    bigfloat lead_den = pow(f2_lead * pow(s2_lead, deg_f2), l_);

    if (deg_num < deg_den) {
      return {LimitResult::FINITE, 0};
    } else if (deg_num == deg_den) {
      return {LimitResult::FINITE, lead_num / lead_den};
    } else {
      int result_sign = sign(lead_num / lead_den);
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
    size_t deg_s1 = s1_.degree();
    size_t deg_f1 = f1_.degree();
    size_t deg_s2 = s2_.degree();
    size_t deg_f2 = f2_.degree();

    bigfloat s1_lead = s1_.coefficients()[deg_s1];
    bigfloat f1_lead = f1_.coefficients()[deg_f1];
    bigfloat s2_lead = s2_.coefficients()[deg_s2];
    bigfloat f2_lead = f2_.coefficients()[deg_f2];

    int deg_num = static_cast<int>(deg_f1 * deg_s1 * k_);
    int deg_den = static_cast<int>(deg_f2 * deg_s2 * l_);

    int parity_num = (deg_f1 * deg_s1 * k_) % 2;
    int parity_den = (deg_f2 * deg_s2 * l_) % 2;

    bigfloat lead_num = pow(f1_lead * pow(s1_lead, deg_f1), k_);
    bigfloat lead_den = pow(f2_lead * pow(s2_lead, deg_f2), l_);

    if (parity_num == 1) {
      lead_num = -lead_num;
    }
    if (parity_den == 1) {
      lead_den = -lead_den;
    }

    if (deg_num < deg_den) {
      return {LimitResult::FINITE, 0};
    } else if (deg_num == deg_den) {
      return {LimitResult::FINITE, lead_num / lead_den};
    } else {
      int result_sign = sign(lead_num / lead_den);
      if (result_sign > 0) {
        return {LimitResult::PLUS_INFINITY, 0};
      } else if (result_sign < 0) {
        return {LimitResult::MINUS_INFINITY, 0};
      } else {
        return {LimitResult::DOES_NOT_EXIST, 0};
      }
    }
  }

  std::string to_string() const {
    return "T(x) = [f1(s1(x))]^" + std::to_string(k_) + " / [f2(s2(x))]^" +
           std::to_string(l_);
  }
};

#endif

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
    bigfloat num_val = s1_.evaluate(x);
    for (size_t i = 0; i < k_; ++i) {
      num_val = f1_.evaluate(num_val);
    }

    bigfloat den_val = s2_.evaluate(x);
    for (size_t i = 0; i < l_; ++i) {
      den_val = f2_.evaluate(den_val);
    }

    if (den_val == bigfloat(0)) {
      throw std::domain_error("Division by zero");
    }

    return num_val / den_val;
  }

  Limit limit_at_point(const bigfloat &A) const {
    bigfloat s1_val;
    try {
      s1_val = s1_.evaluate(A);
    } catch (...) {
      return {LimitResult::DOES_NOT_EXIST, 0};
    }

    bigfloat num_limit = s1_val;
    for (size_t i = 0; i < k_; ++i) {
      Polynomial f1_expanded = f1_.change_expansion_point(num_limit);
      size_t k_f1 = f1_expanded.zero_order();
      size_t n_f1 = f1_expanded.coefficients().dimension();

      if (k_f1 == n_f1) {
        num_limit = bigfloat(0);
      } else if (k_f1 == 0) {
        num_limit = f1_expanded.coefficients()[0];
      } else {
        try {
          num_limit = f1_.evaluate(num_limit);
        } catch (...) {
          return {LimitResult::DOES_NOT_EXIST, 0};
        }
      }
    }

    bigfloat s2_val;
    try {
      s2_val = s2_.evaluate(A);
    } catch (...) {
      return {LimitResult::DOES_NOT_EXIST, 0};
    }

    bigfloat den_limit = s2_val;
    for (size_t i = 0; i < l_; ++i) {
      Polynomial f2_expanded = f2_.change_expansion_point(den_limit);
      size_t k_f2 = f2_expanded.zero_order();
      size_t n_f2 = f2_expanded.coefficients().dimension();

      if (k_f2 == n_f2) {
        den_limit = bigfloat(0);
      } else if (k_f2 == 0) {
        den_limit = f2_expanded.coefficients()[0];
      } else {
        try {
          den_limit = f2_.evaluate(den_limit);
        } catch (...) {
          return {LimitResult::DOES_NOT_EXIST, 0};
        }
      }
    }

    if (den_limit == bigfloat(0)) {
      if (num_limit == bigfloat(0)) {
        return {LimitResult::DOES_NOT_EXIST, 0};
      }
      bigfloat epsilon = bigfloat::DEFAULT_EPS;
      bigfloat left_val, right_val;

      try {
        left_val = evaluate(A - epsilon);
        right_val = evaluate(A + epsilon);

        int left_sign = sign(left_val);
        int right_sign = sign(right_val);

        if (left_sign == right_sign) {
          if (left_sign > 0) {
            return {LimitResult::PLUS_INFINITY, 0};
          } else if (left_sign < 0) {
            return {LimitResult::MINUS_INFINITY, 0};
          }
        }
      } catch (...) {
      }

      return {LimitResult::DOES_NOT_EXIST, 0};
    }

    return {LimitResult::FINITE, num_limit / den_limit};
  }

  Limit limit_at_plus_infinity() const {

    int current_deg_num = static_cast<int>(s1_.degree());
    bigfloat current_lead_num = s1_.coefficients()[s1_.degree()];

    for (size_t i = 0; i < k_; ++i) {
      int deg_f1 = static_cast<int>(f1_.degree());
      bigfloat lead_f1 = f1_.coefficients()[f1_.degree()];

      if (current_deg_num == 0) {
        current_lead_num = f1_.evaluate(current_lead_num);
        current_deg_num = 0;
      } else {
        int new_deg = deg_f1 * current_deg_num;
        bigfloat new_lead = lead_f1 * pow(current_lead_num, deg_f1);

        current_deg_num = new_deg;
        current_lead_num = new_lead;
      }
    }

    int current_deg_den = static_cast<int>(s2_.degree());
    bigfloat current_lead_den = s2_.coefficients()[s2_.degree()];

    for (size_t i = 0; i < l_; ++i) {
      int deg_f2 = static_cast<int>(f2_.degree());
      bigfloat lead_f2 = f2_.coefficients()[f2_.degree()];

      if (current_deg_den == 0) {
        current_lead_den = f2_.evaluate(current_lead_den);
        current_deg_den = 0;
      } else {
        int new_deg = deg_f2 * current_deg_den;
        bigfloat new_lead = lead_f2 * pow(current_lead_den, deg_f2);

        current_deg_den = new_deg;
        current_lead_den = new_lead;
      }
    }

    if (current_deg_num < current_deg_den) {
      return {LimitResult::FINITE, 0};
    } else if (current_deg_num == current_deg_den) {
      if (current_lead_den == bigfloat(0)) {
        return {LimitResult::DOES_NOT_EXIST, 0};
      }
      return {LimitResult::FINITE, current_lead_num / current_lead_den};
    } else {
      if (current_lead_den == bigfloat(0) || current_lead_num == bigfloat(0)) {
        return {LimitResult::DOES_NOT_EXIST, 0};
      }
      int result_sign = sign(current_lead_num / current_lead_den);
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

    int current_deg_num = static_cast<int>(s1_.degree());
    bigfloat current_lead_num = s1_.coefficients()[s1_.degree()];

    if (current_deg_num % 2 == 1) {
      current_lead_num = -current_lead_num;
    }

    for (size_t i = 0; i < k_; ++i) {
      int deg_f1 = static_cast<int>(f1_.degree());
      bigfloat lead_f1 = f1_.coefficients()[f1_.degree()];

      if (current_deg_num == 0) {
        current_lead_num = f1_.evaluate(current_lead_num);
        current_deg_num = 0;
      } else {
        int new_deg = deg_f1 * current_deg_num;

        bigfloat new_lead = lead_f1 * pow(current_lead_num, deg_f1);

        current_deg_num = new_deg;
        current_lead_num = new_lead;
      }
    }

    int current_deg_den = static_cast<int>(s2_.degree());
    bigfloat current_lead_den = s2_.coefficients()[s2_.degree()];

    if (current_deg_den % 2 == 1) {
      current_lead_den = -current_lead_den;
    }

    for (size_t i = 0; i < l_; ++i) {
      int deg_f2 = static_cast<int>(f2_.degree());
      bigfloat lead_f2 = f2_.coefficients()[f2_.degree()];

      if (current_deg_den == 0) {
        current_lead_den = f2_.evaluate(current_lead_den);
        current_deg_den = 0;
      } else {
        int new_deg = deg_f2 * current_deg_den;
        bigfloat new_lead = lead_f2 * pow(current_lead_den, deg_f2);

        current_deg_den = new_deg;
        current_lead_den = new_lead;
      }
    }

    if (current_deg_num < current_deg_den) {
      return {LimitResult::FINITE, 0};
    } else if (current_deg_num == current_deg_den) {
      if (current_lead_den == bigfloat(0)) {
        return {LimitResult::DOES_NOT_EXIST, 0};
      }
      return {LimitResult::FINITE, current_lead_num / current_lead_den};
    } else {
      if (current_lead_den == bigfloat(0) || current_lead_num == bigfloat(0)) {
        return {LimitResult::DOES_NOT_EXIST, 0};
      }
      int result_sign = sign(current_lead_num / current_lead_den);
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
    std::stringstream out;

    std::string num_str;
    if (k_ > 1) {
      num_str = "(" + f1_.to_string() + ")^" + std::to_string(k_) + "{" +
                s1_.to_string() + "}";
    } else {
      num_str = "(" + f1_.to_string() + "){" + s1_.to_string() + "}";
    }

    std::string den_str;
    if (l_ > 1) {
      den_str = "(" + f2_.to_string() + ")^" + std::to_string(l_) + "{" +
                s2_.to_string() + "}";
    } else {
      den_str = "(" + f2_.to_string() + "){" + s2_.to_string() + "}";
    }

    out << num_str << " / " << den_str;

    return out.str();
  }
};

#endif

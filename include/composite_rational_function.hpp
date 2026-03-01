#ifndef COMPOSITE_RATIONAL_FUNCTION_HPP
#define COMPOSITE_RATIONAL_FUNCTION_HPP

#include "bigfloat.h"
#include "limit.hpp"
#include "polynomial.hpp"
#include "vector.h"
#include <stdexcept>
#include <utility>

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

  static bigfloat power(const bigfloat &base, size_t exp) {
    bigfloat result(1);
    for (size_t i = 0; i < exp; ++i) {
      result = result * base;
    }
    return result;
  }

public:
  CompositeRationalFunction(Polynomial f1, Polynomial s1,
                            size_t k, Polynomial f2,
                            Polynomial s2, size_t l)
    : f1_(std::move(f1)), s1_(std::move(s1)), k_(k), f2_(std::move(f2)), s2_(std::move(s2)), l_(l) {
    if (f2_.is_zero() || s2_.is_zero()) {
      throw std::invalid_argument("Denominator polynomials cannot be zero");
    }
  }

  bigfloat evaluate(const bigfloat &x) const {
    bigfloat num_val = f1_.evaluate(s1_.evaluate(x));
    num_val = power(num_val, k_);

    bigfloat den_val = f2_.evaluate(s2_.evaluate(x));
    den_val = power(den_val, l_);

    if (den_val == bigfloat(0)) {
      throw std::domain_error("Division by zero");
    }

    return num_val / den_val;
  }

  Limit limit_at_point(const bigfloat &A) const {
    bigfloat s1_val = s1_.evaluate(A);
    bigfloat num_limit = f1_.evaluate(s1_val);
    num_limit = power(num_limit, k_);

    bigfloat s2_val = s2_.evaluate(A);
    bigfloat den_limit = f2_.evaluate(s2_val);
    den_limit = power(den_limit, l_);

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
    const int deg_s1 = static_cast<int>(s1_.degree());
    const bigfloat lead_s1 = s1_.coefficients()[s1_.degree()];

    const int deg_f1 = static_cast<int>(f1_.degree());
    const bigfloat lead_f1 = f1_.coefficients()[f1_.degree()];

    int current_deg_num;
    bigfloat current_lead_num;

    if (deg_s1 == 0) {
      current_lead_num = f1_.evaluate(lead_s1);
      current_deg_num = 0;
    } else {
      current_deg_num = deg_f1 * deg_s1;
      current_lead_num = lead_f1 * pow(lead_s1, deg_f1);
    }

    current_deg_num = current_deg_num * static_cast<int>(k_);
    current_lead_num = pow(current_lead_num, static_cast<int>(k_));

    const int deg_s2 = static_cast<int>(s2_.degree());
    const bigfloat lead_s2 = s2_.coefficients()[s2_.degree()];

    const int deg_f2 = static_cast<int>(f2_.degree());
    const bigfloat lead_f2 = f2_.coefficients()[f2_.degree()];

    int current_deg_den;
    bigfloat current_lead_den;

    if (deg_s2 == 0) {
      current_lead_den = f2_.evaluate(lead_s2);
      current_deg_den = 0;
    } else {
      current_deg_den = deg_f2 * deg_s2;
      current_lead_den = lead_f2 * pow(lead_s2, deg_f2);
    }

    current_deg_den = current_deg_den * static_cast<int>(l_);
    current_lead_den = pow(current_lead_den, static_cast<int>(l_));

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
      const int result_sign = sign(current_lead_num / current_lead_den);
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
    const int deg_s1 = static_cast<int>(s1_.degree());
    bigfloat lead_s1 = s1_.coefficients()[s1_.degree()];

    if (deg_s1 % 2 == 1) {
      lead_s1 = -lead_s1;
    }

    const int deg_f1 = static_cast<int>(f1_.degree());
    const bigfloat lead_f1 = f1_.coefficients()[f1_.degree()];

    int current_deg_num;
    bigfloat current_lead_num;

    if (deg_s1 == 0) {
      current_lead_num = f1_.evaluate(lead_s1);
      current_deg_num = 0;
    } else {
      current_deg_num = deg_f1 * deg_s1;
      current_lead_num = lead_f1 * pow(lead_s1, deg_f1);
    }

    current_deg_num = current_deg_num * static_cast<int>(k_);
    current_lead_num = pow(current_lead_num, static_cast<int>(k_));

    const int deg_s2 = static_cast<int>(s2_.degree());
    bigfloat lead_s2 = s2_.coefficients()[s2_.degree()];

    if (deg_s2 % 2 == 1) {
      lead_s2 = -lead_s2;
    }

    int deg_f2 = static_cast<int>(f2_.degree());
    bigfloat lead_f2 = f2_.coefficients()[f2_.degree()];

    int current_deg_den;
    bigfloat current_lead_den;

    if (deg_s2 == 0) {
      current_lead_den = f2_.evaluate(lead_s2);
      current_deg_den = 0;
    } else {
      current_deg_den = deg_f2 * deg_s2;
      current_lead_den = lead_f2 * pow(lead_s2, deg_f2);
    }

    current_deg_den = current_deg_den * static_cast<int>(l_);
    current_lead_den = pow(current_lead_den, static_cast<int>(l_));

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
      const int result_sign = sign(current_lead_num / current_lead_den);
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
#pragma once

#include <complex>
#include <format>
#include <sstream>
#include <string>

#include <gmpxx.h>

inline std::string to_latex(int v) { return std::to_string(v); }
inline std::string to_latex(unsigned int v) { return std::to_string(v); }
inline std::string to_latex(long v) { return std::to_string(v); }
inline std::string to_latex(unsigned long v) { return std::to_string(v); }
inline std::string to_latex(long long v) { return std::to_string(v); }
inline std::string to_latex(unsigned long long v) { return std::to_string(v); }
inline std::string to_latex(float v) { return std::format("{:.3}", v); }
inline std::string to_latex(double v) { return std::format("{:.3}", v); }
inline std::string to_latex(long double v) { return std::format("{:.3}", v); }
inline std::string to_latex(bool v) {
  return v ? "\\text{true}" : "\\text{false}";
}

inline std::string to_latex(const mpz_class& v) {
  return v.get_str();
}

inline std::string to_latex(const mpq_class& v) {
  mpz_class num = v.get_num();
  const mpz_class& den = v.get_den();
  if (den == 1) {
    return num.get_str();
  }
  bool negative = num < 0;
  if (negative) {
    num = -num;
  }
  std::string frac = "\\frac{" + num.get_str() + "}{" + den.get_str() + "}";
  return negative ? "-" + frac : frac;
}

inline std::string to_latex(const mpf_class& v) {
  std::ostringstream oss;
  oss << v;
  return oss.str();
}

template <typename T>
std::string to_latex(const std::complex<T>& v) {
  if (v.imag() == T{}) {
    return to_latex(v.real());
  }

  auto imag_abs = v.imag() < T{} ? -v.imag() : v.imag();
  bool imag_negative = v.imag() < T{};
  std::string imag_str = (imag_abs == T{1}) ? "" : to_latex(imag_abs);

  if (v.real() == T{}) {
    return (imag_negative ? "-" : "") + imag_str + "i";
  }
  return to_latex(v.real())
    + (imag_negative ? " - " : " + ")
    + imag_str + "i";
}

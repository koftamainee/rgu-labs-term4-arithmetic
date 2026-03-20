#pragma once

#include <complex>
#include <format>
#include <sstream>
#include <gmpxx.h>

template <>
struct std::formatter<mpz_class> {
  constexpr auto parse(auto& ctx) { return ctx.begin(); }
  auto format(const mpz_class& v, auto& ctx) const {
    return std::format_to(ctx.out(), "{}", v.get_str());
  }
};

template <>
struct std::formatter<mpq_class> {
  constexpr auto parse(auto& ctx) { return ctx.begin(); }
  auto format(const mpq_class& v, auto& ctx) const {
    return std::format_to(ctx.out(), "{}", v.get_str());
  }
};

template <>
struct std::formatter<mpf_class> {
  constexpr auto parse(auto& ctx) { return ctx.begin(); }
  auto format(const mpf_class& v, auto& ctx) const {
    std::ostringstream oss;
    oss << v;
    return std::format_to(ctx.out(), "{}", oss.str());
  }
};

template <typename T>
struct std::formatter<std::complex<T>> {
  constexpr auto parse(auto& ctx) { return ctx.begin(); }
  auto format(const std::complex<T>& v, auto& ctx) const {
    return std::format_to(ctx.out(), "({} + {}i)", v.real(), v.imag());
  }
};
#pragma once

#include <string>

enum class LimitResult {
  Finite,
  PlusInfinity,
  MinusInfinity,
  DoesNotExist
};

template <typename T>
class Limit{
public:
  using value_type = T;

  LimitResult result;
  T value;

  explicit Limit(LimitResult result, T value = T{})
    : result(result), value(std::move(value)) {}

  std::string to_string() const {
    switch (result) {
    case LimitResult::Finite: return value.to_decimal();
    case LimitResult::PlusInfinity: return "+inf";
    case LimitResult::MinusInfinity: return "-inf";
    case LimitResult::DoesNotExist: return "does not exist";
    }
    return "unknown";
  }

  std::string to_latex() const {
    switch (result) {
    case LimitResult::Finite: return value.to_decimal();
    case LimitResult::PlusInfinity: return "+\\infty";
    case LimitResult::MinusInfinity: return "-\\infty";
    case LimitResult::DoesNotExist: return "\\nexists";
    }
    return "\\text{unknown}";
  }
};

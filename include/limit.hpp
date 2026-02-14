#ifndef LIMIT_HPP
#define LIMIT_HPP

#include "bigfloat.h"
#include <string>

enum class LimitResult {
  FINITE,
  PLUS_INFINITY,
  MINUS_INFINITY,
  DOES_NOT_EXIST
};

struct Limit {
  LimitResult type;
  bigfloat value;

  std::string to_string() const {
    switch (type) {
    case LimitResult::FINITE:
      return value.to_decimal();
    case LimitResult::PLUS_INFINITY:
      return "+inf";
    case LimitResult::MINUS_INFINITY:
      return "-inf";
    case LimitResult::DOES_NOT_EXIST:
      return "does not exist";
    }
    return "unknown";
  }
};

#endif

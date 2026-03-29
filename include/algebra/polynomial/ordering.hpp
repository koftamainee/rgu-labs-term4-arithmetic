#pragma once

#include <stdexcept>

#include "monomial_order.hpp"
#include "polynomial.hpp"

namespace order {
  template <typename Order, typename T>
  Monomial lm(const Polynomial<T>& f) {
    if (f.is_zero()) {
      throw std::domain_error("lm: polynomial is zero");
    }
    std::optional<Monomial> best;
    for (const auto& elem : f.terms()) {
      if (Monomial m(elem.key); !best.has_value() || OrderTraits<Order>::compare(*best, m) < 0) {
        best = std::move(m);
      }
    }
    return *best;
  }

  template <typename Order, typename T>
  Monomial multideg(const Polynomial<T>& f) {
    return lm<Order>(f);
  }

  template <typename Order, typename T>
  T lc(const Polynomial<T>& f) {
    if (f.is_zero()) {
      throw std::domain_error("lc: polynomial is zero");
    }
    Monomial m = lm<Order>(f);
    return *f.get(m);
  }

  template <typename Order, typename T>
  Polynomial<T> lt(const Polynomial<T>& f) {
    if (f.is_zero()) {
      throw std::domain_error("lt: polynomial is zero");
    }
    Monomial m = lm<Order>(f);
    T c = *f.get(m);
    Polynomial<T> result(f.ring());
    result.set(m, c);
    return result;
  }
} // namespace order

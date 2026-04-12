#pragma once

#include "algebra/polynomial/division.hpp"
#include "algebra/polynomial/monomial.hpp"
#include "algebra/polynomial/monomial_order.hpp"
#include "algebra/polynomial/ordering.hpp"
#include "algebra/polynomial/polynomial.hpp"

#include <vector>

template <typename Order, typename T>
Monomial lcm_monomial(const Monomial& a, const Monomial& b) {
  if (a.n_vars() != b.n_vars()) {
    throw std::invalid_argument("lcm_monomial: n_vars mismatch");
  }
  Monomial::container exp(a.n_vars());
  for (Monomial::size_type i = 0; i < a.n_vars(); ++i) {
    exp[i] = std::max(a[i], b[i]);
  }
  return Monomial(exp);
}

template <typename Order, typename T>
Polynomial<T> s_polynomial(const Polynomial<T>& f, const Polynomial<T>& g) {
  auto ring = f.ring();

  const Monomial lm_f = order::lm<Order>(f);
  const Monomial lm_g = order::lm<Order>(g);
  const T lc_f = order::lc<Order>(f);
  const T lc_g = order::lc<Order>(g);

  const Monomial gamma = lcm_monomial<Order, T>(lm_f, lm_g);

  Monomial::container exp_f(ring->n_vars());
  Monomial::container exp_g(ring->n_vars());
  for (Monomial::size_type i = 0; i < ring->n_vars(); ++i) {
    exp_f[i] = gamma[i] - lm_f[i];
    exp_g[i] = gamma[i] - lm_g[i];
  }

  Polynomial<T> term_f(ring);
  term_f.set(Monomial(exp_f), T(1) / lc_f);

  Polynomial<T> term_g(ring);
  term_g.set(Monomial(exp_g), T(1) / lc_g);

  return term_f * f - term_g * g;
}

template <typename Order, typename T>
bool is_groebner_basis(const std::vector<Polynomial<T>>& basis) {
  for (std::size_t i = 0; i < basis.size(); ++i) {
    for (std::size_t j = i + 1; j < basis.size(); ++j) {
      const Polynomial<T> s = s_polynomial<Order>(basis[i], basis[j]);
      if (s.is_zero()) continue;
      const auto res = order::divide<Order>(s, basis);
      if (!res.remainder.is_zero()) return false;
    }
  }
  return true;
}

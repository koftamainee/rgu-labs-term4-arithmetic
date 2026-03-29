#pragma once

#include <vector>

#include "ordering.hpp"
#include "monomial_order.hpp"
#include "polynomial.hpp"

namespace order {

  template <typename T>
  struct DivisionResult {
    std::vector<Polynomial<T>> quotients;
    Polynomial<T> remainder;

    explicit DivisionResult(typename Polynomial<T>::ring_ptr ring, std::size_t s)
        : remainder(ring) {
      quotients.reserve(s);
      for (std::size_t i = 0; i < s; ++i) {
        quotients.emplace_back(ring);
      }
    }
  };

  template <typename Order, typename T>
  DivisionResult<T> divide(const Polynomial<T>& f,
                            const std::vector<Polynomial<T>>& divisors) {
    const auto s   = divisors.size();
    auto       ring = f.ring();

    DivisionResult<T> res(ring, s);

    Polynomial<T> p = f;

    while (!p.is_zero()) {
      Monomial lm_p = lm<Order>(p);
      T        lc_p = lc<Order>(p);

      bool division_occurred = false;

      for (std::size_t i = 0; i < s; ++i) {
        const Monomial lm_fi = lm<Order>(divisors[i]);
        const T        lc_fi = lc<Order>(divisors[i]);

        if (!lm_fi.divides(lm_p)) {
          continue;
        }

        Monomial::container quot_exp(ring->n_vars());
        for (std::size_t j = 0; j < ring->n_vars(); ++j) {
          quot_exp[j] = lm_p[j] - lm_fi[j];
        }
        Monomial quot_mono(quot_exp);
        T        quot_coeff = lc_p / lc_fi;

        Polynomial<T> term(ring);
        term.set(quot_mono, quot_coeff);

        res.quotients[i] += term;
        p -= term * divisors[i];

        division_occurred = true;
        break;
      }

      if (!division_occurred) {
        Polynomial<T> lt_p(ring);
        lt_p.set(lm_p, lc_p);

        res.remainder += lt_p;
        p -= lt_p;
      }
    }

    return res;
  }

} // namespace order
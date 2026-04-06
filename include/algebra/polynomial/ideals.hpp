#pragma once

#include "io/latex_serializables.hpp"
#include "algebra/polynomial/division.hpp"
#include "algebra/polynomial/monomial.hpp"
#include "algebra/polynomial/monomial_order.hpp"
#include "algebra/polynomial/poly_ring.hpp"
#include "algebra/polynomial/polynomial.hpp"

#include <vector>
#include <utility>
#include <memory>
#include <string>

using Ring = PolyRing<double>;
using Poly = Polynomial<double>;
using RingPtr = std::shared_ptr<const Ring>;

struct Gen {
  int a = 0, b = 0;
};

struct DivisionResult {
  std::vector<std::string> quotient_strs;
  std::string remainder_str;
  bool remainder_zero = false;
};

inline bool in_ideal(const std::vector<Gen>& gens, int m, int n) {
  for (const auto& g : gens)
    if (m >= g.a && n >= g.b) return true;
  return false;
}

inline Monomial make_mono2(int a, int b) {
  return Monomial(Monomial::container{a, b});
}

inline Poly make_gen_poly(RingPtr ring, int a, int b) {
  Poly p(std::move(ring));
  p.set(make_mono2(a, b), 1.0);
  return p;
}

inline std::vector<std::pair<int, int>> remainder_monomials(
  const std::vector<Gen>& gens, int N) {
  std::vector<std::pair<int, int>> res;
  for (int m = 0; m <= N; ++m)
    for (int n = 0; n <= N; ++n)
      if (!in_ideal(gens, m, n))
        res.emplace_back(m, n);
  return res;
}

inline DivisionResult compute_division(
  const std::vector<Gen>& gens,
  const RingPtr& ring,
  int dm, int dn) {
  Poly f(ring);
  f.set(make_mono2(dm, dn), 1.0);

  std::vector<Poly> divisors;
  for (const auto& g : gens)
    divisors.push_back(make_gen_poly(ring, g.a, g.b));

  auto res = order::divide<order::Lex>(f, divisors);

  DivisionResult out;
  out.remainder_zero = res.remainder.is_zero();
  out.remainder_str = out.remainder_zero ? "0" : to_latex(res.remainder);

  for (size_t i = 0; i < gens.size(); ++i) {
    bool z = res.quotients[i].is_zero();
    out.quotient_strs.push_back(z ? "0" : to_latex(res.quotients[i]));
  }
  return out;
}

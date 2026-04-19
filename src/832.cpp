#include <gmpxx.h>

#include "algebra/polynomial/polynomial.hpp"
#include "algebra/polynomial/groebner.hpp"
#include <chalkboard/task.h>
#include <chalkboard/latex_list.h>
#include <algebra/polynomial/poly_ring_pool.hpp>

inline Monomial make_mono(std::initializer_list<int> exps) {
  return Monomial(Monomial::container(exps));
}

void ideal_membership(chalkboard::Reporter& r) {
  r.section("Exercise 1: " +
    r.math("f = xy^3 - z^2 + y^5 - z^3") +
    ",  " +
    r.math(R"(I = \langle -x^3 + y,\, x^2 y - z \rangle)"));

  {
    auto R3 = make_ring<mpq_class>({"x", "y", "z"});

    Polynomial<mpq_class> f(R3);
    f.set(make_mono({1, 3, 0}), mpq_class(1));
    f.set(make_mono({0, 0, 2}), mpq_class(-1));
    f.set(make_mono({0, 5, 0}), mpq_class(1));
    f.set(make_mono({0, 0, 3}), mpq_class(-1));

    Polynomial<mpq_class> g1(R3);
    g1.set(make_mono({3, 0, 0}), mpq_class(-1));
    g1.set(make_mono({0, 1, 0}), mpq_class(1));

    Polynomial<mpq_class> g2(R3);
    g2.set(make_mono({2, 1, 0}), mpq_class(1));
    g2.set(make_mono({0, 0, 1}), mpq_class(-1));

    std::vector<Polynomial<mpq_class>> gens = {g1, g2};

    r.text("f = {}", f);
    r.text("Generators: {} = {},  {} = {}",
           r.math("g_1"), g1, r.math("g_2"), g2);

    auto gb = groebner_basis<order::Lex>(gens);

    r.subsection("Groebner basis G (lex)");
    {
      auto lst = chalkboard::LatexList::enumerate();
      for (std::size_t i = 0; i < gb.size(); i++) {
        lst.item("{} = {}", r.math("g_{" + std::to_string(i + 1) + "}"), gb[i]);
      }
      r.add(lst);
    }

    const auto res = order::divide<order::Lex>(f, gb);

    r.subsection("Reduction of f by G");
    r.text("remainder = {}", res.remainder.is_zero() ? r.math("0") : std::format("{}", res.remainder));

    if (res.remainder.is_zero()) {
      r.text("Conclusion: {} {} {} {} {}",
             r.math("f"), r.math("\\in"), r.math("I"),
             " -- remainder is ", r.math("0"));
    }
    else {
      r.text("Conclusion: {} {} {} {} {}",
             r.math("f"), r.math("\\notin"), r.math("I"),
             " -- remainder is non-zero: ", res.remainder);
    }
  }

  r.section("Exercise 2: " +
    r.math("f = x^3 z - 2y^2") +
    ",  " +
    r.math(R"(I = \langle xz - y,\, xy + 2z^2,\, y - z \rangle)"));

  {
    auto R3 = make_ring<mpq_class>({"x", "y", "z"});

    Polynomial<mpq_class> f(R3);
    f.set(make_mono({3, 0, 1}), mpq_class(1));
    f.set(make_mono({0, 2, 0}), mpq_class(-2));

    Polynomial<mpq_class> g1(R3);
    g1.set(make_mono({1, 0, 1}), mpq_class(1));
    g1.set(make_mono({0, 1, 0}), mpq_class(-1));

    Polynomial<mpq_class> g2(R3);
    g2.set(make_mono({1, 1, 0}), mpq_class(1));
    g2.set(make_mono({0, 0, 2}), mpq_class(2));

    Polynomial<mpq_class> g3(R3);
    g3.set(make_mono({0, 1, 0}), mpq_class(1));
    g3.set(make_mono({0, 0, 1}), mpq_class(-1));

    std::vector<Polynomial<mpq_class>> gens = {g1, g2, g3};

    r.text("f = {}", f);
    r.text("Generators: {} = {},  {} = {},  {} = {}",
           r.math("g_1"), g1, r.math("g_2"), g2, r.math("g_3"), g3);

    auto gb = groebner_basis<order::Lex>(gens);

    r.subsection("Groebner basis G (lex)");
    {
      auto lst = chalkboard::LatexList::enumerate();
      for (std::size_t i = 0; i < gb.size(); i++) {
        lst.item("{} = {}", r.math("g_{" + std::to_string(i + 1) + "}"), gb[i]);
      }
      r.add(lst);
    }

    const auto res = order::divide<order::Lex>(f, gb);

    r.subsection("Reduction of f by G");
    r.text("remainder = {}", res.remainder.is_zero() ? r.math("0") : std::format("{}", res.remainder));

    if (res.remainder.is_zero()) {
      r.text("Conclusion: {} {} {}{}",
             r.math("f"), r.math("\\in"), r.math("I"),
             " -- remainder is 0");
    }
    else {
      r.text("Conclusion: {} {} {}{}",
             r.math("f"), r.math("\\notin"), r.math("I"),
             " -- remainder is non-zero");
    }
  }
}

int main() {
  chalkboard::Task task("Ideal Membership");
  task.build_and_publish(ideal_membership);
  return 0;
}

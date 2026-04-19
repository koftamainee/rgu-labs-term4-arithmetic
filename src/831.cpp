#include <gmpxx.h>

#include "algebra/polynomial/polynomial.hpp"
#include "algebra/polynomial/groebner.hpp"
#include <chalkboard/task.h>
#include <chalkboard/latex_list.h>
#include <algebra/polynomial/poly_ring_pool.hpp>

inline Monomial make_mono(std::initializer_list<int> exps) {
  return Monomial(Monomial::container(exps));
}

template <typename T>
void show_basis(chalkboard::Reporter& r,
                const std::vector<Polynomial<T>>& basis,
                const std::string& label) {
  auto lst = chalkboard::LatexList::enumerate();
  for (std::size_t i = 0; i < basis.size(); i++) {
    lst.item("{} = {}", r.math(label + "_{" + std::to_string(i + 1) + "}"),
             basis[i]);
  }
  r.add(lst);
}

template <typename Order, typename T>
void show_s_polynomials(chalkboard::Reporter& r,
                        const std::vector<Polynomial<T>>& basis) {
  for (std::size_t i = 0; i < basis.size(); i++) {
    for (std::size_t j = i + 1; j < basis.size(); j++) {
      const Polynomial<T> s = s_polynomial<Order>(basis[i], basis[j]);
      const auto res = order::divide<Order>(s, basis);
      r.text("\\(S(f_{{{}}}, f_{{{}}})\\) = {},  remainder = {}",
             i + 1, j + 1,
             s,
             res.remainder.is_zero() ? r.math("0") : std::format("{}", res.remainder));
    }
  }
}

void groebner_exercises(chalkboard::Reporter& r) {
  r.section("Exercise 1: Verify " + r.math("\\overline{S(f_i,f_j)}^{F} = 0") +
    " for all pairs " + r.math("1 \\le i < j \\le s"));

  {
    auto R2 = make_ring<mpq_class>({"x", "y"});

    Polynomial<mpq_class> f1(R2);
    f1.set(make_mono({3, 0}), mpq_class(1));
    f1.set(make_mono({1, 1}), mpq_class(-2));

    Polynomial<mpq_class> f2(R2);
    f2.set(make_mono({2, 1}), mpq_class(1));
    f2.set(make_mono({0, 2}), mpq_class(-2));
    f2.set(make_mono({1, 0}), mpq_class(1));

    std::vector<Polynomial<mpq_class>> F = {f1, f2};
    F = groebner_basis<order::Grlex>(F);

    r.text("F = {{ {}, {} }}", f1, f2);
    r.subsection("S-polynomials and remainders (grlex)");
    show_s_polynomials<order::Grlex>(r, F);
  }

  r.section("Exercise 2(a): " +
    r.math(R"(I = \langle x^2 y - 1,\, xy^2 - x \rangle)"));

  {
    auto R2 = make_ring<mpq_class>({"x", "y"});

    Polynomial<mpq_class> g1(R2);
    g1.set(make_mono({2, 1}), mpq_class(1));
    g1.set(make_mono({0, 0}), mpq_class(-1));

    Polynomial<mpq_class> g2(R2);
    g2.set(make_mono({1, 2}), mpq_class(1));
    g2.set(make_mono({1, 0}), mpq_class(-1));

    std::vector<Polynomial<mpq_class>> gens = {g1, g2};

    r.text("Generators: {} = {},  {} = {}",
           r.math("f_1"), g1, r.math("f_2"), g2);

    r.subsection("Groebner basis (lex)");
    auto gb_lex = groebner_basis<order::Lex>(gens);
    show_basis(r, gb_lex, "g");

    r.subsection("Groebner basis (grlex)");
    auto gb_grlex = groebner_basis<order::Grlex>(gens);
    show_basis(r, gb_grlex, "g");

    r.subsection("Reduced basis (lex)");
    auto rb_lex = reduced_basis<order::Lex>(buchberger<order::Lex>(gens));
    show_basis(r, rb_lex, "h");

    r.subsection("Reduced basis (grlex)");
    auto rb_grlex = reduced_basis<order::Grlex>(buchberger<order::Grlex>(gens));
    show_basis(r, rb_grlex, "h");
  }

  r.section("Exercise 2(b): " +
    r.math(R"(I = \langle x^2 + y,\, x^4 + 2x^2 y + y^2 + 3 \rangle)"));

  {
    auto R2 = make_ring<mpq_class>({"x", "y"});

    Polynomial<mpq_class> g1(R2);
    g1.set(make_mono({2, 0}), mpq_class(1));
    g1.set(make_mono({0, 1}), mpq_class(1));

    Polynomial<mpq_class> g2(R2);
    g2.set(make_mono({4, 0}), mpq_class(1));
    g2.set(make_mono({2, 1}), mpq_class(2));
    g2.set(make_mono({0, 2}), mpq_class(1));
    g2.set(make_mono({0, 0}), mpq_class(3));

    std::vector<Polynomial<mpq_class>> gens = {g1, g2};

    r.text("Generators: {} = {},  {} = {}",
           r.math("f_1"), g1, r.math("f_2"), g2);

    r.subsection("Groebner basis (lex)");
    auto gb_lex = groebner_basis<order::Lex>(gens);
    show_basis(r, gb_lex, "g");

    r.subsection("Groebner basis (grlex)");
    auto gb_grlex = groebner_basis<order::Grlex>(gens);
    show_basis(r, gb_grlex, "g");

    r.subsection("Reduced basis (lex)");
    auto rb_lex = reduced_basis<order::Lex>(buchberger<order::Lex>(gens));
    show_basis(r, rb_lex, "h");

    r.subsection("Reduced basis (grlex)");
    auto rb_grlex = reduced_basis<order::Grlex>(buchberger<order::Grlex>(gens));
    show_basis(r, rb_grlex, "h");
  }

  r.section("Exercise 2(c): " +
    r.math(R"(I = \langle x - z^4,\, y - z^5 \rangle)"));

  {
    auto R3 = make_ring<mpq_class>({"x", "y", "z"});

    Polynomial<mpq_class> g1(R3);
    g1.set(make_mono({1, 0, 0}), mpq_class(1));
    g1.set(make_mono({0, 0, 4}), mpq_class(-1));

    Polynomial<mpq_class> g2(R3);
    g2.set(make_mono({0, 1, 0}), mpq_class(1));
    g2.set(make_mono({0, 0, 5}), mpq_class(-1));

    std::vector<Polynomial<mpq_class>> gens = {g1, g2};

    r.text("Generators: {} = {},  {} = {}",
           r.math("f_1"), g1, r.math("f_2"), g2);

    r.subsection("Groebner basis (lex)");
    auto gb_lex = buchberger<order::Lex>(gens);
    show_basis(r, gb_lex, "g");

    r.subsection("Reduced basis (lex)");
    auto rb_lex = reduced_basis<order::Lex>(buchberger<order::Lex>(gens));
    show_basis(r, rb_lex, "h");

    r.subsection("Groebner basis (grlex)");
    auto gb_grlex = buchberger<order::Grlex>(gens);
    show_basis(r, gb_grlex, "g");

    r.subsection("Reduced basis (grlex)");
    auto rb_grlex = reduced_basis<order::Grlex>(buchberger<order::Grlex>(gens));
    show_basis(r, rb_grlex, "h");
  }
}

int main() {
  chalkboard::Task task("Groebner Bases 1");
  task.build_and_publish(groebner_exercises);
  return 0;
}

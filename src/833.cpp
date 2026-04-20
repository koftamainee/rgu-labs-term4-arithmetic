#include <gmpxx.h>

#include "algebra/polynomial/polynomial.hpp"
#include "algebra/polynomial/groebner.hpp"
#include <chalkboard/task.h>
#include <chalkboard/latex_list.h>
#include <algebra/polynomial/poly_ring_pool.hpp>

inline Monomial make_mono(std::initializer_list<int> exps) {
  return Monomial(Monomial::container(exps));
}

void groebner_tasks(chalkboard::Reporter& r) {
  r.section("Exercise 21.20");

  {
    auto R2 = make_ring<mpq_class>({"x", "y"});

    Polynomial<mpq_class> g1(R2);
    g1.set(make_mono({2, 0}), mpq_class(1));
    g1.set(make_mono({0, 1}), mpq_class(1));
    g1.set(make_mono({0, 0}), mpq_class(-1));

    Polynomial<mpq_class> g2(R2);
    g2.set(make_mono({1, 1}), mpq_class(1));
    g2.set(make_mono({1, 0}), mpq_class(-1));

    std::vector<Polynomial<mpq_class>> gens = {g1, g2};

    r.text("Generators: {} = {},  {} = {}",
           r.math("g_1"), g1,
           r.math("g_2"), g2);

    auto gb = groebner_basis<order::Lex>(gens);

    r.subsection("Groebner basis G (lex, x > y)");
    {
      auto lst = chalkboard::LatexList::enumerate();
      for (std::size_t i = 0; i < gb.size(); i++) {
        lst.item("{} = {}", r.math("g_{" + std::to_string(i + 1) + "}"), gb[i]);
      }
      r.add(lst);
    }

    Polynomial<mpq_class> f1(R2);
    f1.set(make_mono({2, 0}), mpq_class(1));
    f1.set(make_mono({0, 2}), mpq_class(1));
    f1.set(make_mono({0, 1}), mpq_class(-1));

    Polynomial<mpq_class> f2(R2);
    f2.set(make_mono({1, 2}), mpq_class(3));
    f2.set(make_mono({1, 1}), mpq_class(-4));
    f2.set(make_mono({1, 0}), mpq_class(1));
    f2.set(make_mono({0, 0}), mpq_class(1));

    r.subsection("Membership test");

    {
      const auto res = order::divide<order::Lex>(f1, gb);
      r.text("{} remainder = {}", r.math("f_1"),
             res.remainder.is_zero() ? r.math("0") : std::format("{}", res.remainder));
    }

    {
      const auto res = order::divide<order::Lex>(f2, gb);
      r.text("{} remainder = {}", r.math("f_2"),
             res.remainder.is_zero() ? r.math("0") : std::format("{}", res.remainder));
    }
  }

  r.section("Exercise 21.23");

  {
    auto R3 = make_ring<mpq_class>({"x", "y", "z"});

    Polynomial<mpq_class> f1(R3);
    f1.set(make_mono({2, 1, 0}), mpq_class(1));
    f1.set(make_mono({0, 1, 1}), mpq_class(-2));
    f1.set(make_mono({0, 0, 0}), mpq_class(1));

    Polynomial<mpq_class> f2(R3);
    f2.set(make_mono({1, 2, 0}), mpq_class(1));
    f2.set(make_mono({0, 0, 2}), mpq_class(-1));
    f2.set(make_mono({1, 0, 0}), mpq_class(2));

    Polynomial<mpq_class> f3(R3);
    f3.set(make_mono({0, 2, 1}), mpq_class(1));
    f3.set(make_mono({2, 0, 0}), mpq_class(-1));
    f3.set(make_mono({0, 0, 0}), mpq_class(5));

    std::vector<Polynomial<mpq_class>> gens = {f1, f2, f3};

    r.text("Generators: {} = {},  {} = {},  {} = {}",
           r.math("f_1"), f1,
           r.math("f_2"), f2,
           r.math("f_3"), f3);


    auto gb = buchberger<order::Grlex>(gens);

    auto print_basis = [&]() {
      auto lst = chalkboard::LatexList::enumerate();
      for (std::size_t i = 0; i < gb.size(); i++) {
        lst.item("{} = {}", r.math("g_{" + std::to_string(i + 1) + "}"), gb[i]);
      }
      r.add(lst);
    };

    r.subsection("Groebner basis G (grlex, x < y < z)");
    print_basis();

    gb = minimal_basis<order::Grlex>(gb);

    r.subsection("Minimal Groebner basis G (grlex, x < y < z)");
    print_basis();

    gb = reduced_basis<order::Grlex>(gb);
    r.subsection("Reduced Groebner basis G (grlex, x < y < z)");
    print_basis();
  }
}

int main() {
  chalkboard::Task task("Computer algebra");
  task.build_and_publish(groebner_tasks);
  return 0;
}

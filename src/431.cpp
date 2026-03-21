#include <complex>
#include <iostream>
#include <gmpxx.h>

#include "polynomial.hpp"
#include <chalkboard/task.h>
#include <chalkboard/latex_list.h>

inline Monomial make_mono(std::initializer_list<int> exps) {
  return Monomial(Monomial::container(exps));
}

void poly_demo(chalkboard::Reporter& r) {
  auto R3 = make_ring<mpq_class>({"x", "y", "z"});

  Polynomial<mpq_class> f(R3);
  f.set(make_mono({2, 1, 0}), mpq_class(3));
  f.set(make_mono({1, 0, 1}), mpq_class(-2));
  f.set(make_mono({0, 3, 0}), mpq_class(1, 2));
  f.set(make_mono({0, 0, 0}), mpq_class(5, 3));

  Polynomial<mpq_class> g(R3);
  g.set(make_mono({2, 1, 0}), mpq_class(1, 4));
  g.set(make_mono({1, 0, 1}), mpq_class(7, 3));
  g.set(make_mono({0, 0, 0}), mpq_class(-4));

  Polynomial<mpq_class> h_hom(R3);
  h_hom.set(make_mono({2, 1, 0}), mpq_class(1));
  h_hom.set(make_mono({0, 0, 3}), mpq_class(7, 2));

  std::vector<mpq_class> A = {1.0, 2.0, -1.0};


  using CQ = std::complex<mpq_class>;
  auto CQ2 = make_ring<CQ>({"u", "v"});

  Polynomial<CQ> p(CQ2);
  p.set(make_mono({2, 0}), CQ{mpq_class(1), mpq_class(1)});
  p.set(make_mono({0, 1}), CQ{mpq_class(0), mpq_class(-2)});
  p.set(make_mono({0, 0}), CQ{mpq_class(3, 2), mpq_class(0)});

  Polynomial<CQ> q(CQ2);
  q.set(make_mono({1, 1}), CQ{mpq_class(1), mpq_class(0)});
  q.set(make_mono({0, 0}), CQ{mpq_class(-1, 1), mpq_class(1, 3)});



  std::vector<CQ> B = {
    CQ{mpq_class(1), mpq_class(1)},
    CQ{mpq_class(0), mpq_class(-1)}
  };


  auto Z2 = make_ring<mpz_class>({"\\alpha", "\\beta"});

  Polynomial<mpz_class> poly_r(Z2);
  poly_r.set(make_mono({3, 0}), mpz_class("7"));
  poly_r.set(make_mono({1, 2}), mpz_class("-5"));
  poly_r.set(make_mono({0, 0}), mpz_class("42"));

  Polynomial<mpz_class> poly_s(Z2);
  poly_s.set(make_mono({2, 1}), mpz_class("3"));
  poly_s.set(make_mono({0, 1}), mpz_class("-2"));
  poly_s.set(make_mono({0, 0}), mpz_class("10"));

  std::vector<mpz_class> C_pt = {mpz_class(2), mpz_class(3)};

  r.section(r.math(r.ring(r.R(), {"x", "y", "z"})));
  r.text("{} = {}", r.math("f"), f);
  r.text("{} = {}", r.math("g"), g);

  r.subsection("(a) Arithmetic");
  r.text("{} = {}", r.math("f + g"), f + g);
  r.text("{} = {}", r.math("f - g"), f - g);
  r.text("{} = {}", r.math(r.cdot("f", "g")), f * g);

  r.subsection("(b) " + r.math(r.mrm("supp") + "(f)"));
  {
    auto lst = chalkboard::LatexList::itemize();
    for (const auto& m : f.supp())
      lst.item("{}", m);
    r.add(lst);
  }

  r.subsection("(c) Equality");
  r.text("{} {}", r.math(r.eq("f", "g")), f == g);
  r.text("{} {}", r.math(r.eq("f", "f")), f == f);
  r.text("{} {}", r.math(r.neq("f", "g")), f != g);

  r.subsection("(d) " + r.math("f(1, 2, -1)") + ",  " + r.math("g(1, 2, -1)"));
  r.text("{} = {}", r.math("f(A)"), f.evaluate(A));
  r.text("{} = {}", r.math("g(A)"), g.evaluate(A));

  r.subsection("(e) Homogeneous degree");
  r.text("{} is not homogeneous", r.math("f"));
  r.text("{} = {},  {} = {}", r.math("h"), h_hom,
         r.math(r.sub_script(r.mrm("deg"), r.mrm("hom")) + "(h)"),
         *h_hom.homogeneous_degree());

  r.subsection("(f) Degree-3 component of " + r.math("f"));
  r.math_block(f.homogeneous_component(3));


  r.section(r.math(r.ring(r.C() + "(" + r.Q() + ")", {"u", "v"})));
  r.text("{} = {}", r.math("p"), p);
  r.text("{} = {}", r.math("q"), q);

  r.subsection("(a) Arithmetic");
  r.text("{} = {}", r.math("p + q"), p + q);
  r.text("{} = {}", r.math(r.cdot("p", "q")), p * q);

  r.subsection("(b) " + r.math(r.mrm("supp") + "(p)"));
  {
    auto lst = chalkboard::LatexList::itemize();
    for (const auto& m : p.supp())
      lst.item("{}", m);
    r.add(lst);
  }

  r.subsection("(c) Equality");
  r.text("{} {}", r.math(r.eq("p", "q")), p == q);
  r.text("{} {}", r.math(r.eq("p", "p")), p == p);

  r.subsection("(d) " + r.math("p(1+i,\\, -i)"));
  r.text("{} = {}", r.math("p(B)"), p.evaluate(B));

  r.subsection("(e, f) Degree-2 component of " + r.math("p"));
  r.math_block(p.homogeneous_component(2));


  r.section(r.math(r.ring(r.Z(), {"\\alpha", "\\beta"})));
  r.text("{} = {}", r.math("r"), poly_r);
  r.text("{} = {}", r.math("s"), poly_s);

  r.subsection("(a) Arithmetic");
  r.text("{} = {}", r.math("r + s"), poly_r + poly_s);
  r.text("{} = {}", r.math(r.cdot("r", "s")), poly_r * poly_s);

  r.subsection("(b) " + r.math(r.mrm("supp") + "(r)"));
  {
    auto lst = chalkboard::LatexList::itemize();
    for (const auto& m : poly_r.supp())
      lst.item("{}", m);
    r.add(lst);
  }

  r.subsection("(c) Equality");
  r.text("{} {}", r.math(r.eq("r", "s")), poly_r == poly_s);
  r.text("{} {}", r.math(r.eq("r", "r")), poly_r == poly_r);

  r.subsection("(d) " + r.math("r(2, 3)"));
  r.text("{} = {}", r.math("r(C)"), poly_r.evaluate(C_pt));

  r.subsection("(e, f) Degree-3 component of " + r.math("r"));
  r.math_block(poly_r.homogeneous_component(3));
}

int main() {
  chalkboard::Task task("Polynomial operations");
  task.build_and_publish(poly_demo);

  return 0;
}

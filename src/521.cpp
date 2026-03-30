#include <gmpxx.h>

#include "algebra/polynomial/monomial_order.hpp"
#include "algebra/polynomial/polynomial.hpp"
#include "algebra/polynomial/ordering.hpp"
#include <chalkboard/task.h>

inline Monomial make_mono(const std::initializer_list<int> exps) {
  return Monomial(Monomial::container(exps));
}

using Q = mpq_class;

template <typename Order, typename T>
void show_order_full(chalkboard::Reporter& r,
                     const std::string& order_name,
                     const Polynomial<T>& f) {
  r.subsubsection(order_name);

  std::vector<Monomial> monomials = f.supp();
  std::sort(monomials.begin(), monomials.end(), order::Comparator<Order>{});

  std::string sorted_str;
  bool first = true;
  for (const auto& m : monomials) {
    if (!first) sorted_str += " \\prec_{" + order_name + "} ";
    sorted_str += to_latex(m);
    first = false;
  }
  r.text("Sorted monomials: {}", r.math(sorted_str));

  r.text("{} = {}",
         r.math(r.sub_script(r.mrm("mdeg"), order_name) + "(f)"),
         order::multideg<Order>(f));
  r.text("{} = {}",
         r.math(r.sub_script(r.mrm("lc"), order_name) + "(f)"),
         order::lc<Order>(f));
  r.text("{} = {}",
         r.math(r.sub_script(r.mrm("lm"), order_name) + "(f)"),
         order::lm<Order>(f));
  r.text("{} = {}",
         r.math(r.sub_script(r.mrm("lt"), order_name) + "(f)"),
         order::lt<Order>(f));
}

void ex21_6(chalkboard::Reporter& r) {
  r.section("Exercise 21.6");

  // f = 2x^4y^2z - 6x^4yz^2 + 4xy^4z^2 - 3xy^2z^4 + x^2y^4z - 5x^2yz^4
  const auto R3 = make_ring<Q>({"x", "y", "z"});

  Polynomial<Q> f(R3);
  f.set(make_mono({4, 2, 1}), Q(2));
  f.set(make_mono({4, 1, 2}), Q(-6));
  f.set(make_mono({1, 4, 2}), Q(4));
  f.set(make_mono({1, 2, 4}), Q(-3));
  f.set(make_mono({2, 4, 1}), Q(1));
  f.set(make_mono({2, 1, 4}), Q(-5));

  r.text("Let {}", r.math("f = " + to_latex(f) + " \\in " + r.Q() + "[x,y,z]"));
  r.text("with {} in all cases.", r.math("x \\succ y \\succ z"));

  r.subsection("(i) & (ii) — lex, grlex, grevlex");

  show_order_full<order::Lex>(r, "lex", f);
  show_order_full<order::Grlex>(r, "grlex", f);
  show_order_full<order::Grevlex>(r, "grevlex", f);
}

int main() {
  const chalkboard::Task task("Exercise 21.6");
  task.build_and_publish(ex21_6);
  return 0;
}

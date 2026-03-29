#include <gmpxx.h>

#include "algebra/polynomial/monomial_order.hpp"
#include "algebra/polynomial/polynomial.hpp"
#include <chalkboard/task.h>

#include "algebra/polynomial/ordering.hpp"

inline Monomial make_mono(const std::initializer_list<int> exps) {
  return Monomial(Monomial::container(exps));
}

using Q = mpq_class;

template <typename T>
void show_lm_lt_multideg(chalkboard::Reporter& r,
                         const std::string& order_name,
                         const Polynomial<T>& f,
                         [[maybe_unused]] auto tag) {
  using Order = decltype(tag);
  r.subsubsection(order_name);
  r.text("{} = {}",
         r.math(r.sub_script(r.mrm("multideg"), order_name) + "(f)"),
         order::multideg<Order>(f));
  r.text("{} = {}",
         r.math(r.sub_script(r.mrm("lm"), order_name) + "(f)"),
         order::lm<Order>(f));
  r.text("{} = {}",
         r.math(r.sub_script(r.mrm("lt"), order_name) + "(f)"),
         order::lt<Order>(f));
}

template <typename T>
void show_all_orders(chalkboard::Reporter& r, const Polynomial<T>& f) {
  r.math_block(f);
  show_lm_lt_multideg(r, "lex", f, order::Lex{});
  show_lm_lt_multideg(r, "grlex", f, order::Grlex{});
  show_lm_lt_multideg(r, "grevlex", f, order::Grevlex{});
}

void identify_order(chalkboard::Reporter& r,
                    const std::string& label,
                    const Polynomial<Q>& f,
                    const Monomial& first_term_in_book) {
  r.subsection(label);
  r.math_block(f);

  auto lm_lex = order::lm<order::Lex>(f);
  auto lm_grlex = order::lm<order::Grlex>(f);
  auto lm_grevlex = order::lm<order::Grevlex>(f);

  r.text("lex     {} = {}",
         r.math(r.sub_script(r.mrm("lm"), "lex") + "(f)"), lm_lex);
  r.text("grlex   {} = {}",
         r.math(r.sub_script(r.mrm("lm"), "grlex") + "(f)"), lm_grlex);
  r.text("grevlex {} = {}",
         r.math(r.sub_script(r.mrm("lm"), "grevlex") + "(f)"), lm_grevlex);

  r.text("first term in book: {}", first_term_in_book);

  if (lm_lex == first_term_in_book)
    r.text("ordering: lex");
  else if (lm_grlex == first_term_in_book)
    r.text("ordering: grlex");
  else if (lm_grevlex == first_term_in_book)
    r.text("ordering: grevlex");
  else
    r.text("ordering: none of lex/grlex/grevlex matches");
}

void ordering_demo(chalkboard::Reporter& r) {
  r.section("Ch.2 " + r.math("\\S 2") + " Exercise 1");
  {
    const auto R3 = make_ring<Q>({"x", "y", "z"});

    r.subsection("(a) " + r.math("f = 2x + 3y + z + x^2 - z^2 + x^3"));
    Polynomial<Q> fa(R3);
    fa.set(make_mono({1, 0, 0}), Q(2));
    fa.set(make_mono({0, 1, 0}), Q(3));
    fa.set(make_mono({0, 0, 1}), Q(1));
    fa.set(make_mono({2, 0, 0}), Q(1));
    fa.set(make_mono({0, 0, 2}), Q(-1));
    fa.set(make_mono({3, 0, 0}), Q(1));
    show_all_orders(r, fa);

    r.subsection("(b) " + r.math("f = 2x^2y^8 - 3x^5yz^4 + xyz^3 - xy^4"));
    Polynomial<Q> fb(R3);
    fb.set(make_mono({2, 8, 0}), Q(2));
    fb.set(make_mono({5, 1, 4}), Q(-3));
    fb.set(make_mono({1, 1, 3}), Q(1));
    fb.set(make_mono({1, 4, 0}), Q(-1));
    show_all_orders(r, fb);
  }

  r.section("Ch.2 " + r.math("\\S 2") + " Exercise 2 — identify the ordering");
  {
    const auto R3 = make_ring<Q>({"x", "y", "z"});

    // (a) 7x^2y^4z - 2xy^6 + x^2y^2  — first term x^2y^4z
    Polynomial<Q> fa(R3);
    fa.set(make_mono({2, 4, 1}), Q(7));
    fa.set(make_mono({1, 6, 0}), Q(-2));
    fa.set(make_mono({2, 2, 0}), Q(1));
    identify_order(r,
                   "(a) " + r.math("7x^2y^4z - 2xy^6 + x^2y^2"),
                   fa,
                   make_mono({2, 4, 1}));

    // (b) xy^3z + xy^2z^2 + x^2z^3  — first term xy^3z
    Polynomial<Q> fb(R3);
    fb.set(make_mono({1, 3, 1}), Q(1));
    fb.set(make_mono({1, 2, 2}), Q(1));
    fb.set(make_mono({2, 0, 3}), Q(1));
    identify_order(r,
                   "(b) " + r.math("xy^3z + xy^2z^2 + x^2z^3"),
                   fb,
                   make_mono({1, 3, 1}));

    // (c) x^4y^5z + 2x^3y^2z - 4xy^2z^4  — first term x^4y^5z
    Polynomial<Q> fc(R3);
    fc.set(make_mono({4, 5, 1}), Q(1));
    fc.set(make_mono({3, 2, 1}), Q(2));
    fc.set(make_mono({1, 2, 4}), Q(-4));
    identify_order(r,
                   "(c) " + r.math("x^4y^5z + 2x^3y^2z - 4xy^2z^4"),
                   fc,
                   make_mono({4, 5, 1}));
  }

  r.section("Ch.2 " + r.math("\\S 2") + " Exercise 3 — ordering "
    + r.math("z > y > x"));
  {
    const auto R3 = make_ring<Q>({"z", "y", "x"});

    r.subsection("(a) " + r.math("f = 2x + 3y + z + x^2 - z^2 + x^3"));
    Polynomial<Q> fa(R3);
    fa.set(make_mono({0, 0, 1}), Q(2));
    fa.set(make_mono({0, 1, 0}), Q(3));
    fa.set(make_mono({1, 0, 0}), Q(1));
    fa.set(make_mono({0, 0, 2}), Q(1));
    fa.set(make_mono({2, 0, 0}), Q(-1));
    fa.set(make_mono({0, 0, 3}), Q(1));
    show_all_orders(r, fa);

    r.subsection("(b) " + r.math("f = 2x^2y^8 - 3x^5yz^4 + xyz^3 - xy^4"));
    Polynomial<Q> fb(R3);
    fb.set(make_mono({0, 8, 2}), Q(2));
    fb.set(make_mono({4, 1, 5}), Q(-3));
    fb.set(make_mono({3, 1, 1}), Q(1));
    fb.set(make_mono({0, 4, 1}), Q(-1));
    show_all_orders(r, fb);
  }
}

int main() {
  const chalkboard::Task task("Monomial orderings");
  task.build_and_publish(ordering_demo);
  return 0;
}

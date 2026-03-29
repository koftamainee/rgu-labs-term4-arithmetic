#include <gmpxx.h>

#include "algebra/polynomial/division.hpp"
#include "algebra/polynomial/ordering.hpp"
#include "algebra/polynomial/monomial_order.hpp"
#include "algebra/polynomial/polynomial.hpp"
#include <chalkboard/latex_list.h>
#include <chalkboard/task.h>

inline Monomial make_mono(std::initializer_list<int> exps) {
    return Monomial(Monomial::container(exps));
}

using Q = mpq_class;

template <typename Order, typename T>
void show_division(chalkboard::Reporter& r,
                   const std::string& order_name,
                   const Polynomial<T>& f,
                   const std::vector<Polynomial<T>>& divisors) {
    auto res = order::divide<Order>(f, divisors);
    auto ring = f.ring();

    r.subsubsection(order_name);

    r.text("{} = {}", r.math("f"), f);

    {
        auto lst = chalkboard::LatexList::itemize();
        for (std::size_t i = 0; i < divisors.size(); ++i) {
            lst.item("{} = {}", r.math("f_" + std::to_string(i + 1)), divisors[i]);
        }
        r.add(lst);
    }

    {
        auto lst = chalkboard::LatexList::itemize();
        for (std::size_t i = 0; i < res.quotients.size(); ++i) {
            lst.item("{} = {}", r.math("a_" + std::to_string(i + 1)), res.quotients[i]);
        }
        r.add(lst);
    }

    r.text("{} = {}", r.math("r"), res.remainder);

    Polynomial<T> check(ring);
    for (std::size_t i = 0; i < divisors.size(); ++i) {
        check += res.quotients[i] * divisors[i];
    }
    check += res.remainder;
    r.text("check {} = {}", r.math("f - (a_1 f_1 + \\ldots + a_s f_s + r)"),
           f - check);
}

void division_demo(chalkboard::Reporter& r) {
    r.section("Ch.2 " + r.math("\\S 3") + " Exercise 1");

    r.subsection("(a) " + r.math("F = (xy^2 - x,\\; x - y^3)"));
    {
        auto R2 = make_ring<Q>({"x", "y"});

        Polynomial<Q> f(R2);
        f.set(make_mono({7, 2}), Q(1));
        f.set(make_mono({3, 2}), Q(1));
        f.set(make_mono({0, 1}), Q(-1));
        f.set(make_mono({0, 0}), Q(1));

        Polynomial<Q> f1(R2);
        f1.set(make_mono({1, 2}), Q(1));
        f1.set(make_mono({1, 0}), Q(-1));

        Polynomial<Q> f2(R2);
        f2.set(make_mono({1, 0}), Q(1));
        f2.set(make_mono({0, 3}), Q(-1));

        std::vector<Polynomial<Q>> F = {f1, f2};

        show_division<order::Grlex>(r, "grlex", f, F);
        show_division<order::Lex>  (r, "lex",   f, F);
    }

    r.subsection("(b) " + r.math("F = (x - y^3,\\; xy^2 - x)") + " (reversed)");
    {
        auto R2 = make_ring<Q>({"x", "y"});

        Polynomial<Q> f(R2);
        f.set(make_mono({7, 2}), Q(1));
        f.set(make_mono({3, 2}), Q(1));
        f.set(make_mono({0, 1}), Q(-1));
        f.set(make_mono({0, 0}), Q(1));

        Polynomial<Q> f1(R2);
        f1.set(make_mono({1, 0}), Q(1));
        f1.set(make_mono({0, 3}), Q(-1));

        Polynomial<Q> f2(R2);
        f2.set(make_mono({1, 2}), Q(1));
        f2.set(make_mono({1, 0}), Q(-1));

        std::vector<Polynomial<Q>> F = {f1, f2};

        show_division<order::Grlex>(r, "grlex", f, F);
        show_division<order::Lex>  (r, "lex",   f, F);
    }

    r.section("Ch.2 " + r.math("\\S 3") + " Exercise 2");

    r.subsection("(a) " + r.math("F = (x - y^2,\\; y - z^3,\\; z^2 - 1)"));
    {
        auto R3 = make_ring<Q>({"x", "y", "z"});

        Polynomial<Q> f(R3);
        f.set(make_mono({1, 2, 2}), Q(1));
        f.set(make_mono({1, 1, 0}), Q(1));
        f.set(make_mono({0, 1, 1}), Q(-1));

        Polynomial<Q> f1(R3);
        f1.set(make_mono({1, 0, 0}), Q(1));
        f1.set(make_mono({0, 2, 0}), Q(-1));

        Polynomial<Q> f2(R3);
        f2.set(make_mono({0, 1, 0}), Q(1));
        f2.set(make_mono({0, 0, 3}), Q(-1));

        Polynomial<Q> f3(R3);
        f3.set(make_mono({0, 0, 2}), Q(1));
        f3.set(make_mono({0, 0, 0}), Q(-1));

        std::vector<Polynomial<Q>> F = {f1, f2, f3};

        show_division<order::Grlex>(r, "grlex", f, F);
        show_division<order::Lex>  (r, "lex",   f, F);
    }

    r.subsection("(b) " + r.math("F = (y - z^3,\\; z^2 - 1,\\; x - y^2)") + " (cyclic)");
    {
        auto R3 = make_ring<Q>({"x", "y", "z"});

        Polynomial<Q> f(R3);
        f.set(make_mono({1, 2, 2}), Q(1));
        f.set(make_mono({1, 1, 0}), Q(1));
        f.set(make_mono({0, 1, 1}), Q(-1));

        Polynomial<Q> f1(R3);
        f1.set(make_mono({0, 1, 0}), Q(1));
        f1.set(make_mono({0, 0, 3}), Q(-1));

        Polynomial<Q> f2(R3);
        f2.set(make_mono({0, 0, 2}), Q(1));
        f2.set(make_mono({0, 0, 0}), Q(-1));

        Polynomial<Q> f3(R3);
        f3.set(make_mono({1, 0, 0}), Q(1));
        f3.set(make_mono({0, 2, 0}), Q(-1));

        std::vector<Polynomial<Q>> F = {f1, f2, f3};

        show_division<order::Grlex>(r, "grlex", f, F);
        show_division<order::Lex>  (r, "lex",   f, F);
    }
}

int main() {
    const chalkboard::Task task("Polynomial division");
    task.build_and_publish(division_demo);
    return 0;
}
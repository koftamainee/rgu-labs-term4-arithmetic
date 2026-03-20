#include <complex>
#include <iostream>
#include <gmpxx.h>

#include "polynomial.hpp"
#include "reporter.hpp"
#include "publisher.hpp"
#include "latex_list.hpp"

static Monomial make_mono(std::initializer_list<int> exps) {
    return Monomial(Monomial::container(exps));
}

int main() {

    // Q[x, y, z]

    auto Q3 = make_ring<mpq_class>({"x", "y", "z"});

    Polynomial<mpq_class> f(Q3);
    f.set(make_mono({2, 1, 0}), mpq_class(3, 1));
    f.set(make_mono({1, 0, 1}), mpq_class(-2, 1));
    f.set(make_mono({0, 3, 0}), mpq_class(1, 2));
    f.set(make_mono({0, 0, 0}), mpq_class(5, 3));

    Polynomial<mpq_class> g(Q3);
    g.set(make_mono({2, 1, 0}), mpq_class(1, 4));
    g.set(make_mono({1, 0, 1}), mpq_class(7, 3));
    g.set(make_mono({0, 0, 0}), mpq_class(-4, 1));

    Polynomial<mpq_class> h_hom(Q3);
    h_hom.set(make_mono({2, 1, 0}), mpq_class(1, 1));
    h_hom.set(make_mono({0, 0, 3}), mpq_class(7, 2));

    std::vector<mpq_class> A = {mpq_class(1), mpq_class(2), mpq_class(-1)};


    // C(Q)[u, v]

    using CQ = std::complex<mpq_class>;
    auto CQ2 = make_ring<CQ>({"u", "v"});

    Polynomial<CQ> p(CQ2);
    p.set(make_mono({2, 0}), CQ{mpq_class(1),  mpq_class(1)});
    p.set(make_mono({0, 1}), CQ{mpq_class(0),  mpq_class(-2)});
    p.set(make_mono({0, 0}), CQ{mpq_class(3, 2), mpq_class(0)});

    Polynomial<CQ> q(CQ2);
    q.set(make_mono({1, 1}), CQ{mpq_class(1),    mpq_class(0)});
    q.set(make_mono({0, 0}), CQ{mpq_class(-1, 1), mpq_class(1, 3)});

    std::vector<CQ> B = {
        CQ{mpq_class(1), mpq_class(1)},
        CQ{mpq_class(0), mpq_class(-1)}
    };


    // Z[alpha, beta]

    auto Z2 = make_ring<mpz_class>({"\\alpha", "\\beta"});

    Polynomial<mpz_class> r(Z2);
    r.set(make_mono({3, 0}), mpz_class("12345678901234567890"));
    r.set(make_mono({1, 2}), mpz_class("-9876543210987654321"));
    r.set(make_mono({0, 0}), mpz_class("42"));

    Polynomial<mpz_class> s(Z2);
    s.set(make_mono({2, 1}), mpz_class("7"));
    s.set(make_mono({0, 1}), mpz_class("-3"));
    s.set(make_mono({0, 0}), mpz_class("100"));

    std::vector<mpz_class> C_pt = {mpz_class(2), mpz_class(3)};


    // Z[i][x]

    using CZ = std::complex<mpz_class>;
    auto CZ1 = make_ring<CZ>({"x"});

    Polynomial<CZ> cx(CZ1);
    cx.set(make_mono({4}), CZ{mpz_class(1),  mpz_class(0)});
    cx.set(make_mono({2}), CZ{mpz_class(0),  mpz_class(3)});
    cx.set(make_mono({1}), CZ{mpz_class(-1), mpz_class(-1)});
    cx.set(make_mono({0}), CZ{mpz_class(2),  mpz_class(0)});

    Polynomial<CZ> cy(CZ1);
    cy.set(make_mono({2}), CZ{mpz_class(1),  mpz_class(1)});
    cy.set(make_mono({0}), CZ{mpz_class(-1), mpz_class(0)});


    // Report

    Reporter rep("Polynomial Ring Operations");

    rep.section("Q[x, y, z]");
    rep.text("f = {}", f);
    rep.text("g = {}", g);

    rep.subsection("(a) Arithmetic");
    rep.text("f + g = {}", f + g);
    rep.text("f - g = {}", f - g);
    rep.text("f * g = {}", f * g);

    rep.subsection("(b) supp(f)");
    {
        auto lst = LatexList::itemize();
        for (const auto& m : f.supp())
            lst.item("{}", m);
        rep.add(lst);
    }

    rep.subsection("(c) Equality");
    rep.text("f == g  ->  {}", f == g);
    rep.text("f == f  ->  {}", f == f);
    rep.text("f != g  ->  {}", f != g);

    rep.subsection("(d) f(1, 2, -1),  g(1, 2, -1)");
    rep.text("f(A) = {}", f.evaluate(A));
    rep.text("g(A) = {}", g.evaluate(A));

    rep.subsection("(e) Homogeneous degree");
    rep.text("f: not homogeneous");
    rep.text("h = {},  deg_hom(h) = {}", h_hom, *h_hom.homogeneous_degree());

    rep.subsection("(f) Degree-3 component of f");
    rep.math(f.homogeneous_component(3));


    rep.section("C(Q)[u, v]");
    rep.text("p = {}", p);
    rep.text("q = {}", q);

    rep.subsection("(a) Arithmetic");
    rep.text("p + q = {}", p + q);
    rep.text("p * q = {}", p * q);

    rep.subsection("(b) supp(p)");
    {
        auto lst = LatexList::itemize();
        for (const auto& m : p.supp())
            lst.item("{}", m);
        rep.add(lst);
    }

    rep.subsection("(c) Equality");
    rep.text("p == q  ->  {}", p == q);
    rep.text("p == p  ->  {}", p == p);

    rep.subsection("(d) p(1+i, -i)");
    rep.text("p(B) = {}", p.evaluate(B));

    rep.subsection("(e,f) Degree-2 component of p");
    rep.math(p.homogeneous_component(2));


    rep.section("Z[alpha, beta]  (arbitrary precision)");
    rep.text("r = {}", r);
    rep.text("s = {}", s);

    rep.subsection("(a) Arithmetic");
    rep.text("r + s = {}", r + s);
    rep.text("r * s = {}", r * s);

    rep.subsection("(b) supp(r)");
    {
        auto lst = LatexList::itemize();
        for (const auto& m : r.supp())
            lst.item("{}", m);
        rep.add(lst);
    }

    rep.subsection("(c) Equality");
    rep.text("r == s  ->  {}", r == s);
    rep.text("r == r  ->  {}", r == r);

    rep.subsection("(d) r(2, 3)");
    rep.text("r(C) = {}", r.evaluate(C_pt));

    rep.subsection("(e,f) Degree-3 component of r");
    rep.math(r.homogeneous_component(3));


    rep.section("Z[i][x]  (Gaussian integers)");
    rep.text("cx = {}", cx);
    rep.text("cy = {}", cy);

    rep.subsection("(a) Arithmetic");
    rep.text("cx + cy = {}", cx + cy);
    rep.text("cx * cy = {}", cx * cy);

    rep.subsection("(b) supp(cx)");
    {
        auto lst = LatexList::itemize();
        for (const auto& m : cx.supp())
            lst.item("{}", m);
        rep.add(lst);
    }

    rep.subsection("(c) Equality");
    rep.text("cx == cy  ->  {}", cx == cy);
    rep.text("cx == cx  ->  {}", cx == cx);

    rep.subsection("(e) Homogeneous degree");
    rep.text("cx: not homogeneous");

    rep.subsection("(f) Degree-4 component of cx");
    rep.math(cx.homogeneous_component(4));


    std::string artifact_dir = rep.build();
    std::cout << "Report built at: " << artifact_dir << "\n";

    try {
        Publisher pub("localhost", 8080);
        pub.publish(artifact_dir);
        std::cout << "Published to http://localhost:8080\n";
    } catch (const std::exception& e) {
        std::cerr << "Publisher: " << e.what() << "\n";
    }

    return 0;
}
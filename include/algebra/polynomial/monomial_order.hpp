#pragma once

#include <compare>
#include <stdexcept>

#include "monomial.hpp"

namespace order {

template <typename Order>
struct OrderTraits;

template <typename Order>
std::strong_ordering compare(const Monomial& a, const Monomial& b) {
    return OrderTraits<Order>::compare(a, b);
}

namespace detail {

inline void check_compat(const Monomial& a, const Monomial& b) {
    if (a.n_vars() != b.n_vars()) {
        throw std::invalid_argument("monomial order compare: n_vars mismatch");
    }
}

} // namespace detail

struct Lex {};
struct Grlex {};
struct Grevlex {};
struct Invlex {};
struct Rinvlex {};

template <>
struct OrderTraits<Lex> {
    static std::strong_ordering compare(const Monomial& a, const Monomial& b) {
        detail::check_compat(a, b);
        for (Monomial::size_type i = 0; i < a.n_vars(); ++i) {
            if (const auto cmp = a[i] <=> b[i]; cmp != 0) {
                return cmp;
            }
        }
        return std::strong_ordering::equal;
    }
};

template <>
struct OrderTraits<Grlex> {
    static std::strong_ordering compare(const Monomial& a, const Monomial& b) {
        detail::check_compat(a, b);
        if (const auto cmp = a.total_degree() <=> b.total_degree(); cmp != 0) {
            return cmp;
        }
        return OrderTraits<Lex>::compare(a, b);
    }
};

template <>
struct OrderTraits<Grevlex> {
    static std::strong_ordering compare(const Monomial& a, const Monomial& b) {
        detail::check_compat(a, b);
        if (const auto cmp = a.total_degree() <=> b.total_degree(); cmp != 0) {
            return cmp;
        }
        for (Monomial::size_type i = a.n_vars(); i-- > 0;) {
            if (const auto cmp = b[i] <=> a[i]; cmp != 0) {
                return cmp;
            }
        }
        return std::strong_ordering::equal;
    }
};

template <>
struct OrderTraits<Invlex> {
    static std::strong_ordering compare(const Monomial& a, const Monomial& b) {
        detail::check_compat(a, b);
        for (Monomial::size_type i = a.n_vars(); i-- > 0;) {
            if (auto cmp = a[i] <=> b[i]; cmp != 0) {
                return cmp;
            }
        }
        return std::strong_ordering::equal;
    }
};

template <>
struct OrderTraits<Rinvlex> {
    static std::strong_ordering compare(const Monomial& a, const Monomial& b) {
        detail::check_compat(a, b);
        for (Monomial::size_type i = 0; i < a.n_vars(); ++i) {
            if (const auto cmp = b[i] <=> a[i]; cmp != 0) {
                return cmp;
            }
        }
        return std::strong_ordering::equal;
    }
};

template <typename Order>
struct Comparator {
    bool operator()(const Monomial& a, const Monomial& b) const {
        return OrderTraits<Order>::compare(a, b) < 0;
    }
};

} // namespace order
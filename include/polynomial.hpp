#ifndef RGU_LABS_TERM4_ARITHMETIC_POLYNOMIAL_HPP
#define RGU_LABS_TERM4_ARITHMETIC_POLYNOMIAL_HPP

#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

#include "TTrie.hpp"
#include "monomial.hpp"
#include "poly_ring.hpp"

template <typename T>
class TPolynomial {
public:
    using coeff_type = T;
    using ring_ptr   = std::shared_ptr<const TPolyRing<T>>;
    using monomial   = Monomial;
    using key_type   = Monomial::container;
    using size_type  = std::size_t;
    using point_type = std::vector<T>;

    explicit TPolynomial(ring_ptr ring);

    ring_ptr ring() const noexcept;
    size_type n_vars() const noexcept;

    void set(const monomial& m, const T& coeff);
    const T* get(const monomial& m) const;
    bool is_zero() const noexcept;

    TPolynomial& operator+=(const TPolynomial& other);
    TPolynomial& operator-=(const TPolynomial& other);
    TPolynomial& operator*=(const TPolynomial& other);

    TPolynomial operator+(const TPolynomial& other) const;
    TPolynomial operator-(const TPolynomial& other) const;
    TPolynomial operator*(const TPolynomial& other) const;

    std::vector<monomial> supp() const;

    bool operator==(const TPolynomial& other) const;
    bool operator!=(const TPolynomial& other) const;

    T evaluate(const point_type& point) const;

    std::optional<int> homogeneous_degree() const;

    TPolynomial homogeneous_component(int degree) const;

    std::string to_string() const;

private:
    void check_compat(const TPolynomial& other) const;

    ring_ptr m_ring;
    TTrie<key_type, T> m_terms;
};

template <typename T>
TPolynomial<T>::TPolynomial(ring_ptr ring)
    : m_ring(std::move(ring)) {
    if (m_ring == nullptr) {
        throw std::invalid_argument("TPolynomial: ring must not be null");
    }
}

template <typename T>
typename TPolynomial<T>::ring_ptr TPolynomial<T>::ring() const noexcept {
    return m_ring;
}

template <typename T>
typename TPolynomial<T>::size_type TPolynomial<T>::n_vars() const noexcept {
    return m_ring->n_vars();
}

template <typename T>
void TPolynomial<T>::set(const monomial& m, const T& coeff) {
    if (m.n_vars() != n_vars()) {
        throw std::invalid_argument("TPolynomial::set: n_vars mismatch");
    }
    if (coeff == T{}) {
        m_terms.erase(m.exponents());
    } else {
        m_terms.insert(m.exponents(), coeff);
    }
}

template <typename T>
const T* TPolynomial<T>::get(const monomial& m) const {
    if (m.n_vars() != n_vars()) {
        throw std::invalid_argument("TPolynomial::get: n_vars mismatch");
    }
    auto it = m_terms.find(m.exponents());
    if (it == m_terms.cend()) {
        return nullptr;
    }
    return &it.value();
}

template <typename T>
bool TPolynomial<T>::is_zero() const noexcept {
    return m_terms.empty();
}

template <typename T>
void TPolynomial<T>::check_compat(const TPolynomial& other) const {
    if (m_ring != other.m_ring) {
        throw std::invalid_argument("TPolynomial: ring mismatch");
    }
}

template <typename T>
TPolynomial<T>& TPolynomial<T>::operator+=(const TPolynomial& other) {
    check_compat(other);
    for (const auto& elem : other.m_terms) {
        key_type key = elem.key;
        T val = m_terms.contains(key)
            ? m_terms.at(key) + elem.value
            : elem.value;
        if (val == T{}) {
            m_terms.erase(key);
        } else {
            m_terms.insert(key, std::move(val));
        }
    }
    return *this;
}

template <typename T>
TPolynomial<T>& TPolynomial<T>::operator-=(const TPolynomial& other) {
    check_compat(other);
    for (const auto& elem : other.m_terms) {
        key_type key = elem.key;
        T val = m_terms.contains(key)
            ? m_terms.at(key) - elem.value
            : T{} - elem.value;
        if (val == T{}) {
            m_terms.erase(key);
        } else {
            m_terms.insert(key, std::move(val));
        }
    }
    return *this;
}

template <typename T>
TPolynomial<T>& TPolynomial<T>::operator*=(const TPolynomial& other) {
    check_compat(other);
    TPolynomial result(m_ring);
    for (const auto& a : m_terms) {
        for (const auto& b : other.m_terms) {
            key_type key(n_vars());
            for (size_type i = 0; i < n_vars(); ++i) {
                key[i] = a.key[i] + b.key[i];
            }
            T val = a.value * b.value;
            if (result.m_terms.contains(key)) {
                val = result.m_terms.at(key) + val;
            }
            if (val == T{}) {
                result.m_terms.erase(key);
            } else {
                result.m_terms.insert(key, std::move(val));
            }
        }
    }
    *this = std::move(result);
    return *this;
}

template <typename T>
TPolynomial<T> TPolynomial<T>::operator+(const TPolynomial& other) const {
    TPolynomial result = *this;
    result += other;
    return result;
}

template <typename T>
TPolynomial<T> TPolynomial<T>::operator-(const TPolynomial& other) const {
    TPolynomial result = *this;
    result -= other;
    return result;
}

template <typename T>
TPolynomial<T> TPolynomial<T>::operator*(const TPolynomial& other) const {
    TPolynomial result = *this;
    result *= other;
    return result;
}

template <typename T>
std::vector<Monomial> TPolynomial<T>::supp() const {
    std::vector<monomial> result;
    for (const auto& elem : m_terms) {
        result.emplace_back(elem.key);
    }
    return result;
}

template <typename T>
bool TPolynomial<T>::operator==(const TPolynomial& other) const {
    check_compat(other);
    if (m_terms.size() != other.m_terms.size()) {
        return false;
    }
    for (const auto& elem : m_terms) {
        auto it = other.m_terms.find(elem.key);
        if (it == other.m_terms.cend()) {
            return false;
        }
        if (it.value() != elem.value) {
            return false;
        }
    }
    return true;
}

template <typename T>
bool TPolynomial<T>::operator!=(const TPolynomial& other) const {
    return !(*this == other);
}

template <typename T>
T TPolynomial<T>::evaluate(const point_type& point) const {
    if (point.size() != n_vars()) {
        throw std::invalid_argument("TPolynomial::evaluate: point size mismatch");
    }
    T result{};
    for (const auto& elem : m_terms) {
        T term = elem.value;
        for (size_type i = 0; i < n_vars(); ++i) {
            for (int p = 0; p < elem.key[i]; ++p) {
                term = term * point[i];
            }
        }
        result = result + term;
    }
    return result;
}

template <typename T>
std::optional<int> TPolynomial<T>::homogeneous_degree() const {
    if (m_terms.empty()) {
        return 0;
    }
    std::optional<int> degree;
    for (const auto& elem : m_terms) {
        int d = 0;
        for (size_type i = 0; i < n_vars(); ++i) {
            d += elem.key[i];
        }
        if (!degree.has_value()) {
            degree = d;
        } else if (degree.value() != d) {
            return std::nullopt;
        }
    }
    return degree;
}

template <typename T>
TPolynomial<T> TPolynomial<T>::homogeneous_component(int degree) const {
    TPolynomial result(m_ring);
    for (const auto& elem : m_terms) {
        int d = 0;
        for (size_type i = 0; i < n_vars(); ++i) {
            d += elem.key[i];
        }
        if (d == degree) {
            result.m_terms.insert(elem.key, elem.value);
        }
    }
    return result;
}

template <typename T>
std::string TPolynomial<T>::to_string() const {
    if (m_terms.empty()) {
        return "0";
    }
    std::string result;
    bool first = true;
    for (const auto& elem : m_terms) {
        bool has_vars = false;
        for (size_type i = 0; i < n_vars(); ++i) {
            if (elem.key[i] != 0) {
                has_vars = true;
                break;
            }
        }
        std::string coeff_str = elem.value.get_str();
        bool negative = coeff_str[0] == '-';
        if (!first) {
            result += negative ? " - " : " + ";
        } else if (negative) {
            result += "-";
        }
        first = false;
        std::string abs_coeff_str = negative ? coeff_str.substr(1) : coeff_str;
        if (!has_vars || abs_coeff_str != "1") {
            result += abs_coeff_str;
        }
        for (size_type i = 0; i < n_vars(); ++i) {
            if (elem.key[i] != 0) {
                result += m_ring->var_name(i);
                if (elem.key[i] != 1) {
                    result += "^" + std::to_string(elem.key[i]);
                }
            }
        }
    }
    return result;
}

#endif //RGU_LABS_TERM4_ARITHMETIC_POLYNOMIAL_HPP
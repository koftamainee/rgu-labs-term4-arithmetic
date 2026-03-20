#pragma once

#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

#include "latex_serializable.hpp"

class Monomial : public ILatexSerializable {
public:
    using exponent_type  = int;
    using container      = std::vector<exponent_type>;
    using value_type     = exponent_type;
    using size_type      = std::size_t;
    using iterator       = container::iterator;
    using const_iterator = container::const_iterator;

    explicit Monomial(size_type n_vars) : m_exponents(n_vars, 0) {}

    explicit Monomial(container exponents) : m_exponents(std::move(exponents)) {
        for (const auto& e : m_exponents) {
            if (e < 0) {
                throw std::invalid_argument("Monomial: exponents must be non-negative");
            }
        }
    }

    size_type n_vars() const noexcept { return m_exponents.size(); }

    exponent_type& operator[](size_type i) { return m_exponents[i]; }
    const exponent_type& operator[](size_type i) const { return m_exponents[i]; }

    iterator begin() { return m_exponents.begin(); }
    iterator end() { return m_exponents.end(); }
    const_iterator begin() const { return m_exponents.begin(); }
    const_iterator end() const { return m_exponents.end(); }
    const_iterator cbegin() const { return m_exponents.cbegin(); }
    const_iterator cend() const { return m_exponents.cend(); }

    const container& exponents() const noexcept { return m_exponents; }

    exponent_type total_degree() const noexcept {
        return std::accumulate(m_exponents.begin(), m_exponents.end(), 0);
    }

    Monomial& operator*=(const Monomial& other) {
        if (m_exponents.size() != other.m_exponents.size()) {
            throw std::invalid_argument("Monomial::operator*=: n_vars mismatch");
        }
        for (size_type i = 0; i < m_exponents.size(); ++i) {
            m_exponents[i] += other.m_exponents[i];
        }
        return *this;
    }

    Monomial operator*(const Monomial& other) const {
        Monomial result = *this;
        result *= other;
        return result;
    }

    bool divides(const Monomial& other) const {
        if (m_exponents.size() != other.m_exponents.size()) {
            throw std::invalid_argument("Monomial::divides: n_vars mismatch");
        }
        for (size_type i = 0; i < m_exponents.size(); ++i) {
            if (m_exponents[i] > other.m_exponents[i]) {
                return false;
            }
        }
        return true;
    }

    bool operator==(const Monomial& other) const noexcept {
        return m_exponents == other.m_exponents;
    }

    bool operator!=(const Monomial& other) const noexcept {
        return !(*this == other);
    }

    bool operator<(const Monomial& other) const noexcept {
        return m_exponents < other.m_exponents;
    }

    bool operator<=(const Monomial& other) const noexcept {
        return !(other < *this);
    }

    bool operator>(const Monomial& other) const noexcept {
        return other < *this;
    }

    bool operator>=(const Monomial& other) const noexcept {
        return !(*this < other);
    }

    std::string to_latex() const override {
        std::string result;
        for (size_type i = 0; i < m_exponents.size(); ++i) {
            if (m_exponents[i] == 0) {
                continue;
            }
            result += "x_{" + std::to_string(i + 1) + "}";
            if (m_exponents[i] != 1) {
                result += "^{" + std::to_string(m_exponents[i]) + "}";
            }
        }
        if (result.empty()) {
            return "1";
        }
        return result;
    }

private:
    container m_exponents;
};
#pragma once

#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

#include "trie.hpp"
#include "monomial.hpp"
#include "poly_ring.hpp"

template <typename T>
class Polynomial {
public:
  using coeff_type = T;
  using ring_ptr = std::shared_ptr<const PolyRing<T>>;
  using key_type = Monomial::container;
  using size_type = std::size_t;
  using point_type = std::vector<T>;

  explicit Polynomial(ring_ptr ring);

  ring_ptr ring() const noexcept;
  size_type n_vars() const noexcept;

  void set(const Monomial& m, const T& coeff);
  const T* get(const Monomial& m) const;
  bool is_zero() const noexcept;

  Polynomial& operator+=(const Polynomial& other);
  Polynomial& operator-=(const Polynomial& other);
  Polynomial& operator*=(const Polynomial& other);

  Polynomial operator+(const Polynomial& other) const;
  Polynomial operator-(const Polynomial& other) const;
  Polynomial operator*(const Polynomial& other) const;

  bool operator==(const Polynomial& other) const;
  bool operator!=(const Polynomial& other) const;

  std::vector<Monomial> supp() const;

  T evaluate(const point_type& point) const;

  std::optional<int> homogeneous_degree() const;
  Polynomial homogeneous_component(int degree) const;

  const Trie<key_type, T>& terms() const noexcept;

private:
  void check_compat(const Polynomial& other) const;
  static int total_degree(const key_type& key);
  static T horner_impl(const typename Trie<key_type, T>::Node* node,
                       const point_type& point,
                       size_type var,
                       size_type n_vars);

  ring_ptr m_ring;
  Trie<key_type, T> m_terms;
};

template <typename T>
Polynomial<T>::Polynomial(ring_ptr ring)
  : m_ring(std::move(ring)) {
  if (m_ring == nullptr) {
    throw std::invalid_argument("Polynomial: ring must not be null");
  }
}

template <typename T>
typename Polynomial<T>::ring_ptr Polynomial<T>::ring() const noexcept {
  return m_ring;
}

template <typename T>
typename Polynomial<T>::size_type Polynomial<T>::n_vars() const noexcept {
  return m_ring->n_vars();
}

template <typename T>
void Polynomial<T>::set(const Monomial& m, const T& coeff) {
  if (m.n_vars() != n_vars()) {
    throw std::invalid_argument("Polynomial::set: n_vars mismatch");
  }
  if (coeff == T{}) {
    m_terms.erase(m.exponents());
  }
  else {
    m_terms.insert(m.exponents(), coeff);
  }
}

template <typename T>
const T* Polynomial<T>::get(const Monomial& m) const {
  if (m.n_vars() != n_vars()) {
    throw std::invalid_argument("Polynomial::get: n_vars mismatch");
  }
  auto it = m_terms.find(m.exponents());
  if (it == m_terms.cend()) {
    return nullptr;
  }
  return &it.value();
}

template <typename T>
bool Polynomial<T>::is_zero() const noexcept {
  return m_terms.empty();
}

template <typename T>
Polynomial<T>& Polynomial<T>::operator+=(const Polynomial& other) {
  check_compat(other);
  for (const auto& elem : other.m_terms) {
    key_type key = elem.key;
    auto it = m_terms.find(key);
    T val = (it != m_terms.end())
              ? it.value() + elem.value
              : elem.value;
    if (val == T{}) {
      m_terms.erase(key);
    }
    else {
      m_terms.insert(key, std::move(val));
    }
  }
  return *this;
}

template <typename T>
Polynomial<T>& Polynomial<T>::operator-=(const Polynomial& other) {
  check_compat(other);
  for (const auto& elem : other.m_terms) {
    key_type key = elem.key;
    auto it = m_terms.find(key);
    T val = (it != m_terms.end())
              ? it.value() - elem.value
              : T{} - elem.value;
    if (val == T{}) {
      m_terms.erase(key);
    }
    else {
      m_terms.insert(key, std::move(val));
    }
  }
  return *this;
}

template <typename T>
Polynomial<T>& Polynomial<T>::operator*=(const Polynomial& other) {
  check_compat(other);
  Polynomial result(m_ring);
  for (const auto& a : m_terms) {
    for (const auto& b : other.m_terms) {
      key_type key(n_vars());
      for (size_type i = 0; i < n_vars(); ++i) {
        key[i] = a.key[i] + b.key[i];
      }
      T val = a.value * b.value;
      auto it = result.m_terms.find(key);
      if (it != result.m_terms.end()) {
        val = it.value() + val;
      }
      if (val == T{}) {
        result.m_terms.erase(key);
      }
      else {
        result.m_terms.insert(key, std::move(val));
      }
    }
  }
  *this = std::move(result);
  return *this;
}

template <typename T>
Polynomial<T> Polynomial<T>::operator+(const Polynomial& other) const {
  Polynomial result = *this;
  result += other;
  return result;
}

template <typename T>
Polynomial<T> Polynomial<T>::operator-(const Polynomial& other) const {
  Polynomial result = *this;
  result -= other;
  return result;
}

template <typename T>
Polynomial<T> Polynomial<T>::operator*(const Polynomial& other) const {
  Polynomial result = *this;
  result *= other;
  return result;
}

template <typename T>
bool Polynomial<T>::operator==(const Polynomial& other) const {
  check_compat(other);
  return m_terms == other.m_terms;
}

template <typename T>
bool Polynomial<T>::operator!=(const Polynomial& other) const {
  check_compat(other);
  return m_terms != other.m_terms;
}

template <typename T>
std::vector<Monomial> Polynomial<T>::supp() const {
  std::vector<Monomial> result;
  for (const auto& elem : m_terms) {
    result.emplace_back(elem.key);
  }
  return result;
}

template <typename T>
T Polynomial<T>::evaluate(const point_type& point) const {
  if (point.size() != n_vars()) {
    throw std::invalid_argument("Polynomial::evaluate: point size mismatch");
  }
  if (m_terms.empty()) {
    return T{};
  }
  return horner_impl(m_terms.root(), point, 0, n_vars());
}

template <typename T>
T Polynomial<T>::horner_impl(
  const typename Trie<key_type, T>::Node* node,
  const point_type& point,
  size_type var,
  size_type n_vars) {
  if (!node) return T{};
  if (var == n_vars) {
    return node->value.value_or(T{});
  }
  if (node->children.empty()) {
    return node->value.value_or(T{});
  }

  int prev_deg = node->children.rbegin()->first;
  T result{};

  for (auto it = node->children.rbegin(); it != node->children.rend(); ++it) {
    int deg = it->first;
    int diff = prev_deg - deg;
    for (int i = 0; i < diff; ++i) {
      result = result * point[var];
    }
    result = result + horner_impl(it->second.get(), point, var + 1, n_vars);
    prev_deg = deg;
  }

  for (int i = 0; i < prev_deg; ++i) {
    result = result * point[var];
  }

  if (node->value.has_value()) {
    result = result + node->value.value();
  }

  return result;
}

template <typename T>
std::optional<int> Polynomial<T>::homogeneous_degree() const {
  if (m_terms.empty()) return 0;
  std::optional<int> degree;
  for (const auto& elem : m_terms) {
    int d = total_degree(elem.key);
    if (!degree.has_value()) {
      degree = d;
    }
    else if (degree.value() != d) {
      return std::nullopt;
    }
  }
  return degree;
}

template <typename T>
Polynomial<T> Polynomial<T>::homogeneous_component(int degree) const {
  Polynomial result(m_ring);
  for (const auto& elem : m_terms) {
    if (total_degree(elem.key) == degree) {
      result.m_terms.insert(elem.key, elem.value);
    }
  }
  return result;
}

template <typename T>
const Trie<typename Polynomial<T>::key_type, T>& Polynomial<T>::terms() const noexcept {
  return m_terms;
}

template <typename T>
void Polynomial<T>::check_compat(const Polynomial& other) const {
  if (m_ring != other.m_ring) {
    throw std::invalid_argument("Polynomial::check_compat: ring mismatch");
  }
}

template <typename T>
int Polynomial<T>::total_degree(const key_type& key) {
  int d = 0;
  for (const auto& e : key) d += e;
  return d;
}

template <typename T>
std::string to_latex(const Polynomial<T>& p) {
  if (p.is_zero()) {
    return "0";
  }
  std::string result;
  bool first = true;
  for (const auto& elem : p.terms()) {
    bool has_vars = false;
    for (typename Polynomial<T>::size_type i = 0; i < p.n_vars(); ++i) {
      if (elem.key[i] != 0) {
        has_vars = true;
        break;
      }
    }
    std::string coeff_str = to_latex(elem.value);

    bool needs_parens = has_vars && coeff_str.find(' ') != std::string::npos;
    if (needs_parens) {
      coeff_str.insert(0, 1, '(');
      coeff_str += ')';
    }

    bool negative = !needs_parens && !coeff_str.empty() && coeff_str[0] == '-';
    if (!first) {
      result += negative ? " - " : " + ";
    }
    else if (negative) {
      result += "-";
    }
    first = false;
    std::string abs_coeff_str = negative ? coeff_str.substr(1) : coeff_str;
    if (!has_vars || abs_coeff_str != "1") {
      result += abs_coeff_str;
    }
    for (typename Polynomial<T>::size_type i = 0; i < p.n_vars(); ++i) {
      if (elem.key[i] != 0) {
        result += p.ring()->var_name(i);
        if (elem.key[i] != 1) {
          result += "^{" + std::to_string(elem.key[i]) + "}";
        }
      }
    }
  }
  return result;
}

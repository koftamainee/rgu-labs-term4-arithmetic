#pragma once
#include <gmpxx.h>
#include <sstream>
#include <stdexcept>
#include <string>

template <typename T, int D>
class OQD {
public:
  using value_type = T;

  OQD() : m_a(0), m_b(0) {}
  OQD(T a, T b) : m_a(std::move(a)), m_b(std::move(b)) {}
  explicit OQD(T a) : m_a(std::move(a)), m_b(0) {}

  const T& a() const { return m_a; }
  const T& b() const { return m_b; }

  OQD operator-() const { return OQD(-m_a, -m_b); }

  OQD& operator+=(const OQD& rhs) {
    m_a += rhs.m_a;
    m_b += rhs.m_b;
    return *this;
  }

  OQD& operator-=(const OQD& rhs) {
    m_a -= rhs.m_a;
    m_b -= rhs.m_b;
    return *this;
  }

  OQD& operator*=(const OQD& rhs) {
    T new_a = m_a * rhs.m_a + m_b * rhs.m_b * T(D);
    T new_b = m_a * rhs.m_b + m_b * rhs.m_a;
    m_a = std::move(new_a);
    m_b = std::move(new_b);
    return *this;
  }

  OQD& operator/=(const OQD& rhs) {
    T n = rhs.norm();
    if (n == T(0)) {
      throw std::domain_error("OQD::operator/=: norm is zero");
    }
    T new_a = (m_a * rhs.m_a - m_b * rhs.m_b * T(D)) / n;
    T new_b = (m_b * rhs.m_a - m_a * rhs.m_b) / n;
    m_a = std::move(new_a);
    m_b = std::move(new_b);
    return *this;
  }

  OQD operator+(const OQD& rhs) const {
    auto r = *this;
    r += rhs;
    return r;
  }

  OQD operator-(const OQD& rhs) const {
    auto r = *this;
    r -= rhs;
    return r;
  }

  OQD operator*(const OQD& rhs) const {
    auto r = *this;
    r *= rhs;
    return r;
  }

  OQD operator/(const OQD& rhs) const {
    auto r = *this;
    r /= rhs;
    return r;
  }

  bool operator==(const OQD& rhs) const { return m_a == rhs.m_a && m_b == rhs.m_b; }
  bool operator!=(const OQD& rhs) const { return !(*this == rhs); }

  T norm() const { return m_a * m_a - m_b * m_b * T(D); }

  OQD conjugate() const { return OQD(m_a, -m_b); }

private:
  T m_a;
  T m_b;
};

template <typename T, int D>
std::string to_latex(const OQD<T, D>& x) {
  std::ostringstream oss;
  oss << x.a();
  std::string sa = oss.str();
  oss.str("");
  oss << x.b();
  std::string sb = oss.str();

  if (x.b() == T(0)) {
    return sa;
  }
  if (x.a() == T(0)) {
    if (sb == "1") {
      return "\\sqrt{" + std::to_string(D) + "}";
    }
    if (sb == "-1") {
      return "-\\sqrt{" + std::to_string(D) + "}";
    }
    return sb + "\\sqrt{" + std::to_string(D) + "}";
  }
  if (sb == "1") {
    return sa + " + \\sqrt{" + std::to_string(D) + "}";
  }
  if (sb == "-1") {
    return sa + " - \\sqrt{" + std::to_string(D) + "}";
  }
  return sa + " + " + sb + "\\sqrt{" + std::to_string(D) + "}";
}

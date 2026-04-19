#pragma once
#include <format>
#include <stdexcept>
#include <string>

template <typename T, T N>
class Zmod {
public:
  using value_type = T;

  Zmod() : m_value(0) {}

  explicit Zmod(T v) : m_value(((v % N) + N) % N) {}

  Zmod& operator+=(const Zmod& rhs) {
    m_value = (m_value + rhs.m_value) % N;
    return *this;
  }

  Zmod& operator-=(const Zmod& rhs) {
    m_value = ((m_value - rhs.m_value) % N + N) % N;
    return *this;
  }

  Zmod& operator*=(const Zmod& rhs) {
    m_value = (m_value * rhs.m_value) % N;
    return *this;
  }

  Zmod& operator/=(const Zmod& rhs) {
    *this *= rhs.inverse();
    return *this;
  }

  Zmod operator+(const Zmod& rhs) const {
    auto r = *this;
    r += rhs;
    return r;
  }

  Zmod operator-(const Zmod& rhs) const {
    auto r = *this;
    r -= rhs;
    return r;
  }

  Zmod operator*(const Zmod& rhs) const {
    auto r = *this;
    r *= rhs;
    return r;
  }

  Zmod operator/(const Zmod& rhs) const {
    auto r = *this;
    r /= rhs;
    return r;
  }

  Zmod operator-() const { return Zmod((N - m_value) % N); }

  bool operator==(const Zmod& rhs) const { return m_value == rhs.m_value; }
  bool operator!=(const Zmod& rhs) const { return m_value != rhs.m_value; }

  const T& value() const { return m_value; }

  Zmod inverse() const {
    if (m_value == 0) {
      throw std::domain_error("Zmod::inverse: division by zero");
    }
    T g, x, y;
    egcd(m_value, N, g, x, y);
    if (g != 1) {
      throw std::domain_error("Zmod::inverse: element not invertible");
    }
    return Zmod(x);
  }

private:
  static void egcd(T a, T b, T& g, T& x, T& y) {
    if (b == 0) {
      g = a;
      x = 1;
      y = 0;
      return;
    }
    T g1, x1, y1;
    egcd(b, a % b, g1, x1, y1);
    g = g1;
    x = y1;
    y = x1 - (a / b) * y1;
  }

  T m_value;
};

template <typename T, T N>
std::string to_latex(const Zmod<T, N>& z) {
  return std::format("{}", z.value());
}
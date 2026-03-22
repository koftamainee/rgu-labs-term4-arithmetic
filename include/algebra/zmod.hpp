#pragma once
#include <format>
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


  bool operator==(const Zmod& rhs) const {
    return m_value == rhs.m_value;
  }


  bool operator!=(const Zmod& rhs) const {
    return m_value != rhs.m_value;
  }


  const T& value() const { return m_value; }

private:
  T m_value;
};

template<typename T, T N>
std::string to_latex(const Zmod<T, N>& z) {
  return std::format("{}", z.value());
}

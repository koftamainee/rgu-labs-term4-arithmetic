#pragma once

#include <cmath>
#include <cstddef>
#include <initializer_list>
#include <stdexcept>
#include <string>
#include <vector>

#include "latex_serializable.hpp"
#include "utils.hpp"

template <typename T>
class Vector : public ILatexSerializable {
public:
  using value_type = T;
  using size_type = std::size_t;
  using iterator = typename std::vector<T>::iterator;
  using const_iterator = typename std::vector<T>::const_iterator;

  Vector() = default;

  explicit Vector(size_type dimension) : m_components(dimension, T{}) {}

  explicit Vector(std::vector<T> components) : m_components(std::move(components)) {}

  Vector(std::initializer_list<T> init) : m_components(init) {}

  Vector(const Vector&) = default;
  Vector(Vector&&) noexcept = default;
  Vector& operator=(const Vector&) = default;
  Vector& operator=(Vector&&) = default;

  size_type dimension() const noexcept { return m_components.size(); }

  const std::vector<T>& components() const noexcept { return m_components; }

  T& operator[](size_type index) {
    if (index >= dimension()) {
      throw std::out_of_range("Vector::operator[]: index out of range");
    }
    return m_components[index];
  }

  const T& operator[](size_type index) const {
    if (index >= dimension()) {
      throw std::out_of_range("Vector::operator[]: index out of range");
    }
    return m_components[index];
  }

  iterator begin() { return m_components.begin(); }
  iterator end() { return m_components.end(); }
  const_iterator begin() const { return m_components.begin(); }
  const_iterator end() const { return m_components.end(); }
  const_iterator cbegin() const { return m_components.cbegin(); }
  const_iterator cend() const { return m_components.cend(); }

  Vector& operator+=(const Vector& other) {
    check_dimension(other.dimension(), "operator+=");
    for (size_type i = 0; i < dimension(); ++i) {
      m_components[i] += other.m_components[i];
    }
    return *this;
  }

  Vector& operator-=(const Vector& other) {
    check_dimension(other.dimension(), "operator-=");
    for (size_type i = 0; i < dimension(); ++i) {
      m_components[i] -= other.m_components[i];
    }
    return *this;
  }

  Vector& operator*=(const T& scalar) {
    for (T& c : m_components) {
      c *= scalar;
    }
    return *this;
  }

  Vector& operator/=(const T& scalar) {
    if (scalar == T{}) {
      throw std::domain_error("Vector::operator/=: division by zero");
    }
    for (T& c : m_components) {
      c /= scalar;
    }
    return *this;
  }

  Vector operator+() const { return *this; }

  Vector operator-() const {
    Vector copy = *this;
    copy *= T{-1};
    return copy;
  }

  Vector operator+(const Vector& other) const {
    Vector r = *this;
    r += other;
    return r;
  }

  Vector operator-(const Vector& other) const {
    Vector r = *this;
    r -= other;
    return r;
  }

  Vector operator*(const T& scalar) const {
    Vector r = *this;
    r *= scalar;
    return r;
  }

  Vector operator/(const T& scalar) const {
    Vector r = *this;
    r /= scalar;
    return r;
  }

  friend Vector operator*(const T& scalar, const Vector& v) { return v * scalar; }

  bool operator==(const Vector& other) const {
    if (dimension() != other.dimension()) {
      return false;
    }
    for (size_type i = 0; i < dimension(); ++i) {
      if (m_components[i] != other.m_components[i]) {
        return false;
      }
    }
    return true;
  }

  bool operator!=(const Vector& other) const { return !(*this == other); }

  T dot(const Vector& other) const {
    check_dimension(other.dimension(), "dot");
    T result{};
    for (size_type i = 0; i < dimension(); ++i) {
      result += m_components[i] * other.m_components[i];
    }
    return result;
  }

  T norm() const { return std::sqrt(dot(*this)); }

  Vector normalize() const {
    check_non_zero();
    return *this / norm();
  }

  Vector cross_3d(const Vector& other) const {
    check_dimension(3, "cross_3d");
    other.check_dimension(3, "cross_3d");
    const auto& a = m_components;
    const auto& b = other.m_components;
    return Vector{
      a[1] * b[2] - a[2] * b[1],
      a[2] * b[0] - a[0] * b[2],
      a[0] * b[1] - a[1] * b[0]
    };
  }

  Vector cross_7d(const Vector& other) const {
    check_dimension(7, "cross_7d");
    other.check_dimension(7, "cross_7d");
    const auto& a = m_components;
    const auto& b = other.m_components;
    return Vector{
      a[1] * b[3] - a[3] * b[1] + a[2] * b[6] - a[6] * b[2] + a[4] * b[5] - a[5] * b[4],
      a[2] * b[4] - a[4] * b[2] + a[3] * b[0] - a[0] * b[3] + a[5] * b[6] - a[6] * b[5],
      a[3] * b[5] - a[5] * b[3] + a[4] * b[1] - a[1] * b[4] + a[6] * b[0] - a[0] * b[6],
      a[4] * b[6] - a[6] * b[4] + a[5] * b[2] - a[2] * b[5] + a[0] * b[1] - a[1] * b[0],
      a[5] * b[0] - a[0] * b[5] + a[6] * b[3] - a[3] * b[6] + a[1] * b[2] - a[2] * b[1],
      a[6] * b[1] - a[1] * b[6] + a[0] * b[4] - a[4] * b[0] + a[2] * b[3] - a[3] * b[2],
      a[0] * b[2] - a[2] * b[0] + a[1] * b[5] - a[5] * b[1] + a[3] * b[4] - a[4] * b[3]
    };
  }

  static T triple_product_3d(const Vector& a, const Vector& b, const Vector& c) {
    return a.dot(b.cross_3d(c));
  }

  static T triple_product_7d(const Vector& a, const Vector& b, const Vector& c) {
    return a.dot(b.cross_7d(c));
  }

  static Vector zero(size_type dimension) {
    return Vector(dimension);
  }

  static Vector basis_vector(size_type dimension, size_type index) {
    if (index >= dimension) {
      throw std::out_of_range("Vector::basis_vector: index out of range");
    }
    Vector result(dimension);
    result[index] = T{1};
    return result;
  }

  bool is_zero() const {
    for (const T& c : m_components) {
      if (c != T{}) {
        return false;
      }
    }
    return true;
  }

  bool is_orthogonal_to(const Vector& other) const {
    return dot(other) == T{};
  }

  friend T angle_between(const Vector& a, const Vector& b) {
    a.check_dimension(b.dimension(), "angle_between");
    if (a.is_zero() || b.is_zero()) {
      throw std::domain_error("Vector::angle_between: angle with zero vector is undefined");
    }
    T cosine = a.dot(b) / (a.norm() * b.norm());
    cosine = std::fmax(T{-1}, std::fmin(T{1}, cosine));
    return std::acos(cosine);
  }

  friend bool are_orthogonal(const Vector& a, const Vector& b) {
    return a.is_orthogonal_to(b);
  }

  friend bool are_collinear(const Vector& a, const Vector& b) {
    if (a.dimension() != b.dimension()) {
      return false;
    }
    if (a.is_zero() || b.is_zero()) {
      return true;
    }
    for (size_type i = 0; i < a.dimension(); ++i) {
      for (size_type j = i + 1; j < a.dimension(); ++j) {
        if (a[i] * b[j] - a[j] * b[i] != T{}) {
          return false;
        }
      }
    }
    return true;
  }

  friend bool is_point_on_segment(const Vector& pt, const Vector& a, const Vector& b, T eps = T{1e-10}) {
    pt.check_dimension(a.dimension(), "is_point_on_segment");
    pt.check_dimension(b.dimension(), "is_point_on_segment");
    Vector ab = b - a;
    Vector ap = pt - a;
    for (size_type i = 0; i < ab.dimension(); ++i) {
      for (size_type j = i + 1; j < ab.dimension(); ++j) {
        if (std::fabs(ap[i] * ab[j] - ap[j] * ab[i]) > eps) {
          return false;
        }
      }
    }
    T ab2 = ab.dot(ab);
    if (ab2 == T{}) {
      return pt == a;
    }
    T t = ap.dot(ab) / ab2;
    return t >= -eps && t <= T{1} + eps;
  }

  friend T point_to_segment_distance(const Vector& pt, const Vector& a, const Vector& b) {
    pt.check_dimension(a.dimension(), "point_to_segment_distance");
    pt.check_dimension(b.dimension(), "point_to_segment_distance");
    Vector ab = b - a;
    Vector ap = pt - a;
    T ab2 = ab.dot(ab);
    if (ab2 == T{}) {
      return (pt - a).norm();
    }
    T t = std::fmax(T{}, std::fmin(T{1}, ap.dot(ab) / ab2));
    return (pt - (a + ab * t)).norm();
  }

  std::string to_string() const {
    if (m_components.empty()) {
      return "[]";
    }
    std::string result = "[";
    for (size_type i = 0; i < m_components.size(); ++i) {
      if (i != 0) {
        result += ", ";
      }
      result += to_str(m_components[i]);
    }
    result += "]";
    return result;
  }

  std::string to_latex() const override {
    if (m_components.empty()) {
      return "\\begin{pmatrix}\\end{pmatrix}";
    }
    std::string result = "\\begin{pmatrix}";
    for (size_type i = 0; i < m_components.size(); ++i) {
      result += to_str(m_components[i]);
      if (i + 1 < m_components.size()) {
        result += " \\\\ ";
      }
    }
    result += "\\end{pmatrix}";
    return result;
  }

private:
  void check_dimension(size_type expected, const std::string& op) const {
    if (dimension() != expected) {
      throw std::invalid_argument("Vector::" + op + ": dimension mismatch");
    }
  }

  void check_non_zero() const {
    if (is_zero()) {
      throw std::domain_error("Vector: operation on zero vector");
    }
  }

  std::vector<T> m_components;
};

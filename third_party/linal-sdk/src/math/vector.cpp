#include "vector.h"

#include <cmath>
#include <cstdio>
#include <sstream>
#include <stdexcept>
#include <string>

#include "bigfloat.h"

void Vector::check_dimension(size_t expected,
                             const std::string &operation) const {
  if (dimension() != expected) {
    throw std::invalid_argument(
        "Vector::" + operation + " - dimension mismatch: " +
        std::to_string(dimension()) + " != " + std::to_string(expected));
  }
}

void Vector::check_non_zero() const {
  if (is_zero()) {
    throw std::domain_error("Vector operation on zero vector");
  }
}

Vector::Vector(size_t dimension) : components_(dimension) {}

Vector::Vector(const std::vector<bigfloat> &components)
    : components_(components) {}

Vector::Vector(std::initializer_list<bigfloat> init) : components_(init) {}

size_t Vector::dimension() const noexcept { return components_.size(); }

const bigfloat &Vector::operator[](size_t index) const {
  if (index >= dimension()) {
    throw std::out_of_range("Vector index out of range");
  }
  return components_[index];
}

bigfloat &Vector::operator[](size_t index) {
  if (index >= dimension()) {
    throw std::out_of_range("Vector index out of range");
  }
  return components_[index];
}

Vector &Vector::operator+=(const Vector &other) & {
  check_dimension(other.dimension(), "operator+=");

  for (size_t i = 0; i < dimension(); ++i) {
    components_[i] += other.components_[i];
  }
  return *this;
}

Vector &Vector::operator-=(const Vector &other) & {
  check_dimension(other.dimension(), "operator-=");

  for (size_t i = 0; i < dimension(); ++i) {
    components_[i] -= other.components_[i];
  }
  return *this;
}

Vector &Vector::operator*=(const bigfloat &scalar) & {
  for (auto &c : components_) {
    c *= scalar;
  }
  return *this;
}

Vector &Vector::operator/=(const bigfloat &scalar) & {
  if (scalar == 0) {
    throw std::domain_error("Division by zero");
  }
  for (auto &component : components_) {
    component /= scalar;
  }
  return *this;
}

Vector Vector::operator+() const { return *this; }

Vector Vector::operator-() const {
  auto copy = *this;
  return copy *= -1;
}

Vector operator+(Vector first, const Vector &second) { return first += second; }
Vector operator-(Vector first, const Vector &second) { return first -= second; }
Vector operator*(Vector vec, const bigfloat &scalar) { return vec *= scalar; }
Vector operator*(const bigfloat &scalar, Vector vec) { return vec *= scalar; }
Vector operator/(Vector vec, const bigfloat &scalar) { return vec /= scalar; }

bool operator==(const Vector &first, const Vector &second) {
  if (first.dimension() != second.dimension()) {
    return false;
  }
  for (size_t i = 0; i < first.dimension(); ++i) {
    if (first[i] != second[i]) {
      return false;
    }
  }
  return true;
}

bool operator!=(const Vector &first, const Vector &second) {
  return !(first == second);
}

bigfloat Vector::dot(const Vector &other) const {
  check_dimension(other.dimension(), "dot");

  bigfloat result = 0;
  for (size_t i = 0; i < dimension(); ++i) {
    result += components_[i] * other.components_[i];
  }
  return result;
}

bigfloat Vector::norm() const { return sqrt(dot(*this)); }

Vector Vector::normalize(const bigfloat &EPS) const {
  check_non_zero();
  const bigfloat n = norm();
  if (n < EPS) {
    throw std::domain_error("Cannot normalize zero vector");
  }
  return *this / n;
}

Vector Vector::cross_3d(const Vector &other) const {
  check_dimension(3, "cross_3d");
  other.check_dimension(3, "cross_3d");

  bigfloat x = components_[1] * other.components_[2] -
               components_[2] * other.components_[1];
  bigfloat y = components_[2] * other.components_[0] -
               components_[0] * other.components_[2];
  bigfloat z = components_[0] * other.components_[1] -
               components_[1] * other.components_[0];

  return Vector{x, y, z};
}

Vector Vector::cross_7d(const Vector &other) const {
  check_dimension(7, "cross_7d");
  other.check_dimension(7, "cross_7d");

  std::vector<bigfloat> r(7);

  r[0] = components_[1] * other.components_[3] -
         components_[3] * other.components_[1] +
         components_[2] * other.components_[6] -
         components_[6] * other.components_[2] +
         components_[4] * other.components_[5] -
         components_[5] * other.components_[4];

  r[1] = components_[2] * other.components_[4] -
         components_[4] * other.components_[2] +
         components_[3] * other.components_[0] -
         components_[0] * other.components_[3] +
         components_[5] * other.components_[6] -
         components_[6] * other.components_[5];

  r[2] = components_[3] * other.components_[5] -
         components_[5] * other.components_[3] +
         components_[4] * other.components_[1] -
         components_[1] * other.components_[4] +
         components_[6] * other.components_[0] -
         components_[0] * other.components_[6];

  r[3] = components_[4] * other.components_[6] -
         components_[6] * other.components_[4] +
         components_[5] * other.components_[2] -
         components_[2] * other.components_[5] +
         components_[0] * other.components_[1] -
         components_[1] * other.components_[0];

  r[4] = components_[5] * other.components_[0] -
         components_[0] * other.components_[5] +
         components_[6] * other.components_[3] -
         components_[3] * other.components_[6] +
         components_[1] * other.components_[2] -
         components_[2] * other.components_[1];

  r[5] = components_[6] * other.components_[1] -
         components_[1] * other.components_[6] +
         components_[0] * other.components_[4] -
         components_[4] * other.components_[0] +
         components_[2] * other.components_[3] -
         components_[3] * other.components_[2];

  r[6] = components_[0] * other.components_[2] -
         components_[2] * other.components_[0] +
         components_[1] * other.components_[5] -
         components_[5] * other.components_[1] +
         components_[3] * other.components_[4] -
         components_[4] * other.components_[3];

  return Vector{r[0], r[1], r[2], r[3], r[4], r[5], r[6]};
}

// Static methods
bigfloat Vector::triple_product_3d(const Vector &a, const Vector &b,
                                   const Vector &c) {
  return a.dot(b.cross_3d(c));
}

bigfloat Vector::triple_product_7d(const Vector &a, const Vector &b,
                                   const Vector &c) {
  return a.dot(b.cross_7d(c));
}

bool Vector::is_zero() const {
  for (const auto &component : components_) {
    if (component != 0) {
      return false;
    }
  }
  return true;
}

bool Vector::is_orthogonal_to(const Vector &other) const {
  return dot(other) == 0;
}

Vector Vector::zero(size_t dimension) { return Vector(dimension); }

Vector Vector::basis_vector(size_t dimension, size_t index) {
  if (index >= dimension) {
    throw std::out_of_range("Basis vector index out of range");
  }
  Vector result(dimension);
  result[index] = 1;
  return result;
}

bigfloat angle_between(const Vector &a, const Vector &b, const bigfloat &EPS) {
  a.check_dimension(b.dimension(), "angle_between");

  if (a.is_zero() || b.is_zero()) {
    throw std::domain_error("Angle with zero vector is undefined");
  }

  const bigfloat dot_product = a.dot(b);
  const bigfloat norms_product = a.norm() * b.norm();

  if (norms_product < EPS) {
    throw std::domain_error("Vectors too small for angle calculation");
  }

  return arccos(dot_product / norms_product, EPS * 1000);
}

bool are_orthogonal(const Vector &a, const Vector &b) {
  return a.is_orthogonal_to(b);
}

bool are_collinear(const Vector &a, const Vector &b) {
  if (a.dimension() != b.dimension()) {
    return false;
  }
  if (a.is_zero() || b.is_zero()) {
    return true;
  }

  const bigfloat ratio = a[0] / b[0];
  for (size_t i = 1; i < a.dimension(); ++i) {
    if (a[i] / b[i] != ratio) {
      return false;
    }
  }
  return true;
}

std::string Vector::to_string() const {
  if (components_.empty()) {
    return "[]";
  }

  std::string result = "[";
  for (size_t i = 0; i < components_.size(); ++i) {
    if (i != 0) {
      result += ", ";
    }
    result += components_[i].to_decimal();
  }
  result += "]";
  return result;
}

Vector::Vector(std::string const &str) {
  std::string cleaned = str;

  size_t start = cleaned.find('(');
  size_t end = cleaned.rfind(')');

  if (start != std::string::npos && end != std::string::npos && start < end) {
    cleaned = cleaned.substr(start + 1, end - start - 1);
  }

  std::stringstream in(cleaned);
  std::string num;
  while (in >> num) {
    components_.emplace_back(num);
  }
}

const std::vector<bigfloat> &Vector::components() const { return components_; }

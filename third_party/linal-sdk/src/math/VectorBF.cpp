#include "VectorBF.h"

#include <cmath>
#include <sstream>
#include <stdexcept>
#include <string>

void VectorBF::check_dimension(size_t expected,
                             const std::string &operation) const {
  if (dimension() != expected) {
    throw std::invalid_argument(
        "Vector::" + operation + " - dimension mismatch: " +
        std::to_string(dimension()) + " != " + std::to_string(expected));
  }
}

void VectorBF::check_non_zero() const {
  if (is_zero()) {
    throw std::domain_error("Vector operation on zero vector");
  }
}

VectorBF::VectorBF(size_t dimension) : components_(dimension) {}

VectorBF::VectorBF(const std::vector<bigfloat> &components)
    : components_(components) {}

VectorBF::VectorBF(std::initializer_list<bigfloat> init) : components_(init) {}

size_t VectorBF::dimension() const noexcept { return components_.size(); }

const bigfloat &VectorBF::operator[](size_t index) const {
  if (index >= dimension()) {
    throw std::out_of_range("Vector index out of range");
  }
  return components_[index];
}

bigfloat &VectorBF::operator[](size_t index) {
  if (index >= dimension()) {
    throw std::out_of_range("Vector index out of range");
  }
  return components_[index];
}

VectorBF &VectorBF::operator+=(const VectorBF &other) & {
  check_dimension(other.dimension(), "operator+=");

  for (size_t i = 0; i < dimension(); ++i) {
    components_[i] += other.components_[i];
  }
  return *this;
}

VectorBF &VectorBF::operator-=(const VectorBF &other) & {
  check_dimension(other.dimension(), "operator-=");

  for (size_t i = 0; i < dimension(); ++i) {
    components_[i] -= other.components_[i];
  }
  return *this;
}

VectorBF &VectorBF::operator*=(const bigfloat &scalar) & {
  for (auto &c : components_) {
    c *= scalar;
  }
  return *this;
}

VectorBF &VectorBF::operator/=(const bigfloat &scalar) & {
  if (scalar == bigfloat(0UL)) {
    throw std::domain_error("Division by zero");
  }
  for (auto &component : components_) {
    component /= scalar;
  }
  return *this;
}

VectorBF VectorBF::operator+() const { return *this; }

VectorBF VectorBF::operator-() const {
  auto copy = *this;
  return copy *= bigfloat(-1UL);
}

VectorBF operator+(VectorBF first, const VectorBF &second) { return first += second; }
VectorBF operator-(VectorBF first, const VectorBF &second) { return first -= second; }
VectorBF operator*(VectorBF vec, const bigfloat &scalar) { return vec *= scalar; }
VectorBF operator*(const bigfloat &scalar, VectorBF vec) { return vec *= scalar; }
VectorBF operator/(VectorBF vec, const bigfloat &scalar) { return vec /= scalar; }

bool operator==(const VectorBF &first, const VectorBF &second) {
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

bool operator!=(const VectorBF &first, const VectorBF &second) {
  return !(first == second);
}

bigfloat VectorBF::dot(const VectorBF &other) const {
  check_dimension(other.dimension(), "dot");

  bigfloat result(0UL);
  for (size_t i = 0; i < dimension(); ++i) {
    result += components_[i] * other.components_[i];
  }
  return result;
}

bigfloat VectorBF::norm() const { return sqrt(dot(*this)); }

VectorBF VectorBF::normalize(const bigfloat &EPS) const {
  check_non_zero();
  const bigfloat n = norm();
  if (n < EPS) {
    throw std::domain_error("Cannot normalize zero vector");
  }
  return *this / n;
}

VectorBF VectorBF::cross_3d(const VectorBF &other) const {
  check_dimension(3, "cross_3d");
  other.check_dimension(3, "cross_3d");

  bigfloat x = components_[1] * other.components_[2] -
               components_[2] * other.components_[1];
  bigfloat y = components_[2] * other.components_[0] -
               components_[0] * other.components_[2];
  bigfloat z = components_[0] * other.components_[1] -
               components_[1] * other.components_[0];

  return VectorBF{x, y, z};
}

VectorBF VectorBF::cross_7d(const VectorBF &other) const {
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

  return VectorBF{r[0], r[1], r[2], r[3], r[4], r[5], r[6]};
}

// Static methods
bigfloat VectorBF::triple_product_3d(const VectorBF &a, const VectorBF &b,
                                   const VectorBF &c) {
  return a.dot(b.cross_3d(c));
}

bigfloat VectorBF::triple_product_7d(const VectorBF &a, const VectorBF &b,
                                   const VectorBF &c) {
  return a.dot(b.cross_7d(c));
}

bool VectorBF::is_zero() const {
  for (const auto &component : components_) {
    if (component != bigfloat(0UL)) {
      return false;
    }
  }
  return true;
}

bool VectorBF::is_orthogonal_to(const VectorBF &other) const {
  return dot(other) == bigfloat(0UL);
}

VectorBF VectorBF::zero(size_t dimension) { return VectorBF(dimension); }

VectorBF VectorBF::basis_vector(size_t dimension, size_t index) {
  if (index >= dimension) {
    throw std::out_of_range("Basis vector index out of range");
  }
  VectorBF result(dimension);
  result[index] = bigfloat(1UL);
  return result;
}

bigfloat angle_between(const VectorBF &a, const VectorBF &b, const bigfloat &EPS) {
  a.check_dimension(b.dimension(), "angle_between");

  if (a.is_zero() || b.is_zero()) {
    throw std::domain_error("Angle with zero vector is undefined");
  }

  const bigfloat dot_product = a.dot(b);
  const bigfloat norms_product = a.norm() * b.norm();

  if (norms_product < EPS) {
    throw std::domain_error("Vectors too small for angle calculation");
  }

  return arccos(dot_product / norms_product, EPS * bigfloat(1000UL));
}

bool are_orthogonal(const VectorBF &a, const VectorBF &b) {
  return a.is_orthogonal_to(b);
}

bool are_collinear(const VectorBF &a, const VectorBF &b) {
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

std::string VectorBF::to_string() const {
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


const std::vector<bigfloat> &VectorBF::components() const { return components_; }

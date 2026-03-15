#pragma once

#include <string>
#include <vector>

#include "bigmath.hpp"

class VectorBF {
private:
  std::vector<bigfloat> components_;

  void check_dimension(size_t expected, const std::string &operation) const;
  void check_non_zero() const;

public:
  VectorBF() = default;
  explicit VectorBF(size_t dimension);
  VectorBF(const std::vector<bigfloat> &components);
  VectorBF(std::initializer_list<bigfloat> init);

  VectorBF(const VectorBF &) = default;
  VectorBF(VectorBF &&) noexcept = default;
  VectorBF &operator=(const VectorBF &) = default;
  VectorBF &operator=(VectorBF &&) noexcept = default;
  ~VectorBF() = default;

  size_t dimension() const noexcept;
  const bigfloat &operator[](size_t index) const;
  bigfloat &operator[](size_t index);

  VectorBF &operator+=(const VectorBF &other) &;
  VectorBF &operator-=(const VectorBF &other) &;
  VectorBF &operator*=(const bigfloat &scalar) &;
  VectorBF &operator/=(const bigfloat &scalar) &;

  VectorBF operator+() const;
  VectorBF operator-() const;

  friend VectorBF operator+(VectorBF first, const VectorBF &second);
  friend VectorBF operator-(VectorBF first, const VectorBF &second);
  friend VectorBF operator*(VectorBF vec, const bigfloat &scalar);
  friend VectorBF operator*(const bigfloat &scalar, VectorBF vec);
  friend VectorBF operator/(VectorBF vec, const bigfloat &scalar);

  friend bool operator==(const VectorBF &first, const VectorBF &second);
  friend bool operator!=(const VectorBF &first, const VectorBF &second);

  const std::vector<bigfloat> &components() const;

  bigfloat dot(const VectorBF &other) const;
  bigfloat norm() const;

  VectorBF normalize(const bigfloat &EPS = bigfloat::DEFAULT_EPS) const;

  VectorBF cross_3d(const VectorBF &other) const;
  VectorBF cross_7d(const VectorBF &other) const;

  static bigfloat triple_product_3d(const VectorBF &a, const VectorBF &b,
                                    const VectorBF &c);
  static bigfloat triple_product_7d(const VectorBF &a, const VectorBF &b,
                                    const VectorBF &c);

  bool is_zero() const;
  bool is_orthogonal_to(const VectorBF &other) const;

  std::string to_string() const;

  static VectorBF zero(size_t dimension);
  static VectorBF basis_vector(size_t dimension, size_t index);

  friend bigfloat angle_between(const VectorBF &a, const VectorBF &b,
                                const bigfloat &EPS);
  friend bool are_orthogonal(const VectorBF &a, const VectorBF &b);
  friend bool are_collinear(const VectorBF &a, const VectorBF &b);
};

bigfloat angle_between(const VectorBF &a, const VectorBF &b,
                       const bigfloat &EPS = bigfloat::DEFAULT_EPS);
bool are_orthogonal(const VectorBF &a, const VectorBF &b);
bool are_collinear(const VectorBF &a, const VectorBF &b);

bool is_point_on_segment(const VectorBF &pt, const VectorBF &a, const VectorBF &b,
                         const bigfloat &EPS = bigfloat::DEFAULT_EPS);

bigfloat point_to_segment_distance(const VectorBF &pt, const VectorBF &a,
                                   const VectorBF &b);

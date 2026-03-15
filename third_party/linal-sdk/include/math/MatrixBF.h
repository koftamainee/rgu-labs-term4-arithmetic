#pragma once

#include <string>
#include <vector>

#include "bigmath/bigfloat.hpp"
#include "VectorBF.h"

class MatrixBF {
 private:
  std::vector<std::vector<bigfloat>> data_;
  size_t rows_{}, cols_{};

  void check_same_size(const MatrixBF& other, const std::string& op) const;
  void check_square(const std::string& op) const;

 public:
  MatrixBF() = default;
  MatrixBF(size_t rows, size_t cols);
  MatrixBF(std::vector<std::vector<bigfloat>> const& data);

  size_t rows() const noexcept;
  size_t cols() const noexcept;
  bigfloat& at(size_t row, size_t col);
  const bigfloat& at(size_t row, size_t col) const;

  MatrixBF& operator+=(const MatrixBF& other);
  MatrixBF& operator-=(const MatrixBF& other);
  MatrixBF& operator*=(const bigfloat& scalar);
  MatrixBF& operator*=(const MatrixBF& other);

  MatrixBF operator+(const MatrixBF& other) const;
  MatrixBF operator-(const MatrixBF& other) const;
  MatrixBF operator*(const bigfloat& scalar) const;
  MatrixBF operator*(const MatrixBF& other) const;

  bool operator==(const MatrixBF& other) const;
  bool operator!=(const MatrixBF& other) const;

  bigfloat determinant() const;
  MatrixBF inverse() const;
  MatrixBF transpose() const;
  std::vector<bigfloat> solve_gauss(std::vector<bigfloat> const& b) const;
  std::vector<bigfloat> solve_gauss_jordan(
      std::vector<bigfloat> const& b) const;

  // std::vector<bigfloat> eigenvalues(
  // bigfloat const& EPS = bigfloat::DEFAULT_EPS) const;
  // std::vector<Vector> eigenvectors(
  // bigfloat const& EPS = bigfloat::DEFAULT_EPS) const;

  size_t rank() const;


  static size_t span_dimension(
      const std::vector<std::vector<bigfloat>>& vectors);
  static bool is_in_span(const std::vector<std::vector<bigfloat>>& basis,
                         const std::vector<bigfloat>& vector);

  std::string to_string() const;
};

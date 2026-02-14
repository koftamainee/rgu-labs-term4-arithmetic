#include "matrix.h"

#include <exception>
#include <sstream>
#include <stdexcept>
#include <string>

#include "vector.h"

void Matrix::check_same_size(const Matrix& other, const std::string& op) const {
  if (rows_ != other.rows_ || cols_ != other.cols_) {
    throw std::runtime_error("Matrix size mismatch in operation: " + op);
  }
}

void Matrix::check_square(const std::string& op) const {
  if (rows_ != cols_) {
    throw std::runtime_error("Matrix must be square for operation: " + op);
  }
}

Matrix::Matrix(size_t rows, size_t cols)
    : data_(rows, std::vector<bigfloat>(cols, 0)), rows_(rows), cols_(cols) {}

Matrix::Matrix(const std::vector<std::vector<bigfloat>>& data)
    : data_(data),
      rows_(data.size()),
      cols_(data.empty() ? 0 : data[0].size()) {
  for (const auto& row : data) {
    if (row.size() != cols_) {
      throw std::runtime_error("Inconsistent row sizes in matrix");
    }
  }
}

Matrix::Matrix(const std::string& str) {
  std::string clean_str = str;

  if (!clean_str.empty() && clean_str.front() == '[') {
    clean_str.erase(0, 1);
  }
  if (!clean_str.empty() && clean_str.back() == ']') {
    clean_str.pop_back();
  }

  std::vector<std::vector<bigfloat>> result;
  std::istringstream stream(clean_str);
  std::string token;

  while (std::getline(stream, token, ')')) {
    size_t start = token.find('(');
    if (start == std::string::npos) {
      continue;
    }

    std::string row_str = token.substr(start + 1);
    std::istringstream row_stream(row_str);
    std::vector<bigfloat> row;
    std::string input;
    while (row_stream >> input) {
      row.emplace_back(input);
    }

    if (!row.empty()) {
      result.push_back(std::move(row));
    }
  }

  if (result.empty()) {
    throw std::runtime_error("Failed to parse matrix from string");
  }

  size_t cols = result[0].size();
  for (const auto& row : result) {
    if (row.size() != cols) {
      throw std::runtime_error("Inconsistent row sizes in parsed matrix");
    }
  }

  rows_ = result.size();
  cols_ = cols;
  data_ = std::move(result);
}

size_t Matrix::rows() const noexcept { return rows_; }
size_t Matrix::cols() const noexcept { return cols_; }

bigfloat& Matrix::at(size_t row, size_t col) { return data_.at(row).at(col); }
const bigfloat& Matrix::at(size_t row, size_t col) const {
  return data_.at(row).at(col);
}

Matrix& Matrix::operator+=(const Matrix& other) {
  check_same_size(other, "+=");
  for (size_t i = 0; i < rows_; ++i) {
    for (size_t j = 0; j < cols_; ++j) {
      data_[i][j] += other.data_[i][j];
    }
  }
  return *this;
}

Matrix Matrix::operator+(const Matrix& other) const {
  Matrix result = *this;
  result += other;
  return result;
}

Matrix& Matrix::operator-=(const Matrix& other) {
  check_same_size(other, "-=");
  for (size_t i = 0; i < rows_; ++i) {
    for (size_t j = 0; j < cols_; ++j) {
      data_[i][j] -= other.data_[i][j];
    }
  }
  return *this;
}

Matrix Matrix::operator-(const Matrix& other) const {
  Matrix result = *this;
  result -= other;
  return result;
}

Matrix& Matrix::operator*=(const bigfloat& scalar) {
  for (size_t i = 0; i < rows_; ++i) {
    for (size_t j = 0; j < cols_; ++j) {
      data_[i][j] *= scalar;
    }
  }
  return *this;
}

Matrix Matrix::operator*(const bigfloat& scalar) const {
  Matrix result = *this;
  result *= scalar;
  return result;
}

Matrix& Matrix::operator*=(const Matrix& other) {
  if (cols_ != other.rows_) {
    throw std::runtime_error("Matrix multiplication dimension mismatch");
  }

  Matrix result(rows_, other.cols_);

  for (size_t i = 0; i < rows_; ++i) {
    for (size_t j = 0; j < other.cols_; ++j) {
      for (size_t k = 0; k < cols_; ++k) {
        bigfloat product = data_[i][k] * other.data_[k][j];
        result.data_[i][j] += product;
      }
    }
  }

  *this = std::move(result);
  return *this;
}

Matrix Matrix::operator*(const Matrix& other) const {
  Matrix result = *this;
  result *= other;
  return result;
}
bool Matrix::operator==(const Matrix& other) const {
  return data_ == other.data_;
}

bool Matrix::operator!=(const Matrix& other) const { return !(*this == other); }

std::string Matrix::to_string() const {
  std::string result;
  for (const auto& row : data_) {
    result += "(";
    for (size_t j = 0; j < row.size(); ++j) {
      result += row[j].to_decimal();
      if (j + 1 < row.size()) {
        result += " ";
      }
    }
    result += ") ";
  }
  return result;
}

bigfloat Matrix::determinant() const {
  check_square("determinant");
  size_t n = rows_;
  Matrix temp = *this;
  bigfloat det = 1;

  for (size_t i = 0; i < n; ++i) {
    size_t pivot = i;
    while (pivot < n && temp.at(pivot, i) == 0) {
      ++pivot;
    }
    if (pivot == n) {
      return 0;
    }

    if (pivot != i) {
      std::swap(temp.data_[i], temp.data_[pivot]);
      det = -det;
    }

    det *= temp.at(i, i);

    for (size_t j = i + 1; j < n; ++j) {
      bigfloat factor = temp.at(j, i) / temp.at(i, i);

      for (size_t k = i; k < n; ++k) {
        bigfloat old_value = temp.at(j, k);
        temp.at(j, k) -= factor * temp.at(i, k);
      }
    }
  }

  return det;
}

Matrix Matrix::inverse() const {
  check_square("inverse");
  size_t n = rows_;
  Matrix a = *this;
  Matrix inv(n, n);
  for (size_t i = 0; i < n; ++i) {
    inv.at(i, i) = 1;
  }

  for (size_t i = 0; i < n; ++i) {
    size_t pivot = i;
    while (pivot < n && a.at(pivot, i) == 0) {
      ++pivot;
    }
    if (pivot == n) {
      throw std::runtime_error("Singular matrix");
    }

    if (pivot != i) {
      std::swap(a.data_[i], a.data_[pivot]);
      std::swap(inv.data_[i], inv.data_[pivot]);
    }

    bigfloat div = a.at(i, i);

    for (size_t j = 0; j < n; ++j) {
      a.at(i, j) /= div;
      inv.at(i, j) /= div;
    }

    for (size_t j = 0; j < n; ++j) {
      if (i == j) {
        continue;
      }

      bigfloat factor = a.at(j, i);

      for (size_t k = 0; k < n; ++k) {
        a.at(j, k) -= factor * a.at(i, k);
        inv.at(j, k) -= factor * inv.at(i, k);
      }
    }
  }

  return inv;
}

std::vector<bigfloat> Matrix::solve_gauss(
    std::vector<bigfloat> const& b) const {
  check_square("solve_gauss");
  size_t n = rows_;
  Matrix a = *this;
  std::vector<bigfloat> x = b;

  for (size_t i = 0; i < n; ++i) {
    size_t pivot = i;
    while (pivot < n && a.at(pivot, i) == 0) {
      ++pivot;
    }

    if (pivot == n) {
      throw std::runtime_error("No unique solution");
    }

    if (pivot != i) {
      std::swap(a.data_[i], a.data_[pivot]);
      std::swap(x[i], x[pivot]);
    }

    for (size_t j = i + 1; j < n; ++j) {
      bigfloat factor = a.at(j, i) / a.at(i, i);

      for (size_t k = i; k < n; ++k) {
        bigfloat old_value = a.at(j, k);
        a.at(j, k) -= factor * a.at(i, k);
      }

      bigfloat old_rhs = x[j];
      x[j] -= factor * x[i];
    }
  }

  std::vector<bigfloat> result(n);
  for (int i = static_cast<int>(n) - 1; i >= 0; --i) {
    bigfloat sum = x[i];
    for (size_t j = i + 1; j < n; ++j) {
      sum -= a.at(i, j) * result[j];
    }

    result[i] = sum / a.at(i, i);
  }

  return result;
}

std::vector<bigfloat> Matrix::solve_gauss_jordan(
    std::vector<bigfloat> const& b) const {
  check_square("solve_gauss_jordan");
  size_t n = rows_;
  Matrix a = *this;
  std::vector<bigfloat> x = b;

  for (size_t i = 0; i < n; ++i) {
    size_t pivot = i;
    while (pivot < n && a.at(pivot, i) == 0) {
      ++pivot;
    }

    if (pivot == n) {
      throw std::runtime_error("No unique solution");
    }

    if (pivot != i) {
      std::swap(a.data_[i], a.data_[pivot]);
      std::swap(x[i], x[pivot]);
    }

    bigfloat div = a.at(i, i);

    for (size_t j = 0; j < n; ++j) {
      bigfloat old_val = a.at(i, j);
      a.at(i, j) /= div;
    }

    bigfloat old_rhs = x[i];
    x[i] /= div;

    for (size_t j = 0; j < n; ++j) {
      if (j == i) {
        continue;
      }

      bigfloat factor = a.at(j, i);

      for (size_t k = 0; k < n; ++k) {
        bigfloat old_value = a.at(j, k);
        a.at(j, k) -= factor * a.at(i, k);
      }

      bigfloat old_rhs_j = x[j];
      x[j] -= factor * x[i];
    }
  }

  return x;
}

size_t Matrix::rank() const {
  Matrix temp = *this;
  size_t rank = 0;
  size_t m = rows_;
  size_t n = cols_;

  for (size_t col = 0, row = 0; col < n && row < m; ++col) {
    size_t sel = row;
    for (size_t i = row + 1; i < m; ++i) {
      if ((temp.at(i, col)).abs() > (temp.at(sel, col)).abs()) {
        sel = i;
      }
    }

    if (temp.at(sel, col) == 0) {
      continue;
    }

    if (sel != row) {
      std::swap(temp.data_[row], temp.data_[sel]);
    }

    for (size_t i = row + 1; i < m; ++i) {
      bigfloat factor = temp.at(i, col) / temp.at(row, col);

      for (size_t j = col; j < n; ++j) {
        bigfloat before = temp.at(i, j);
        temp.at(i, j) -= factor * temp.at(row, j);
      }
    }

    ++rank;
    ++row;
  }

  return rank;
}

// std::vector<bigfloat> Matrix::eigenvalues(bigfloat const& EPS) const {
//   check_square("eigenvalues");
//
//   const size_t n = rows_;
//   Matrix A = *this;
//   const size_t MAX_ITER = 100;
//
//   for (size_t iter = 0; iter < MAX_ITER; ++iter) {
//     std::vector<Vector> Q_vectors;
//     for (size_t i = 0; i < n; ++i) {
//       Vector v(n);
//       for (size_t j = 0; j < n; ++j) {
//         v[j] = A.at(j, i);
//       }
//       Q_vectors.push_back(v);
//     }
//
//     Q_vectors = Vector::gram_schmidt_process(Q_vectors, EPS);
//
//     Matrix Q(n, n);
//     for (size_t i = 0; i < n; ++i) {
//       for (size_t j = 0; j < n; ++j) {
//         Q.at(j, i) = Q_vectors[i][j];
//       }
//     }
//
//     Matrix R = Q.transpose() * A;
//     A = R * Q;
//
//     bigfloat off_diagonal = 0;
//     for (size_t i = 0; i < n; ++i) {
//       for (size_t j = 0; j < n; ++j) {
//         if (i != j) {
//           off_diagonal += A.at(i, j).abs();
//         }
//       }
//     }
//
//     if (off_diagonal < EPS) {
//       break;
//     }
//   }
//
//   std::vector<bigfloat> result(n);
//   for (size_t i = 0; i < n; ++i) {
//     result[i] = A.at(i, i);
//   }
//
//   return result;
// }
//
// std::vector<Vector> Matrix::eigenvectors(bigfloat const& EPS) const {
//   std::vector<bigfloat> eigvals = eigenvalues();
//   std::vector<Vector> eigvecs;
//   const size_t n = rows_;
//   const size_t MAX_ITER = 50;
//
//   for (const bigfloat& lambda : eigvals) {
//     Matrix shifted = *this;
//     for (size_t i = 0; i < n; ++i) {
//       shifted.at(i, i) -= lambda;
//     }
//
//     Vector x(n);
//     for (size_t i = 0; i < n; ++i) {
//       x[i] = bigfloat(1);
//     }
//
//     for (size_t iter = 0; iter < MAX_ITER; ++iter) {
//       std::vector<bigfloat> b(n);
//       for (size_t i = 0; i < n; ++i) {
//         b[i] = x[i];
//       }
//       b = shifted.solve_gauss_jordan(b);
//       Vector new_x(b);
//
//       new_x = new_x.normalize(EPS);
//
//       if ((new_x - x).norm() < EPS) {
//         break;
//       }
//       x = new_x;
//     }
//
//     eigvecs.push_back(x.normalize(EPS));
//   }
//
//   return eigvecs;
// }

size_t Matrix::span_dimension(
    const std::vector<std::vector<bigfloat>>& vectors) {
  Matrix m(vectors.size(), vectors[0].size());
  for (size_t i = 0; i < vectors.size(); ++i) {
    m.data_[i] = vectors[i];
  }

  return m.rank();
}

bool Matrix::is_in_span(const std::vector<std::vector<bigfloat>>& basis,
                        const std::vector<bigfloat>& vector) {
  Matrix m(basis.size(), basis[0].size());
  for (size_t i = 0; i < basis.size(); ++i) {
    m.data_[i] = basis[i];
  }

  try {
    m.solve_gauss(vector);

    return true;
  } catch (std::exception const& e) {
    return false;
  }
}

Matrix Matrix::transpose() const {
  Matrix result(cols_, rows_);

  for (size_t i = 0; i < rows_; ++i) {
    for (size_t j = 0; j < cols_; ++j) {
      result.at(j, i) = at(i, j);
    }
  }

  return result;
}

Vector Matrix::nullspace_vector(const bigfloat& EPS) const {
  size_t m = rows_;
  size_t n = cols_;

  std::vector<std::vector<bigfloat>> mat(m, std::vector<bigfloat>(n));
  for (size_t i = 0; i < m; ++i) {
    for (size_t j = 0; j < n; ++j) {
      mat[i][j] = data_[i][j];
    }
  }

  size_t rank = 0;
  std::vector<int> pivot_col(m, -1);
  for (size_t col = 0; col < n && rank < m; ++col) {
    size_t pivot_row = rank;
    for (size_t i = rank + 1; i < m; ++i) {
      if (mat[i][col].abs() > mat[pivot_row][col].abs()) {
        pivot_row = i;
      }
    }

    if (mat[pivot_row][col].abs() < EPS) {
      continue;
    }

    if (pivot_row != rank) {
      std::swap(mat[pivot_row], mat[rank]);
    }

    bigfloat pivot_val = mat[rank][col];
    for (size_t j = col; j < n; ++j) {
      mat[rank][j] /= pivot_val;
    }

    for (size_t i = 0; i < m; ++i) {
      if (i != rank && mat[i][col].abs() > EPS) {
        bigfloat factor = mat[i][col];
        for (size_t j = col; j < n; ++j) {
          mat[i][j] -= factor * mat[rank][j];
        }
      }
    }

    pivot_col[rank] = static_cast<int>(col);
    ++rank;
  }

  if (rank == n) {
    return Vector(n);
  }

  std::vector<bigfloat> nullvec(n, bigfloat(0));

  std::vector<bool> is_pivot_col(n, false);
  for (size_t i = 0; i < rank; ++i) {
    if (pivot_col[i] >= 0) {
      is_pivot_col[pivot_col[i]] = true;
    }
  }

  size_t free_col = 0;
  for (; free_col < n; ++free_col) {
    if (!is_pivot_col[free_col]) {
      break;
    }
  }
  nullvec[free_col] = bigfloat(1);

  for (int i = static_cast<int>(rank) - 1; i >= 0; --i) {
    int pc = pivot_col[i];
    if (pc == -1) {
      continue;
    }

    bigfloat sum = 0;
    for (size_t j = pc + 1; j < n; ++j) {
      sum += mat[i][j] * nullvec[j];
    }
    nullvec[pc] = -sum;
  }

  return nullvec;
}

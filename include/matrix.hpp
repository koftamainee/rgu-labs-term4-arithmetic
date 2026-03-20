#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <vector>

#include "utils.hpp"

template <typename T>
class Matrix{
public:
  using value_type = T;
  using size_type = std::size_t;

  Matrix() = default;

  Matrix(size_type rows, size_type cols)
    : m_data(rows, std::vector<T>(cols, T{})), m_rows(rows), m_cols(cols) {}

  explicit Matrix(std::vector<std::vector<T>> data)
    : m_data(std::move(data)),
      m_rows(m_data.size()),
      m_cols(m_data.empty() ? 0 : m_data[0].size()) {
    for (const auto& row : m_data) {
      if (row.size() != m_cols) {
        throw std::runtime_error("Matrix: inconsistent row sizes");
      }
    }
  }

  Matrix(const Matrix&) = default;
  Matrix(Matrix&&) noexcept = default;
  Matrix& operator=(const Matrix&) = default;
  Matrix& operator=(Matrix&&) = default;
  ~Matrix() = default;

  size_type rows() const noexcept { return m_rows; }
  size_type cols() const noexcept { return m_cols; }

  T& at(size_type row, size_type col) { return m_data.at(row).at(col); }
  const T& at(size_type row, size_type col) const { return m_data.at(row).at(col); }

  Matrix& operator+=(const Matrix& other) {
    check_same_size(other, "+=");
    for (size_type i = 0; i < m_rows; ++i) {
      for (size_type j = 0; j < m_cols; ++j) {
        m_data[i][j] += other.m_data[i][j];
      }
    }
    return *this;
  }

  Matrix& operator-=(const Matrix& other) {
    check_same_size(other, "-=");
    for (size_type i = 0; i < m_rows; ++i) {
      for (size_type j = 0; j < m_cols; ++j) {
        m_data[i][j] -= other.m_data[i][j];
      }
    }
    return *this;
  }

  Matrix& operator*=(const T& scalar) {
    for (size_type i = 0; i < m_rows; ++i) {
      for (size_type j = 0; j < m_cols; ++j) {
        m_data[i][j] *= scalar;
      }
    }
    return *this;
  }

  Matrix& operator*=(const Matrix& other) {
    if (m_cols != other.m_rows) {
      throw std::runtime_error("Matrix: multiplication dimension mismatch");
    }
    Matrix result(m_rows, other.m_cols);
    for (size_type i = 0; i < m_rows; ++i) {
      for (size_type j = 0; j < other.m_cols; ++j) {
        for (size_type k = 0; k < m_cols; ++k) {
          result.m_data[i][j] += m_data[i][k] * other.m_data[k][j];
        }
      }
    }
    *this = std::move(result);
    return *this;
  }

  Matrix operator+(const Matrix& other) const {
    Matrix r = *this;
    r += other;
    return r;
  }

  Matrix operator-(const Matrix& other) const {
    Matrix r = *this;
    r -= other;
    return r;
  }

  Matrix operator*(const T& scalar) const {
    Matrix r = *this;
    r *= scalar;
    return r;
  }

  Matrix operator*(const Matrix& other) const {
    Matrix r = *this;
    r *= other;
    return r;
  }

  friend Matrix operator*(const T& scalar, const Matrix& m) { return m * scalar; }

  bool operator==(const Matrix& other) const { return m_data == other.m_data; }
  bool operator!=(const Matrix& other) const { return !(*this == other); }

  Matrix transpose() const {
    Matrix result(m_cols, m_rows);
    for (size_type i = 0; i < m_rows; ++i) {
      for (size_type j = 0; j < m_cols; ++j) {
        result.at(j, i) = m_data[i][j];
      }
    }
    return result;
  }

  T determinant() const {
    check_square("determinant");
    size_type n = m_rows;
    Matrix temp = *this;
    T det = T{1};
    for (size_type i = 0; i < n; ++i) {
      size_type pivot = i;
      while (pivot < n && temp.m_data[pivot][i] == T{}) {
        ++pivot;
      }
      if (pivot == n) {
        return T{};
      }
      if (pivot != i) {
        std::swap(temp.m_data[i], temp.m_data[pivot]);
        det = det * T{-1};
      }
      det = det * temp.m_data[i][i];
      for (size_type j = i + 1; j < n; ++j) {
        T factor = temp.m_data[j][i] / temp.m_data[i][i];
        for (size_type k = i; k < n; ++k) {
          temp.m_data[j][k] -= factor * temp.m_data[i][k];
        }
      }
    }
    return det;
  }

  Matrix inverse() const {
    check_square("inverse");
    size_type n = m_rows;
    Matrix a = *this;
    Matrix inv(n, n);
    for (size_type i = 0; i < n; ++i) {
      inv.m_data[i][i] = T{1};
    }
    for (size_type i = 0; i < n; ++i) {
      size_type pivot = i;
      while (pivot < n && a.m_data[pivot][i] == T{}) {
        ++pivot;
      }
      if (pivot == n) {
        throw std::runtime_error("Matrix: singular matrix");
      }
      if (pivot != i) {
        std::swap(a.m_data[i], a.m_data[pivot]);
        std::swap(inv.m_data[i], inv.m_data[pivot]);
      }
      T div = a.m_data[i][i];
      for (size_type j = 0; j < n; ++j) {
        a.m_data[i][j] = a.m_data[i][j] / div;
        inv.m_data[i][j] = inv.m_data[i][j] / div;
      }
      for (size_type j = 0; j < n; ++j) {
        if (j == i) {
          continue;
        }
        T factor = a.m_data[j][i];
        for (size_type k = 0; k < n; ++k) {
          a.m_data[j][k] -= factor * a.m_data[i][k];
          inv.m_data[j][k] -= factor * inv.m_data[i][k];
        }
      }
    }
    return inv;
  }

  std::vector<T> solve_gauss(const std::vector<T>& b) const {
    check_square("solve_gauss");
    size_type n = m_rows;
    Matrix a = *this;
    std::vector<T> x = b;
    for (size_type i = 0; i < n; ++i) {
      size_type pivot = i;
      while (pivot < n && a.m_data[pivot][i] == T{}) {
        ++pivot;
      }
      if (pivot == n) {
        throw std::runtime_error("Matrix: no unique solution");
      }
      if (pivot != i) {
        std::swap(a.m_data[i], a.m_data[pivot]);
        std::swap(x[i], x[pivot]);
      }
      for (size_type j = i + 1; j < n; ++j) {
        T factor = a.m_data[j][i] / a.m_data[i][i];
        for (size_type k = i; k < n; ++k) {
          a.m_data[j][k] -= factor * a.m_data[i][k];
        }
        x[j] -= factor * x[i];
      }
    }
    std::vector<T> result(n);
    for (int i = static_cast<int>(n) - 1; i >= 0; --i) {
      T sum = x[i];
      for (size_type j = i + 1; j < n; ++j) {
        sum -= a.m_data[i][j] * result[j];
      }
      result[i] = sum / a.m_data[i][i];
    }
    return result;
  }

  std::vector<T> solve_gauss_jordan(const std::vector<T>& b) const {
    check_square("solve_gauss_jordan");
    size_type n = m_rows;
    Matrix a = *this;
    std::vector<T> x = b;
    for (size_type i = 0; i < n; ++i) {
      size_type pivot = i;
      while (pivot < n && a.m_data[pivot][i] == T{}) {
        ++pivot;
      }
      if (pivot == n) {
        throw std::runtime_error("Matrix: no unique solution");
      }
      if (pivot != i) {
        std::swap(a.m_data[i], a.m_data[pivot]);
        std::swap(x[i], x[pivot]);
      }
      T div = a.m_data[i][i];
      for (size_type j = 0; j < n; ++j) {
        a.m_data[i][j] = a.m_data[i][j] / div;
      }
      x[i] = x[i] / div;
      for (size_type j = 0; j < n; ++j) {
        if (j == i) {
          continue;
        }
        T factor = a.m_data[j][i];
        for (size_type k = 0; k < n; ++k) {
          a.m_data[j][k] -= factor * a.m_data[i][k];
        }
        x[j] -= factor * x[i];
      }
    }
    return x;
  }

  size_type rank() const {
    Matrix temp = *this;
    size_type rnk = 0;
    for (size_type col = 0, row = 0; col < m_cols && row < m_rows; ++col) {
      size_type sel = row;
      for (size_type i = row + 1; i < m_rows; ++i) {
        if (std::abs(temp.m_data[i][col]) > std::abs(temp.m_data[sel][col])) {
          sel = i;
        }
      }
      if (std::abs(temp.m_data[sel][col]) < T{1e-10}) {
        continue;
      }
      if (sel != row) {
        std::swap(temp.m_data[row], temp.m_data[sel]);
      }
      for (size_type i = row + 1; i < m_rows; ++i) {
        T factor = temp.m_data[i][col] / temp.m_data[row][col];
        for (size_type j = col; j < m_cols; ++j) {
          temp.m_data[i][j] -= factor * temp.m_data[row][j];
        }
      }
      ++rnk;
      ++row;
    }
    return rnk;
  }

  static size_type span_dimension(const std::vector<std::vector<T>>& vectors) {
    Matrix m(vectors.size(), vectors[0].size());
    for (size_type i = 0; i < vectors.size(); ++i) {
      m.m_data[i] = vectors[i];
    }
    return m.rank();
  }

  static bool is_in_span(const std::vector<std::vector<T>>& basis,
                         const std::vector<T>& vector) {
    Matrix m(basis.size(), basis[0].size());
    for (size_type i = 0; i < basis.size(); ++i) {
      m.m_data[i] = basis[i];
    }
    try {
      m.solve_gauss(vector);
      return true;
    }
    catch (...) {
      return false;
    }
  }

  std::string to_string() const {
    std::string result;
    for (const auto& row : m_data) {
      result += "(";
      for (size_type j = 0; j < row.size(); ++j) {
        result += to_str(row[j]);
        if (j + 1 < row.size()) {
          result += " ";
        }
      }
      result += ") ";
    }
    return result;
  }

  std::string to_latex() const {
    std::string result = "\\begin{pmatrix}";
    for (size_type i = 0; i < m_rows; ++i) {
      for (size_type j = 0; j < m_cols; ++j) {
        result += to_str(m_data[i][j]);
        if (j + 1 < m_cols) {
          result += " & ";
        }
      }
      if (i + 1 < m_rows) {
        result += " \\\\ ";
      }
    }
    result += "\\end{pmatrix}";
    return result;
  }

private:
  void check_same_size(const Matrix& other, const std::string& op) const {
    if (m_rows != other.m_rows || m_cols != other.m_cols) {
      throw std::runtime_error("Matrix: size mismatch in operation: " + op);
    }
  }

  void check_square(const std::string& op) const {
    if (m_rows != m_cols) {
      throw std::runtime_error("Matrix: must be square for operation: " + op);
    }
  }

  std::vector<std::vector<T>> m_data;
  size_type m_rows = 0;
  size_type m_cols = 0;
};

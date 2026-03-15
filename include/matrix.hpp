#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <vector>

class Matrix {
private:
    std::vector<std::vector<double>> data_;
    size_t rows_{};
    size_t cols_{};

    void check_same_size(const Matrix &other, const std::string &op) const {
        if (rows_ != other.rows_ || cols_ != other.cols_)
            throw std::runtime_error("Matrix size mismatch in operation: " + op);
    }

    void check_square(const std::string &op) const {
        if (rows_ != cols_)
            throw std::runtime_error("Matrix must be square for operation: " + op);
    }

public:
    Matrix() = default;

    Matrix(size_t rows, size_t cols)
        : data_(rows, std::vector<double>(cols, 0.0)), rows_(rows), cols_(cols) {}

    Matrix(const std::vector<std::vector<double>> &data)
        : data_(data),
          rows_(data.size()),
          cols_(data.empty() ? 0 : data[0].size()) {
        for (const auto &row : data_)
            if (row.size() != cols_)
                throw std::runtime_error("Inconsistent row sizes in matrix");
    }

    Matrix(const Matrix &)            = default;
    Matrix(Matrix &&) noexcept        = default;
    Matrix &operator=(const Matrix &) = default;
    Matrix &operator=(Matrix &&)      = default;
    ~Matrix()                         = default;

    size_t rows() const noexcept { return rows_; }
    size_t cols() const noexcept { return cols_; }

    double &at(size_t row, size_t col) { return data_.at(row).at(col); }
    const double &at(size_t row, size_t col) const { return data_.at(row).at(col); }

    Matrix &operator+=(const Matrix &other) {
        check_same_size(other, "+=");
        for (size_t i = 0; i < rows_; ++i)
            for (size_t j = 0; j < cols_; ++j)
                data_[i][j] += other.data_[i][j];
        return *this;
    }

    Matrix &operator-=(const Matrix &other) {
        check_same_size(other, "-=");
        for (size_t i = 0; i < rows_; ++i)
            for (size_t j = 0; j < cols_; ++j)
                data_[i][j] -= other.data_[i][j];
        return *this;
    }

    Matrix &operator*=(double scalar) {
        for (size_t i = 0; i < rows_; ++i)
            for (size_t j = 0; j < cols_; ++j)
                data_[i][j] *= scalar;
        return *this;
    }

    Matrix &operator*=(const Matrix &other) {
        if (cols_ != other.rows_)
            throw std::runtime_error("Matrix multiplication dimension mismatch");
        Matrix result(rows_, other.cols_);
        for (size_t i = 0; i < rows_; ++i)
            for (size_t j = 0; j < other.cols_; ++j)
                for (size_t k = 0; k < cols_; ++k)
                    result.data_[i][j] += data_[i][k] * other.data_[k][j];
        *this = std::move(result);
        return *this;
    }

    Matrix operator+(const Matrix &other) const { Matrix r = *this; r += other; return r; }
    Matrix operator-(const Matrix &other) const { Matrix r = *this; r -= other; return r; }
    Matrix operator*(double scalar)        const { Matrix r = *this; r *= scalar; return r; }
    Matrix operator*(const Matrix &other)  const { Matrix r = *this; r *= other;  return r; }

    friend Matrix operator*(double scalar, const Matrix &m) { return m * scalar; }

    bool operator==(const Matrix &other) const { return data_ == other.data_; }
    bool operator!=(const Matrix &other) const { return !(*this == other); }

    Matrix transpose() const {
        Matrix result(cols_, rows_);
        for (size_t i = 0; i < rows_; ++i)
            for (size_t j = 0; j < cols_; ++j)
                result.at(j, i) = data_[i][j];
        return result;
    }

    double determinant() const {
        check_square("determinant");
        size_t n = rows_;
        Matrix temp = *this;
        double det = 1.0;
        for (size_t i = 0; i < n; ++i) {
            size_t pivot = i;
            while (pivot < n && temp.data_[pivot][i] == 0.0) ++pivot;
            if (pivot == n) return 0.0;
            if (pivot != i) {
                std::swap(temp.data_[i], temp.data_[pivot]);
                det = -det;
            }
            det *= temp.data_[i][i];
            for (size_t j = i + 1; j < n; ++j) {
                double factor = temp.data_[j][i] / temp.data_[i][i];
                for (size_t k = i; k < n; ++k)
                    temp.data_[j][k] -= factor * temp.data_[i][k];
            }
        }
        return det;
    }

    Matrix inverse() const {
        check_square("inverse");
        size_t n = rows_;
        Matrix a = *this;
        Matrix inv(n, n);
        for (size_t i = 0; i < n; ++i) inv.data_[i][i] = 1.0;
        for (size_t i = 0; i < n; ++i) {
            size_t pivot = i;
            while (pivot < n && a.data_[pivot][i] == 0.0) ++pivot;
            if (pivot == n) throw std::runtime_error("Singular matrix");
            if (pivot != i) {
                std::swap(a.data_[i], a.data_[pivot]);
                std::swap(inv.data_[i], inv.data_[pivot]);
            }
            double div = a.data_[i][i];
            for (size_t j = 0; j < n; ++j) {
                a.data_[i][j]   /= div;
                inv.data_[i][j] /= div;
            }
            for (size_t j = 0; j < n; ++j) {
                if (j == i) continue;
                double factor = a.data_[j][i];
                for (size_t k = 0; k < n; ++k) {
                    a.data_[j][k]   -= factor * a.data_[i][k];
                    inv.data_[j][k] -= factor * inv.data_[i][k];
                }
            }
        }
        return inv;
    }

    std::vector<double> solve_gauss(const std::vector<double> &b) const {
        check_square("solve_gauss");
        size_t n = rows_;
        Matrix a = *this;
        std::vector<double> x = b;
        for (size_t i = 0; i < n; ++i) {
            size_t pivot = i;
            while (pivot < n && a.data_[pivot][i] == 0.0) ++pivot;
            if (pivot == n) throw std::runtime_error("No unique solution");
            if (pivot != i) {
                std::swap(a.data_[i], a.data_[pivot]);
                std::swap(x[i], x[pivot]);
            }
            for (size_t j = i + 1; j < n; ++j) {
                double factor = a.data_[j][i] / a.data_[i][i];
                for (size_t k = i; k < n; ++k)
                    a.data_[j][k] -= factor * a.data_[i][k];
                x[j] -= factor * x[i];
            }
        }
        std::vector<double> result(n);
        for (int i = static_cast<int>(n) - 1; i >= 0; --i) {
            double sum = x[i];
            for (size_t j = i + 1; j < n; ++j)
                sum -= a.data_[i][j] * result[j];
            result[i] = sum / a.data_[i][i];
        }
        return result;
    }

    std::vector<double> solve_gauss_jordan(const std::vector<double> &b) const {
        check_square("solve_gauss_jordan");
        size_t n = rows_;
        Matrix a = *this;
        std::vector<double> x = b;
        for (size_t i = 0; i < n; ++i) {
            size_t pivot = i;
            while (pivot < n && a.data_[pivot][i] == 0.0) ++pivot;
            if (pivot == n) throw std::runtime_error("No unique solution");
            if (pivot != i) {
                std::swap(a.data_[i], a.data_[pivot]);
                std::swap(x[i], x[pivot]);
            }
            double div = a.data_[i][i];
            for (size_t j = 0; j < n; ++j) a.data_[i][j] /= div;
            x[i] /= div;
            for (size_t j = 0; j < n; ++j) {
                if (j == i) continue;
                double factor = a.data_[j][i];
                for (size_t k = 0; k < n; ++k)
                    a.data_[j][k] -= factor * a.data_[i][k];
                x[j] -= factor * x[i];
            }
        }
        return x;
    }

    size_t rank() const {
        Matrix temp = *this;
        size_t rnk = 0;
        for (size_t col = 0, row = 0; col < cols_ && row < rows_; ++col) {
            size_t sel = row;
            for (size_t i = row + 1; i < rows_; ++i)
                if (std::fabs(temp.data_[i][col]) > std::fabs(temp.data_[sel][col]))
                    sel = i;
            if (std::fabs(temp.data_[sel][col]) < 1e-10) continue;
            if (sel != row) std::swap(temp.data_[row], temp.data_[sel]);
            for (size_t i = row + 1; i < rows_; ++i) {
                double factor = temp.data_[i][col] / temp.data_[row][col];
                for (size_t j = col; j < cols_; ++j)
                    temp.data_[i][j] -= factor * temp.data_[row][j];
            }
            ++rnk;
            ++row;
        }
        return rnk;
    }

    static size_t span_dimension(const std::vector<std::vector<double>> &vectors) {
        Matrix m(vectors.size(), vectors[0].size());
        for (size_t i = 0; i < vectors.size(); ++i)
            m.data_[i] = vectors[i];
        return m.rank();
    }

    static bool is_in_span(const std::vector<std::vector<double>> &basis,
                            const std::vector<double> &vector) {
        Matrix m(basis.size(), basis[0].size());
        for (size_t i = 0; i < basis.size(); ++i)
            m.data_[i] = basis[i];
        try {
            m.solve_gauss(vector);
            return true;
        } catch (...) {
            return false;
        }
    }

    std::string to_string() const {
        std::string result;
        for (const auto &row : data_) {
            result += "(";
            for (size_t j = 0; j < row.size(); ++j) {
                result += std::to_string(row[j]);
                if (j + 1 < row.size()) result += " ";
            }
            result += ") ";
        }
        return result;
    }
};
#pragma once

#include <cmath>
#include <cstddef>
#include <initializer_list>
#include <stdexcept>
#include <string>
#include <vector>

class Vector {
private:
    std::vector<double> components_;

    void check_dimension(size_t expected, const std::string &op) const {
        if (dimension() != expected)
            throw std::invalid_argument(
                "VectorDouble::" + op + " - dimension mismatch: " +
                std::to_string(dimension()) + " != " + std::to_string(expected));
    }

    void check_non_zero() const {
        if (is_zero())
            throw std::domain_error("VectorDouble: operation on zero vector");
    }

public:

    Vector() = default;

    explicit Vector(size_t dimension) : components_(dimension, 0.0) {}

    Vector(const std::vector<double> &components)
        : components_(components) {}

    Vector(std::initializer_list<double> init) : components_(init) {}

    Vector(const Vector &)            = default;
    Vector(Vector &&) noexcept        = default;
    Vector &operator=(const Vector &) = default;
    Vector &operator=(Vector &&)      = default;
    ~Vector()                               = default;

    size_t dimension() const noexcept { return components_.size(); }

    const std::vector<double> &components() const { return components_; }

    const double &operator[](size_t index) const {
        if (index >= dimension())
            throw std::out_of_range("VectorDouble: index out of range");
        return components_[index];
    }

    double &operator[](size_t index) {
        if (index >= dimension())
            throw std::out_of_range("VectorDouble: index out of range");
        return components_[index];
    }

    Vector &operator+=(const Vector &other) & {
        check_dimension(other.dimension(), "operator+=");
        for (size_t i = 0; i < dimension(); ++i)
            components_[i] += other.components_[i];
        return *this;
    }

    Vector &operator-=(const Vector &other) & {
        check_dimension(other.dimension(), "operator-=");
        for (size_t i = 0; i < dimension(); ++i)
            components_[i] -= other.components_[i];
        return *this;
    }

    Vector &operator*=(double scalar) & {
        for (double &c : components_)
            c *= scalar;
        return *this;
    }

    Vector &operator/=(double scalar) & {
        if (scalar == 0.0)
            throw std::domain_error("VectorDouble: division by zero");
        for (double &c : components_)
            c /= scalar;
        return *this;
    }

    Vector operator+() const { return *this; }

    Vector operator-() const {
        Vector copy = *this;
        return copy *= -1.0;
    }

    friend Vector operator+(Vector a, const Vector &b) { return a += b; }
    friend Vector operator-(Vector a, const Vector &b) { return a -= b; }
    friend Vector operator*(Vector v, double s)              { return v *= s; }
    friend Vector operator*(double s, Vector v)              { return v *= s; }
    friend Vector operator/(Vector v, double s)              { return v /= s; }

    friend bool operator==(const Vector &a, const Vector &b) {
        if (a.dimension() != b.dimension()) return false;
        for (size_t i = 0; i < a.dimension(); ++i)
            if (a[i] != b[i]) return false;
        return true;
    }

    friend bool operator!=(const Vector &a, const Vector &b) {
        return !(a == b);
    }

    double dot(const Vector &other) const {
        check_dimension(other.dimension(), "dot");
        double result = 0.0;
        for (size_t i = 0; i < dimension(); ++i)
            result += components_[i] * other.components_[i];
        return result;
    }

    double norm() const { return std::sqrt(dot(*this)); }

    Vector normalize() const {
        check_non_zero();
        double n = norm();
        if (n == 0.0)
            throw std::domain_error("VectorDouble: cannot normalize zero vector");
        return *this / n;
    }

    Vector cross_3d(const Vector &other) const {
        check_dimension(3, "cross_3d");
        other.check_dimension(3, "cross_3d");
        return Vector{
            components_[1] * other.components_[2] - components_[2] * other.components_[1],
            components_[2] * other.components_[0] - components_[0] * other.components_[2],
            components_[0] * other.components_[1] - components_[1] * other.components_[0]
        };
    }

    Vector cross_7d(const Vector &other) const {
        check_dimension(7, "cross_7d");
        other.check_dimension(7, "cross_7d");
        const auto &a = components_;
        const auto &b = other.components_;
        return Vector{
            a[1]*b[3] - a[3]*b[1] + a[2]*b[6] - a[6]*b[2] + a[4]*b[5] - a[5]*b[4],
            a[2]*b[4] - a[4]*b[2] + a[3]*b[0] - a[0]*b[3] + a[5]*b[6] - a[6]*b[5],
            a[3]*b[5] - a[5]*b[3] + a[4]*b[1] - a[1]*b[4] + a[6]*b[0] - a[0]*b[6],
            a[4]*b[6] - a[6]*b[4] + a[5]*b[2] - a[2]*b[5] + a[0]*b[1] - a[1]*b[0],
            a[5]*b[0] - a[0]*b[5] + a[6]*b[3] - a[3]*b[6] + a[1]*b[2] - a[2]*b[1],
            a[6]*b[1] - a[1]*b[6] + a[0]*b[4] - a[4]*b[0] + a[2]*b[3] - a[3]*b[2],
            a[0]*b[2] - a[2]*b[0] + a[1]*b[5] - a[5]*b[1] + a[3]*b[4] - a[4]*b[3]
        };
    }

    static double triple_product_3d(const Vector &a,
                                    const Vector &b,
                                    const Vector &c) {
        return a.dot(b.cross_3d(c));
    }

    static double triple_product_7d(const Vector &a,
                                    const Vector &b,
                                    const Vector &c) {
        return a.dot(b.cross_7d(c));
    }

    static Vector zero(size_t dimension) {
        return Vector(dimension);
    }

    static Vector basis_vector(size_t dimension, size_t index) {
        if (index >= dimension)
            throw std::out_of_range("VectorDouble: basis vector index out of range");
        Vector result(dimension);
        result[index] = 1.0;
        return result;
    }

    bool is_zero() const {
        for (double c : components_)
            if (c != 0.0) return false;
        return true;
    }

    bool is_orthogonal_to(const Vector &other) const {
        return dot(other) == 0.0;
    }

    friend double angle_between(const Vector &a, const Vector &b) {
        a.check_dimension(b.dimension(), "angle_between");
        if (a.is_zero() || b.is_zero())
            throw std::domain_error("VectorDouble: angle with zero vector is undefined");
        double denom = a.norm() * b.norm();
        if (denom == 0.0)
            throw std::domain_error("VectorDouble: vectors too small for angle calculation");
        double cosine = a.dot(b) / denom;

        cosine = std::fmax(-1.0, std::fmin(1.0, cosine));
        return std::acos(cosine);
    }

    friend bool are_orthogonal(const Vector &a, const Vector &b) {
        return a.is_orthogonal_to(b);
    }

    friend bool are_collinear(const Vector &a, const Vector &b) {
        if (a.dimension() != b.dimension()) return false;
        if (a.is_zero() || b.is_zero()) return true;

        for (size_t i = 0; i < a.dimension(); ++i)
            for (size_t j = i + 1; j < a.dimension(); ++j)
                if (a[i] * b[j] - a[j] * b[i] != 0.0) return false;
        return true;
    }

    friend bool is_point_on_segment(const Vector &pt,
                                    const Vector &a,
                                    const Vector &b,
                                    double eps) {
        pt.check_dimension(a.dimension(), "is_point_on_segment");
        pt.check_dimension(b.dimension(), "is_point_on_segment");
        Vector ab = b - a;
        Vector ap = pt - a;

        for (size_t i = 0; i < ab.dimension(); ++i)
            for (size_t j = i + 1; j < ab.dimension(); ++j)
                if (std::fabs(ap[i] * ab[j] - ap[j] * ab[i]) > eps) return false;

        double ab2 = ab.dot(ab);
        if (ab2 == 0.0) return (pt == a);
        double t = ap.dot(ab) / ab2;
        return (t >= -eps && t <= 1.0 + eps);
    }

    friend double point_to_segment_distance(const Vector &pt,
                                            const Vector &a,
                                            const Vector &b) {
        pt.check_dimension(a.dimension(), "point_to_segment_distance");
        pt.check_dimension(b.dimension(), "point_to_segment_distance");
        Vector ab = b - a;
        Vector ap = pt - a;
        double ab2 = ab.dot(ab);
        if (ab2 == 0.0) return (pt - a).norm();
        double t = ap.dot(ab) / ab2;
        t = std::fmax(0.0, std::fmin(1.0, t));
        Vector closest = a + ab * t;
        return (pt - closest).norm();
    }

    std::string to_string() const {
        if (components_.empty()) return "[]";
        std::string result = "[";
        for (size_t i = 0; i < components_.size(); ++i) {
            if (i != 0) result += ", ";
            result += std::to_string(components_[i]);
        }
        result += "]";
        return result;
    }
};

inline double angle_between(const Vector &a, const Vector &b);
inline bool   are_orthogonal(const Vector &a, const Vector &b);
inline bool   are_collinear(const Vector &a, const Vector &b);
inline bool   is_point_on_segment(const Vector &pt, const Vector &a, const Vector &b, double eps = 1e-10);
inline double point_to_segment_distance(const Vector &pt, const Vector &a, const Vector &b);
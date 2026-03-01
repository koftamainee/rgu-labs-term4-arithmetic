#include <iostream>
#include <complex>
#include "bigfloat.h"
#include "poly_tostring.hpp"

using CBF = std::complex<bigfloat>;

static std::vector<CBF> compute_omega_powers(const size_t n) {
    const bigfloat two_pi_over_n = bigfloat(2) * bigfloat::PI() / bigfloat(n);
    std::vector<CBF> result(n);
    for (size_t i = 0; i < n; ++i)
        result[i] = CBF(cos(two_pi_over_n * bigfloat(i)), sin(two_pi_over_n * bigfloat(i)));
    return result;
}

static std::vector<CBF> poly_rem_xn_minus_1(const std::vector<CBF>& f, const size_t n) {
    std::vector<CBF> result(n, CBF(0));
    for (size_t i = 0; i < f.size(); ++i) result[i % n] += f[i];
    return result;
}

static std::vector<CBF> poly_rem_xn_plus_1(const std::vector<CBF>& f, const size_t n) {
    std::vector<CBF> result(n, CBF(0));
    for (size_t i = 0; i < f.size(); ++i) {
        if ((i / n) % 2 == 0) result[i % n] += f[i];
        else                   result[i % n] -= f[i];
    }
    return result;
}

static std::vector<CBF> poly_mul_classical(const std::vector<CBF>& f, const std::vector<CBF>& g) {
    std::vector<CBF> result(f.size() + g.size() - 1, CBF(0));
    for (size_t i = 0; i < f.size(); ++i)
        for (size_t j = 0; j < g.size(); ++j)
            result[i + j] += f[i] * g[j];
    return result;
}

static std::vector<CBF> fast_convolution(
    const std::vector<CBF>& f,
    const std::vector<CBF>& g,
    const std::vector<CBF>& omega_powers,
    const size_t n,
    const size_t k,
    const size_t d)
{
    if (k <= d)
        return poly_rem_xn_minus_1(poly_mul_classical(f, g), n);

    const size_t half = n / 2;

    const auto f0 = poly_rem_xn_minus_1(f, half);
    const auto f1 = poly_rem_xn_plus_1(f, half);
    const auto g0 = poly_rem_xn_minus_1(g, half);
    const auto g1 = poly_rem_xn_plus_1(g, half);

    std::vector<CBF> half_omega(half);
    for (size_t i = 0; i < half; ++i) half_omega[i] = omega_powers[i * 2];

    std::vector<CBF> f1_wx(half), g1_wx(half);
    for (size_t i = 0; i < half; ++i) {
        f1_wx[i] = f1[i] * omega_powers[i];
        g1_wx[i] = g1[i] * omega_powers[i];
    }

    const auto h0 = fast_convolution(f0, g0, half_omega, half, k - 1, d);
    const auto h1 = fast_convolution(f1_wx, g1_wx, half_omega, half, k - 1, d);

    std::vector<CBF> result(n, CBF(0));
    for (size_t i = 0; i < half; ++i) {
        result[i]        = (h0[i] + h1[i]) * CBF(bigfloat(1, 2));
        result[i + half] = (h0[i] - h1[i]) * CBF(bigfloat(1, 2));
    }
    return result;
}

int main() {
    std::cout << "=== Fast Convolution ===\n\n";

    constexpr size_t k = 2;
    constexpr size_t n = 1 << k;

    const std::vector<CBF> f = {CBF(1), CBF(2), CBF(3), CBF(4)};
    const std::vector<CBF> g = {CBF(1), CBF(1), CBF(1), CBF(1)};

    const auto omega_powers = compute_omega_powers(n);

    std::cout << "Omega powers (primitive " << n << "th roots of unity):\n";
    for (size_t i = 0; i < n; ++i)
        std::cout << "  omega^" << i << " = "
                  << complex_poly_tostring({omega_powers[i]}, CBF(0)) << "\n";

    const auto classical = poly_rem_xn_minus_1(poly_mul_classical(f, g), n);
    std::cout << "\nClassical (mod x^n - 1) = "
              << complex_poly_tostring(classical, CBF(0)) << "\n";

    for (size_t d = 0; d <= k; ++d) {
        std::cout << "\n--- d = " << d << " ---\n";
        const auto result = fast_convolution(f, g, omega_powers, n, k, d);
        std::cout << "f * g = " << complex_poly_tostring(result, CBF(0)) << "\n";
    }

    return 0;
}
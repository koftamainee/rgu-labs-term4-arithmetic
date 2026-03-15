#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "power_series.hpp"

void print_series(const std::string &label, const PowerSeries &ps, size_t n) {
    std::cout << label << "\n";
    for (size_t i = 0; i < n && i < ps.size(); ++i) {
        std::cout << "  a[" << i << "] = "
                  << std::fixed << std::setprecision(10) << ps[i] << "\n";
    }
    std::cout << "\n";
}

void part_ii(const PowerSeries &f) {
    std::cout << "--- Part ii): first 3 coefficients of f^{-1} ---\n";
    const PowerSeries g = f.inverse(3);
    print_series("  g = f^{-1} (first 3 coefficients)", g, 3);
}

void part_iii_exp(size_t M, size_t n) {
    std::vector<double> coeffs(M + 1);
    double factorial = 1.0;
    for (size_t k = 0; k <= M; ++k) {
        if (k > 0) factorial *= static_cast<double>(k);
        coeffs[k] = 1.0 / factorial;
    }
    const PowerSeries f(coeffs);

    std::cout << "--- Part iii), Example (1): f = sum_{k=0}^{" << M << "} x^k/k! ---\n";
    print_series("  f (first " + std::to_string(M + 1) + " coefficients)", f, M + 1);

    const PowerSeries g = f.inverse(n);
    print_series("  f^{-1} (first " + std::to_string(n) + " coefficients)", g, n);

    const PowerSeries check = (f * g).first_n(n);
    std::cout << "  Verification f * f^{-1} (first " << n << " coefficients):\n";
    for (size_t i = 0; i < n; ++i)
        std::cout << "    [x^" << i << "] = "
                  << std::fixed << std::setprecision(10) << check[i] << "\n";
    std::cout << "\n";
}

void part_iii_quadratic(size_t n) {
    const PowerSeries f(std::vector<double>{-1.0, -1.0, 1.0});

    std::cout << "--- Part iii), Example (2): f = x^2 - x - 1 ---\n";
    print_series("  f", f, f.size());

    const PowerSeries g = f.inverse(n);
    print_series("  f^{-1} (first " + std::to_string(n) + " coefficients)", g, n);

    const PowerSeries check = (f * g).first_n(n);
    std::cout << "  Verification f * f^{-1} (first " << n << " coefficients):\n";
    for (size_t i = 0; i < n; ++i)
        std::cout << "    [x^" << i << "] = "
                  << std::fixed << std::setprecision(10) << check[i] << "\n";
    std::cout << "\n";
}

int main() {
    std::cout << "=== Task 13 – Formal power series inversion ===\n\n";

    {
        const PowerSeries f_example(std::vector<double>{2.0, 3.0, -1.0});
        std::cout << "--- Part ii): example f = 2 + 3x - x^2 ---\n";
        print_series("  f", f_example, 3);
        part_ii(f_example);
    }

    const size_t n = 10;
    const size_t M = 5;

    part_iii_exp(M, n);
    part_iii_quadratic(n);

    return 0;
}
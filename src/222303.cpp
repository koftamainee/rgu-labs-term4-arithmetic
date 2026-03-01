#include <iostream>
#include <vector>
#include <complex>
#include "fft.hpp"

using Array = std::vector<CD>;

Array dft_nd(Array data, const std::vector<size_t>& dims) {
    const size_t n = data.size();
    size_t stride = 1;
    for (int dim = static_cast<int>(dims.size()) - 1; dim >= 0; --dim) {
        const size_t nd = dims[dim];
        const size_t num_transforms = n / nd;
        for (size_t t = 0; t < num_transforms; ++t) {
            const size_t block  = t / stride;
            const size_t offset = t % stride;
            const size_t start  = block * nd * stride + offset;
            std::vector<CD> slice(nd);
            for (size_t i = 0; i < nd; ++i)
                slice[i] = data[start + i * stride];
            fft(slice, false);
            for (size_t i = 0; i < nd; ++i)
                data[start + i * stride] = slice[i];
        }
        stride *= nd;
    }
    return data;
}

int main() {
    const std::vector<size_t> dims = {4, 4};
    const size_t rows = dims[0], cols = dims[1];
    const size_t n = rows * cols;

    Array data(n);
    for (size_t i = 0; i < n; ++i)
        data[i] = CD(static_cast<double>(i + 1), 0.0);

    std::cout << "Multidimensional DFT via sequential 1D FFTs\n";
    std::cout << "n = " << n << ",  dims = " << rows << " x " << cols << "\n\n";

    std::cout << "Input:\n";
    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < cols; ++j)
            std::cout << "  " << static_cast<int>(data[i * cols + j].real());
        std::cout << "\n";
    }

    const Array result = dft_nd(data, dims);

    std::cout << "\nOutput (2D DFT, e^{+2pi*i} convention):\n";
    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < cols; ++j) {
            const CD& v = result[i * cols + j];
            std::cout << "  (" << v.real() << ", " << v.imag() << "i)";
        }
        std::cout << "\n";
    }

    return 0;
}
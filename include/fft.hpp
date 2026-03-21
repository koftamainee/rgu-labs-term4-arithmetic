#pragma once

#include <complex>
#include <cmath>
#include <vector>
#include <stdexcept>

template <typename F = double>
void fft(std::vector<std::complex<F>>& arr, bool inverse, F pi = F(M_PI)) {
  const std::size_t n = arr.size();
  if (n == 1) {
    return;
  }
  if ((n & (n - 1)) != 0) {
    throw std::invalid_argument("fft: size must be a power of two");
  }

  std::vector<std::complex<F>> even(n / 2);
  std::vector<std::complex<F>> odd(n / 2);
  for (std::size_t i = 0; i < n / 2; ++i) {
    even[i] = arr[2 * i];
    odd[i] = arr[2 * i + 1];
  }

  fft(even, inverse, pi);
  fft(odd, inverse, pi);

  const F ang = F(2) * pi / F(n) * (inverse ? F(-1) : F(1));
  std::complex<F> w(1, 0);
  const std::complex<F> wn(std::cos(ang), std::sin(ang));

  for (std::size_t i = 0; i < n / 2; ++i) {
    const std::complex<F> t = w * odd[i];
    arr[i] = even[i] + t;
    arr[i + n / 2] = even[i] - t;
    w *= wn;
  }

  if (inverse) {
    for (auto& x : arr) {
      x /= F(2);
    }
  }
}

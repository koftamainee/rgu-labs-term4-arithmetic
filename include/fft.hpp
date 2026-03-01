#ifndef RGU_LABS_TERM4_ARITHMETIC_FFT_HPP
#define RGU_LABS_TERM4_ARITHMETIC_FFT_HPP
#include <complex>
#include <vector>

using CD = std::complex<double>;
const double PI = acos(-1.0);

inline void fft(std::vector<CD>& arr, bool inverse) {
  const size_t n = arr.size();
  if (n == 1) { return; }

  std::vector<CD> even(n / 2), odd(n / 2);
  for (size_t i = 0; i < n / 2; ++i) {
    even[i] = arr[2 * i];
    odd[i] = arr[2 * i + 1];
  }

  fft(even, inverse);
  fft(odd, inverse);

  const double ang = 2 * PI / n * (inverse ? -1 : 1);
  CD w(1);
  const CD wn(std::cos(ang), std::sin(ang));
  for (size_t i = 0; i < n / 2; ++i) {
    arr[i] = even[i] + w * odd[i];
    arr[i + n / 2] = even[i] - w * odd[i];
    w *= wn;
  }
}

#endif //RGU_LABS_TERM4_ARITHMETIC_FFT_HPP
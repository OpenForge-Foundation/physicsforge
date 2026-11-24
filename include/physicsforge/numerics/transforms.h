// include/physicsforge/numerics/transforms.h
//
// Fourier and Laplace transforms.
//
//   * fft / ifft: discrete Fourier transform of any length in O(n log n). Powers of
//     two use an iterative radix-2 Cooley-Tukey kernel with a precomputed twiddle
//     table; other lengths use Bluestein's chirp-z algorithm on top of it.
//     Convention (matches NumPy): X_k = sum_n x_n exp(-2 pi i k n / N); ifft
//     includes the 1/N factor.
//   * Fourier series: real trigonometric coefficients from one period of samples,
//         f(x) ~ a_0 / 2 + sum_k [a_k cos(2 pi k x / P) + b_k sin(2 pi k x / P)].
//   * Laplace transform: F(s) = int_0^inf f(t) exp(-s t) dt for real s > 0, by
//     adaptive quadrature after mapping [0, inf) onto [0, 1).
//   * Inverse Laplace transform: Gaver-Stehfest algorithm, which needs F only on the
//     positive real axis. Accurate to ~4-6 significant digits for smooth,
//     non-oscillatory f(t); unsuitable for strongly oscillatory or discontinuous f.

#ifndef PHYSICSFORGE_NUMERICS_TRANSFORMS_H
#define PHYSICSFORGE_NUMERICS_TRANSFORMS_H

#include <complex>
#include <span>
#include <vector>

#include "physicsforge/config.h"
#include "physicsforge/numerics/functions.h"

namespace physicsforge::numerics {

    // Forward discrete Fourier transform (any length; empty input yields empty output).
    std::vector<std::complex<double>> fft(std::span<const std::complex<double>> input);

    // Inverse discrete Fourier transform, normalised by 1/N.
    std::vector<std::complex<double>> ifft(std::span<const std::complex<double>> input);

    // Sample frequencies for an FFT of length n and sample spacing d (NumPy fftfreq order).
    std::vector<double> fft_frequencies(int n, double sample_spacing);

    // Real Fourier-series coefficients (a_0..a_K, b_0..b_K with b_0 = 0) of a function sampled
    // at N uniform points over one period, end point excluded. Requires 2K < N.
    struct FourierSeries {
        std::vector<double> a;
        std::vector<double> b;
        double period = 1.0;
    };
    FourierSeries fourier_series_coefficients(std::span<const double> samples, int n_terms, double period);

    // Evaluates a truncated Fourier series at the given points.
    std::vector<double> fourier_series_evaluate(const FourierSeries &series, std::span<const double> points);

    // Laplace transform F(s) of f for real s > 0.
    double laplace_transform(const ScalarFunction &f, double s,
                             double tolerance = config::DEFAULT_INTEGRATION_TOLERANCE);

    // Inverse Laplace transform f(t) for t > 0 via Gaver-Stehfest with an even number of terms
    // in [2, 30] (default 14).
    double inverse_laplace_transform(const ScalarFunction &transform, double t,
                                     int n_terms = config::DEFAULT_STEHFEST_TERMS);

} // namespace physicsforge::numerics

#endif // PHYSICSFORGE_NUMERICS_TRANSFORMS_H

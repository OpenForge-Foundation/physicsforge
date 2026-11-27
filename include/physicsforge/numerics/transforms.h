/**
 * @file transforms.h
 * @brief Fourier and Laplace transforms.
 *
 * - **FFT / IFFT** of any length in O(n log n). Powers of two use an iterative radix-2 Cooley–Tukey kernel
 *   with a precomputed twiddle table; other lengths use Bluestein's chirp-z algorithm on top of it.
 *   Convention (matches NumPy): \f$ X_k = \sum_n x_n e^{-2\pi i k n / N} \f$, and the inverse includes 1/N.
 * - **Fourier series** from one period of samples,
 *   \f$ f(x) \approx \tfrac{a_0}{2} + \sum_k [a_k \cos(2\pi k x / P) + b_k \sin(2\pi k x / P)] \f$.
 * - **Laplace transform** \f$ F(s) = \int_0^\infty f(t) e^{-st}\,dt \f$ for real \f$ s > 0 \f$, by adaptive
 *   quadrature after mapping \f$ [0, \infty) \f$ onto \f$ [0, 1) \f$.
 * - **Inverse Laplace transform** by the Gaver–Stehfest algorithm, which needs \f$ F \f$ only on the positive
 *   real axis. Accurate to about 3–6 significant digits for smooth, non-oscillatory \f$ f(t) \f$.
 *
 * @ingroup numerics
 */

#ifndef PHYSICSFORGE_NUMERICS_TRANSFORMS_H
#define PHYSICSFORGE_NUMERICS_TRANSFORMS_H

#include <complex>
#include <span>
#include <vector>

#include "physicsforge/config.h"
#include "physicsforge/numerics/calculus.h"

namespace physicsforge::numerics {

    /**
     * @brief Forward discrete Fourier transform.
     * @param input Samples \f$ x_0 \ldots x_{N-1} \f$ (any length; empty yields empty).
     * @return Spectrum \f$ X_0 \ldots X_{N-1} \f$.
     */
    std::vector<std::complex<double>> fft(std::span<const std::complex<double>> input);

    /**
     * @brief Inverse discrete Fourier transform, normalised by \f$ 1/N \f$.
     * @param input Spectrum \f$ X_0 \ldots X_{N-1} \f$.
     * @return Samples \f$ x_0 \ldots x_{N-1} \f$.
     */
    std::vector<std::complex<double>> ifft(std::span<const std::complex<double>> input);

    /**
     * @brief Sample frequencies of an FFT in NumPy @c fftfreq order.
     * @param n Transform length; must be non-negative.
     * @param sample_spacing Time between samples \f$ d \f$; must be positive.
     * @return Frequencies in cycles per unit of @p sample_spacing.
     * @throws std::invalid_argument If @p n is negative or @p sample_spacing is not positive.
     */
    std::vector<double> fft_frequencies(int n, double sample_spacing);

    /** @brief Real trigonometric Fourier-series coefficients. */
    struct FourierSeries {
        std::vector<double> a; ///< Cosine coefficients \f$ a_0 \ldots a_K \f$.
        std::vector<double> b; ///< Sine coefficients \f$ b_0 \ldots b_K \f$, with \f$ b_0 = 0 \f$.
        double period = 1.0; ///< Period \f$ P \f$.
    };

    /**
     * @brief Fourier-series coefficients of a function sampled over one period (end point excluded).
     * @param samples \f$ N \f$ uniform samples over one period.
     * @param n_terms Highest harmonic \f$ K \f$; must satisfy \f$ 2K < N \f$ (Nyquist).
     * @param period Period \f$ P \f$; must be positive.
     * @return Coefficients \f$ a_0 \ldots a_K \f$ and \f$ b_0 \ldots b_K \f$.
     * @throws std::invalid_argument If @p n_terms exceeds the Nyquist limit or @p period is not positive.
     */
    FourierSeries fourier_series_coefficients(std::span<const double> samples, int n_terms, double period);

    /**
     * @brief Evaluates a truncated Fourier series.
     * @param series Coefficients from fourier_series_coefficients().
     * @param points Evaluation points (any real values; the series is periodic).
     * @return Series values at @p points.
     */
    std::vector<double> fourier_series_evaluate(const FourierSeries &series, std::span<const double> points);

    /**
     * @brief Laplace transform \f$ F(s) = \int_0^\infty f(t) e^{-st}\,dt \f$.
     * @param f Function of time.
     * @param s Real transform variable; must be positive.
     * @param tolerance Absolute quadrature tolerance; must be positive.
     * @return \f$ F(s) \f$.
     * @throws std::invalid_argument If @p s or @p tolerance is not positive.
     */
    double laplace_transform(const ScalarFunction &f, double s,
                             double tolerance = config::DEFAULT_INTEGRATION_TOLERANCE);

    /**
     * @brief Inverse Laplace transform \f$ f(t) \f$ by the Gaver–Stehfest algorithm.
     * @param transform Laplace-domain function \f$ F(s) \f$, evaluated only for real \f$ s > 0 \f$.
     * @param t Time; must be positive.
     * @param n_terms Even number of terms in [2, 30]. More terms are not always better in double precision.
     * @return Approximation of \f$ f(t) \f$.
     * @throws std::invalid_argument If @p t is not positive or @p n_terms is odd or out of range.
     */
    double inverse_laplace_transform(const ScalarFunction &transform, double t,
                                     int n_terms = config::DEFAULT_STEHFEST_TERMS);

} // namespace physicsforge::numerics

#endif // PHYSICSFORGE_NUMERICS_TRANSFORMS_H

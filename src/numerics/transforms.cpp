// src/numerics/transforms.cpp
#include "physicsforge/numerics/transforms.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <numbers>
#include <utility>

#include "physicsforge/numerics/integrate.h"
#include "physicsforge/validation.h"

namespace physicsforge::numerics {

    using cd = std::complex<double>;

    namespace {
        constexpr double TWO_PI = 2.0 * std::numbers::pi;

        bool is_power_of_two(const std::size_t n) { return n != 0 && (n & (n - 1)) == 0; }

        std::size_t next_power_of_two(std::size_t n) {
            std::size_t p = 1;
            while (p < n) {
                p <<= 1;
            }
            return p;
        }

        // In-place iterative radix-2 Cooley-Tukey FFT. n must be a power of two.
        // Twiddles come from a table computed once per call with direct sin/cos, which is
        // markedly more accurate than accumulating w *= w_len across a stage.
        void fft_radix2_inplace(std::vector<cd> &a, const bool inverse) {
            const std::size_t n = a.size();
            if (n <= 1) {
                return;
            }

            // Bit-reversal permutation.
            for (std::size_t i = 1, j = 0; i < n; ++i) {
                std::size_t bit = n >> 1;
                for (; j & bit; bit >>= 1) {
                    j ^= bit;
                }
                j ^= bit;
                if (i < j) {
                    std::swap(a[i], a[j]);
                }
            }

            const double sign = inverse ? 1.0 : -1.0;
            std::vector<cd> twiddle(n / 2);
            for (std::size_t k = 0; k < n / 2; ++k) {
                const double angle = sign * TWO_PI * static_cast<double>(k) / static_cast<double>(n);
                twiddle[k] = {std::cos(angle), std::sin(angle)};
            }

            for (std::size_t len = 2; len <= n; len <<= 1) {
                const std::size_t half = len >> 1;
                const std::size_t stride = n / len;
                for (std::size_t start = 0; start < n; start += len) {
                    for (std::size_t k = 0; k < half; ++k) {
                        const cd u = a[start + k];
                        const cd v = a[start + k + half] * twiddle[k * stride];
                        a[start + k] = u + v;
                        a[start + k + half] = u - v;
                    }
                }
            }
        }

        // Unnormalised DFT of arbitrary length via Bluestein's chirp-z algorithm:
        // X_k = w_k * sum_n (x_n w_n) conj(w_{k-n}), with w_k = exp(sign * i pi k^2 / n),
        // evaluated as a power-of-two circular convolution.
        std::vector<cd> bluestein(std::span<const cd> input, const bool inverse) {
            const std::size_t n = input.size();
            const std::size_t m = next_power_of_two(2 * n - 1);
            const double sign = inverse ? 1.0 : -1.0;

            // k^2 mod 2n keeps the chirp argument small, preserving accuracy for large n.
            std::vector<cd> chirp(n);
            const auto two_n = static_cast<std::uint64_t>(2 * n);
            for (std::size_t k = 0; k < n; ++k) {
                const auto k64 = static_cast<std::uint64_t>(k);
                const double angle =
                        sign * std::numbers::pi * static_cast<double>((k64 * k64) % two_n) / static_cast<double>(n);
                chirp[k] = {std::cos(angle), std::sin(angle)};
            }

            std::vector<cd> a(m, cd{0.0, 0.0});
            std::vector<cd> b(m, cd{0.0, 0.0});
            for (std::size_t k = 0; k < n; ++k) {
                a[k] = input[k] * chirp[k];
            }
            b[0] = std::conj(chirp[0]);
            for (std::size_t k = 1; k < n; ++k) {
                b[k] = std::conj(chirp[k]);
                b[m - k] = std::conj(chirp[k]);
            }

            fft_radix2_inplace(a, false);
            fft_radix2_inplace(b, false);
            for (std::size_t i = 0; i < m; ++i) {
                a[i] *= b[i];
            }
            fft_radix2_inplace(a, true);

            const double inv_m = 1.0 / static_cast<double>(m);
            std::vector<cd> out(n);
            for (std::size_t k = 0; k < n; ++k) {
                out[k] = a[k] * inv_m * chirp[k];
            }
            return out;
        }

        std::vector<cd> dft(std::span<const cd> input, const bool inverse) {
            if (is_power_of_two(input.size())) {
                std::vector<cd> data(input.begin(), input.end());
                fft_radix2_inplace(data, inverse);
                return data;
            }
            if (input.empty()) {
                return {};
            }
            return bluestein(input, inverse);
        }

        // Gaver-Stehfest weights V_k, k = 1..N.
        std::vector<double> stehfest_weights(const int n_terms) {
            const int half = n_terms / 2;
            std::vector<double> factorial(static_cast<std::size_t>(n_terms) + 1, 1.0);
            for (int i = 1; i <= n_terms; ++i) {
                factorial[static_cast<std::size_t>(i)] = factorial[static_cast<std::size_t>(i - 1)] * i;
            }
            auto fact = [&](const int i) { return factorial[static_cast<std::size_t>(i)]; };

            std::vector<double> weights(static_cast<std::size_t>(n_terms));
            for (int k = 1; k <= n_terms; ++k) {
                double sum = 0.0;
                for (int j = (k + 1) / 2; j <= std::min(k, half); ++j) {
                    sum += std::pow(static_cast<double>(j), half) * fact(2 * j) /
                           (fact(half - j) * fact(j) * fact(j - 1) * fact(k - j) * fact(2 * j - k));
                }
                const double sign = ((k + half) % 2 == 0) ? 1.0 : -1.0;
                weights[static_cast<std::size_t>(k - 1)] = sign * sum;
            }
            return weights;
        }
    } // namespace

    std::vector<cd> fft(const std::span<const cd> input) { return dft(input, false); }

    std::vector<cd> ifft(const std::span<const cd> input) {
        std::vector<cd> out = dft(input, true);
        if (!out.empty()) {
            const double inv_n = 1.0 / static_cast<double>(out.size());
            for (cd &v: out) {
                v *= inv_n;
            }
        }
        return out;
    }

    std::vector<double> fft_frequencies(const int n, const double sample_spacing) {
        if (n < 1) {
            detail::throw_invalid("n", "at least 1", static_cast<double>(n));
        }
        detail::require_positive(sample_spacing, "sample_spacing");
        std::vector<double> freqs(static_cast<std::size_t>(n));
        const double scale = 1.0 / (static_cast<double>(n) * sample_spacing);
        const int positive_count = (n - 1) / 2 + 1;
        for (int i = 0; i < n; ++i) {
            const int k = i < positive_count ? i : i - n;
            freqs[static_cast<std::size_t>(i)] = static_cast<double>(k) * scale;
        }
        return freqs;
    }

    FourierSeries fourier_series_coefficients(const std::span<const double> samples, const int n_terms,
                                              const double period) {
        detail::require_positive(period, "period");
        if (n_terms < 0) {
            detail::throw_invalid("n_terms", "non-negative", static_cast<double>(n_terms));
        }
        const std::size_t n = samples.size();
        if (n == 0 || 2 * static_cast<std::size_t>(n_terms) >= n) {
            detail::throw_invalid("need more than 2 * n_terms samples per period (Nyquist limit)");
        }

        std::vector<cd> data(n);
        for (std::size_t i = 0; i < n; ++i) {
            detail::require_finite(samples[i], "samples");
            data[i] = {samples[i], 0.0};
        }
        const std::vector<cd> spectrum = fft(data);

        FourierSeries series;
        series.period = period;
        series.a.resize(static_cast<std::size_t>(n_terms) + 1);
        series.b.resize(static_cast<std::size_t>(n_terms) + 1);
        const double scale = 2.0 / static_cast<double>(n);
        for (int k = 0; k <= n_terms; ++k) {
            const auto idx = static_cast<std::size_t>(k);
            series.a[idx] = scale * spectrum[idx].real();
            series.b[idx] = -scale * spectrum[idx].imag();
        }
        series.b[0] = 0.0;
        return series;
    }

    std::vector<double> fourier_series_evaluate(const FourierSeries &series, const std::span<const double> points) {
        detail::require_positive(series.period, "period");
        if (series.a.empty() || series.a.size() != series.b.size()) {
            detail::throw_invalid("Fourier series a and b must be non-empty and the same length");
        }
        const double omega = TWO_PI / series.period;
        std::vector<double> out(points.size());
        for (std::size_t p = 0; p < points.size(); ++p) {
            detail::require_finite(points[p], "points");
            const double theta = omega * points[p];
            // Evaluate cos(k theta), sin(k theta) by the angle-addition recurrence: one sin/cos per point.
            const double c1 = std::cos(theta);
            const double s1 = std::sin(theta);
            double ck = 1.0;
            double sk = 0.0;
            double value = 0.5 * series.a[0];
            for (std::size_t k = 1; k < series.a.size(); ++k) {
                const double c_next = ck * c1 - sk * s1;
                sk = sk * c1 + ck * s1;
                ck = c_next;
                value += series.a[k] * ck + series.b[k] * sk;
            }
            out[p] = value;
        }
        return out;
    }

    double laplace_transform(const ScalarFunction &f, const double s, const double tolerance) {
        detail::require_positive(s, "s");
        detail::require_positive(tolerance, "tolerance");
        // Scale time by s (tau = s t) so the kernel is exp(-tau) regardless of s, then map
        // tau in [0, inf) to u in [0, 1) with tau = u / (1 - u):
        //   F(s) = (1/s) int_0^1 f(tau/s) exp(-tau) / (1 - u)^2 du.
        const ScalarFunction integrand = [&f, s](const double u) {
            if (u >= 1.0) {
                return 0.0;
            }
            const double one_minus_u = 1.0 - u;
            const double tau = u / one_minus_u;
            const double weight = std::exp(-tau) / (one_minus_u * one_minus_u);
            if (weight == 0.0) {
                return 0.0; // Kernel underflow: skip evaluating f at enormous t.
            }
            return f(tau / s) * weight;
        };
        return integrate(integrand, 0.0, 1.0, tolerance * s) / s;
    }

    double inverse_laplace_transform(const ScalarFunction &transform, const double t, const int n_terms) {
        detail::require_positive(t, "t");
        if (n_terms < 2 || n_terms > 30 || n_terms % 2 != 0) {
            detail::throw_invalid("n_terms", "an even integer in [2, 30]", static_cast<double>(n_terms));
        }
        const std::vector<double> weights = stehfest_weights(n_terms);
        const double a = std::numbers::ln2 / t;
        double sum = 0.0;
        for (int k = 1; k <= n_terms; ++k) {
            sum += weights[static_cast<std::size_t>(k - 1)] * transform(k * a);
        }
        return a * sum;
    }

} // namespace physicsforge::numerics

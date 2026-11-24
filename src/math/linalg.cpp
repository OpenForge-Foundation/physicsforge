// src/math/linalg.cpp
#include "physicsforge/math/linalg.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>

#include "physicsforge/config.h"

namespace physicsforge::math {

    double dot(const std::span<const double> a, const std::span<const double> b) {
        if (a.size() != b.size()) {
            throw std::invalid_argument("dot: vectors must have the same length");
        }
        double sum = 0.0;
        for (std::size_t i = 0; i < a.size(); ++i) {
            sum += a[i] * b[i];
        }
        return sum;
    }

    double norm(const std::span<const double> a) {
        // Two-pass scaled norm (as in BLAS dnrm2) so that very large or very small
        // components neither overflow nor underflow when squared.
        double scale = 0.0;
        for (const double v: a) {
            scale = std::max(scale, std::abs(v));
        }
        if (scale == 0.0 || !std::isfinite(scale)) {
            return scale;
        }
        double sum = 0.0;
        for (const double v: a) {
            const double r = v / scale;
            sum += r * r;
        }
        return scale * std::sqrt(sum);
    }

    std::vector<double> solve_linear_system(const std::span<const double> matrix, const std::span<const double> rhs) {
        const std::size_t n = rhs.size();
        if (n == 0) {
            throw std::invalid_argument("solve_linear_system: system must be non-empty");
        }
        if (matrix.size() != n * n) {
            throw std::invalid_argument("solve_linear_system: matrix must be n x n where n = len(rhs)");
        }

        std::vector<double> a(matrix.begin(), matrix.end());
        std::vector<double> x(rhs.begin(), rhs.end());

        // Forward elimination with partial pivoting.
        for (std::size_t col = 0; col < n; ++col) {
            std::size_t pivot = col;
            double best = std::abs(a[col * n + col]);
            for (std::size_t row = col + 1; row < n; ++row) {
                const double candidate = std::abs(a[row * n + col]);
                if (candidate > best) {
                    best = candidate;
                    pivot = row;
                }
            }
            if (best < config::SINGULAR_PIVOT_THRESHOLD || !std::isfinite(best)) {
                throw std::domain_error("solve_linear_system: matrix is singular or ill-conditioned");
            }
            if (pivot != col) {
                std::swap_ranges(a.begin() + static_cast<std::ptrdiff_t>(col * n),
                                 a.begin() + static_cast<std::ptrdiff_t>((col + 1) * n),
                                 a.begin() + static_cast<std::ptrdiff_t>(pivot * n));
                std::swap(x[col], x[pivot]);
            }

            const double inv_pivot = 1.0 / a[col * n + col];
            for (std::size_t row = col + 1; row < n; ++row) {
                const double factor = a[row * n + col] * inv_pivot;
                if (factor == 0.0) {
                    continue;
                }
                for (std::size_t k = col; k < n; ++k) {
                    a[row * n + k] -= factor * a[col * n + k];
                }
                x[row] -= factor * x[col];
            }
        }

        // Back substitution.
        for (std::size_t i = n; i-- > 0;) {
            double sum = x[i];
            for (std::size_t k = i + 1; k < n; ++k) {
                sum -= a[i * n + k] * x[k];
            }
            x[i] = sum / a[i * n + i];
        }
        return x;
    }

} // namespace physicsforge::math

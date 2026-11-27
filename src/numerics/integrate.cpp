// src/numerics/integrate.cpp
#include "physicsforge/numerics/integrate.h"

#include <cmath>
#include <cstddef>
#include <string>

#include "physicsforge/validation.h"

namespace physicsforge::numerics {

    namespace {
        void require_samples(const std::span<const double> samples, const double spacing, const std::size_t minimum) {
            if (samples.size() < minimum) {
                detail::throw_invalid("samples must contain at least " + std::to_string(minimum) + " values");
            }
            detail::require_positive(spacing, "spacing");
        }

        // Recursive adaptive Simpson step on [a, b] with midpoint m and cached function values.
        double adaptive_simpson(const ScalarFunction &f, const double a, const double b, const double fa,
                                const double fm, const double fb, const double whole, const double tolerance,
                                const int depth) {
            const double m = 0.5 * (a + b);
            const double lm = 0.5 * (a + m);
            const double rm = 0.5 * (m + b);
            const double flm = f(lm);
            const double frm = f(rm);
            const double left = (m - a) / 6.0 * (fa + 4.0 * flm + fm);
            const double right = (b - m) / 6.0 * (fm + 4.0 * frm + fb);
            const double delta = left + right - whole;
            if (depth <= 0 || std::abs(delta) <= 15.0 * tolerance) {
                return left + right + delta / 15.0; // Richardson extrapolation.
            }
            return adaptive_simpson(f, a, m, fa, flm, fm, left, 0.5 * tolerance, depth - 1) +
                   adaptive_simpson(f, m, b, fm, frm, fb, right, 0.5 * tolerance, depth - 1);
        }
    } // namespace

    double trapezoid(const std::span<const double> samples, const double spacing) {
        require_samples(samples, spacing, 2);
        double interior = 0.0;
        for (std::size_t i = 1; i + 1 < samples.size(); ++i) {
            interior += samples[i];
        }
        return spacing * (0.5 * (samples.front() + samples.back()) + interior);
    }

    double simpson(const std::span<const double> samples, const double spacing) {
        require_samples(samples, spacing, 2);
        const std::size_t intervals = samples.size() - 1;
        if (intervals == 1) {
            return trapezoid(samples, spacing);
        }

        // Use 1/3 rule over an even number of leading intervals; if the total is odd, finish
        // with a 3/8 panel over the last three intervals.
        const std::size_t simpson_intervals = (intervals % 2 == 0) ? intervals : intervals - 3;
        double total = 0.0;
        if (simpson_intervals > 0) {
            double odd = 0.0;
            double even = 0.0;
            for (std::size_t i = 1; i < simpson_intervals; ++i) {
                (i % 2 == 1 ? odd : even) += samples[i];
            }
            total += spacing / 3.0 * (samples[0] + 4.0 * odd + 2.0 * even + samples[simpson_intervals]);
        }
        if (simpson_intervals != intervals) {
            const std::size_t s = simpson_intervals;
            total += 3.0 * spacing / 8.0 * (samples[s] + 3.0 * samples[s + 1] + 3.0 * samples[s + 2] + samples[s + 3]);
        }
        return total;
    }

    double integrate(const ScalarFunction &f, const double a, const double b, const double tolerance,
                     const int max_depth) {
        detail::require_finite(a, "a");
        detail::require_finite(b, "b");
        detail::require_positive(tolerance, "tolerance");
        if (max_depth < 1) {
            detail::throw_invalid("max_depth", "at least 1", static_cast<double>(max_depth));
        }
        if (a == b) {
            return 0.0;
        }
        const double fa = f(a);
        const double fb = f(b);
        const double m = 0.5 * (a + b);
        const double fm = f(m);
        const double whole = (b - a) / 6.0 * (fa + 4.0 * fm + fb);
        return adaptive_simpson(f, a, b, fa, fm, fb, whole, tolerance, max_depth);
    }

} // namespace physicsforge::numerics

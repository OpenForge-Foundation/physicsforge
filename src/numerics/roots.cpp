// src/numerics/roots.cpp
#include "physicsforge/numerics/roots.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "physicsforge/math/linalg.h"
#include "physicsforge/numerics/calculus.h"
#include "physicsforge/validation.h"

namespace physicsforge::numerics {

    namespace {
        void require_solver_settings(const double tolerance, const int max_iterations) {
            detail::require_positive(tolerance, "tolerance");
            if (max_iterations < 1) {
                detail::throw_invalid("max_iterations", "at least 1", static_cast<double>(max_iterations));
            }
        }

        bool step_converged(const double step, const double x, const double tolerance) {
            return std::abs(step) <= tolerance * std::max(1.0, std::abs(x));
        }
    } // namespace

    RootResult newton_raphson(const ScalarFunction &f, const ScalarFunction &derivative, const double x0,
                              const double tolerance, const int max_iterations) {
        detail::require_finite(x0, "x0");
        require_solver_settings(tolerance, max_iterations);

        RootResult result;
        double x = x0;
        double fx = f(x);
        for (int k = 1; k <= max_iterations; ++k) {
            result.iterations = k;
            if (fx == 0.0) {
                result.converged = true;
                break;
            }
            const double dfx = derivative(x);
            if (dfx == 0.0 || !std::isfinite(dfx) || !std::isfinite(fx)) {
                break; // Stationary point or blow-up: Newton step undefined.
            }
            const double step = fx / dfx;
            x -= step;
            fx = f(x);
            if (step_converged(step, x, tolerance)) {
                result.converged = std::isfinite(fx);
                break;
            }
        }
        result.root = x;
        result.residual = fx;
        return result;
    }

    RootResult newton_raphson(const ScalarFunction &f, const double x0, const double tolerance,
                              const int max_iterations) {
        const ScalarFunction numeric_derivative = [&f](const double x) { return derivative(f, x); };
        return newton_raphson(f, numeric_derivative, x0, tolerance, max_iterations);
    }

    RootResult bisection(const ScalarFunction &f, double a, double b, const double tolerance,
                         const int max_iterations) {
        detail::require_finite(a, "a");
        detail::require_finite(b, "b");
        require_solver_settings(tolerance, max_iterations);
        if (a > b) {
            std::swap(a, b);
        }

        double fa = f(a);
        const double fb = f(b);
        RootResult result;
        if (fa == 0.0 || fb == 0.0) {
            result.root = fa == 0.0 ? a : b;
            result.residual = 0.0;
            result.converged = true;
            return result;
        }
        if ((fa > 0.0) == (fb > 0.0)) {
            detail::throw_invalid("f(a) and f(b) must have opposite signs to bracket a root");
        }

        double mid = a;
        double fmid = fa;
        for (int k = 1; k <= max_iterations; ++k) {
            result.iterations = k;
            mid = a + 0.5 * (b - a); // Avoids overflow of (a + b).
            fmid = f(mid);
            if (fmid == 0.0 || 0.5 * (b - a) <= tolerance * std::max(1.0, std::abs(mid))) {
                result.converged = true;
                break;
            }
            if ((fmid > 0.0) == (fa > 0.0)) {
                a = mid;
                fa = fmid;
            } else {
                b = mid;
            }
        }
        result.root = mid;
        result.residual = fmid;
        return result;
    }

    SystemRootResult newton_system(const VectorFunction &f, const std::vector<double> &x0, const double tolerance,
                                   const int max_iterations) {
        require_solver_settings(tolerance, max_iterations);
        if (x0.empty()) {
            detail::throw_invalid("x0 must contain at least one coordinate");
        }

        SystemRootResult result;
        std::vector<double> x = x0;
        std::vector<double> fx = f(x);
        if (fx.size() != x.size()) {
            detail::throw_invalid("newton_system requires a square system: len(F(x)) must equal len(x)");
        }

        for (int k = 1; k <= max_iterations; ++k) {
            result.iterations = k;
            const double residual = math::norm(fx);
            if (residual == 0.0) {
                result.converged = true;
                break;
            }
            if (!std::isfinite(residual)) {
                break;
            }

            std::vector<double> neg_fx(fx.size());
            std::transform(fx.begin(), fx.end(), neg_fx.begin(), [](const double v) { return -v; });
            std::vector<double> step;
            try {
                step = math::solve_linear_system(jacobian(f, x), neg_fx);
            } catch (const std::domain_error &) {
                break; // Singular Jacobian.
            }

            for (std::size_t i = 0; i < x.size(); ++i) {
                x[i] += step[i];
            }
            fx = f(x);
            if (math::norm(step) <= tolerance * std::max(1.0, math::norm(x))) {
                result.converged = std::isfinite(math::norm(fx));
                break;
            }
        }
        result.root = x;
        result.residual_norm = math::norm(fx);
        return result;
    }

} // namespace physicsforge::numerics

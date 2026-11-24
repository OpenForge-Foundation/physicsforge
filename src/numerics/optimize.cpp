// src/numerics/optimize.cpp
#include "physicsforge/numerics/optimize.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "physicsforge/detail/validation.h"
#include "physicsforge/math/linalg.h"
#include "physicsforge/numerics/calculus.h"

namespace physicsforge::numerics {

    namespace {
        constexpr double ARMIJO_C = 1e-4;
        constexpr double MIN_STEP_SCALE = 1e-16;

        // Newton direction p solving H p = -g, or steepest descent if that is not a descent direction.
        std::vector<double> search_direction(const MultivariateFunction &f, const std::vector<double> &x,
                                             const std::vector<double> &grad) {
            std::vector<double> neg_grad(grad.size());
            std::transform(grad.begin(), grad.end(), neg_grad.begin(), [](const double g) { return -g; });
            try {
                std::vector<double> p = math::solve_linear_system(hessian(f, x), neg_grad);
                if (math::dot(grad, p) < 0.0) {
                    return p;
                }
            } catch (const std::domain_error &) {
                // Singular Hessian: fall through to steepest descent.
            }
            return neg_grad;
        }
    } // namespace

    OptimizeResult newton_minimize(const MultivariateFunction &f, const std::vector<double> &x0, const double tolerance,
                                   const int max_iterations) {
        detail::require_positive(tolerance, "tolerance");
        if (max_iterations < 1) {
            detail::throw_invalid("max_iterations", "at least 1", static_cast<double>(max_iterations));
        }
        if (x0.empty()) {
            detail::throw_invalid("x0 must contain at least one coordinate");
        }

        OptimizeResult result;
        std::vector<double> x = x0;
        double fx = f(x);
        if (!std::isfinite(fx)) {
            detail::throw_invalid("f(x0) must be finite");
        }
        std::vector<double> grad = gradient(f, x);
        std::vector<double> candidate(x.size());

        for (int k = 1; k <= max_iterations; ++k) {
            result.iterations = k;
            if (math::norm(grad) <= tolerance) {
                result.converged = true;
                break;
            }

            const std::vector<double> p = search_direction(f, x, grad);
            const double slope = math::dot(grad, p);

            // Armijo backtracking: accept t with f(x + t p) <= f(x) + c t grad.p.
            double t = 1.0;
            double f_candidate = fx;
            bool accepted = false;
            while (t >= MIN_STEP_SCALE) {
                for (std::size_t i = 0; i < x.size(); ++i) {
                    candidate[i] = x[i] + t * p[i];
                }
                f_candidate = f(candidate);
                if (std::isfinite(f_candidate) && f_candidate <= fx + ARMIJO_C * t * slope) {
                    accepted = true;
                    break;
                }
                t *= 0.5;
            }
            if (!accepted) {
                break; // No decrease possible along p: stalled.
            }

            const double step_norm = t * math::norm(p);
            x.swap(candidate);
            fx = f_candidate;
            grad = gradient(f, x);
            if (step_norm <= tolerance * std::max(1.0, math::norm(x))) {
                result.converged = true;
                break;
            }
        }
        result.x = std::move(x);
        result.value = fx;
        result.gradient = std::move(grad);
        return result;
    }

} // namespace physicsforge::numerics

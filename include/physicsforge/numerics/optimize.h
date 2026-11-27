/**
 * @file optimize.h
 * @brief Unconstrained minimisation with Newton's method on the gradient and Hessian.
 *
 * Each iteration solves \f$ H(\mathbf x_k)\, \mathbf p = -\nabla f(\mathbf x_k) \f$. Where the Hessian is
 * singular or not positive definite (so \f$ \mathbf p \f$ is not a descent direction) the step falls back
 * to steepest descent, and every step is safeguarded by Armijo backtracking. Near a strict local minimum
 * this reduces to pure Newton and converges quadratically.
 *
 * @ingroup numerics
 */

#ifndef PHYSICSFORGE_NUMERICS_OPTIMIZE_H
#define PHYSICSFORGE_NUMERICS_OPTIMIZE_H

#include <vector>

#include "physicsforge/config.h"
#include "physicsforge/numerics/calculus.h"

namespace physicsforge::numerics {

    /** @brief Outcome of a minimisation. */
    struct OptimizeResult {
        std::vector<double> x; ///< Final iterate (the minimiser on success).
        double value = 0.0; ///< \f$ f(\mathbf x) \f$.
        std::vector<double> gradient; ///< \f$ \nabla f(\mathbf x) \f$.
        int iterations = 0; ///< Iterations performed.
        bool converged = false; ///< Whether the tolerance was met.
    };

    /**
     * @brief Default gradient-norm tolerance.
     *
     * Looser than config::DEFAULT_TOLERANCE because gradients come from finite differences, which are
     * accurate to roughly \f$ \epsilon^{4/5} \f$.
     */
    inline constexpr double DEFAULT_OPTIMIZE_TOLERANCE = 1e-8;

    /**
     * @brief Minimises \f$ f : \mathbb R^n \to \mathbb R \f$ with a damped Newton method.
     * @param f Objective function.
     * @param x0 Starting point; non-empty and finite.
     * @param tolerance Converged when \f$ \|\nabla f\| \le \text{tol} \f$ or the step is below
     *        \f$ \text{tol} \max(1, \|\mathbf x\|) \f$; must be positive.
     * @param max_iterations Iteration cap; must be positive.
     * @return The final iterate, objective value, gradient and convergence information.
     * @throws std::invalid_argument If an argument is invalid.
     */
    OptimizeResult newton_minimize(const MultivariateFunction &f, const std::vector<double> &x0,
                                   double tolerance = DEFAULT_OPTIMIZE_TOLERANCE,
                                   int max_iterations = config::DEFAULT_MAX_ITERATIONS);

} // namespace physicsforge::numerics

#endif // PHYSICSFORGE_NUMERICS_OPTIMIZE_H

// include/physicsforge/numerics/optimize.h
//
// Unconstrained minimisation with Newton's method on the gradient and Hessian.
//
// Each iteration solves H(x_k) p = -grad f(x_k). Where the Hessian is singular or
// not positive definite (p is not a descent direction) the step falls back to
// steepest descent, and every step is safeguarded with an Armijo backtracking line
// search. Near a strict local minimum this reduces to pure Newton and converges
// quadratically.

#ifndef PHYSICSFORGE_NUMERICS_OPTIMIZE_H
#define PHYSICSFORGE_NUMERICS_OPTIMIZE_H

#include <vector>

#include "physicsforge/config.h"
#include "physicsforge/numerics/functions.h"

namespace physicsforge::numerics {

    // Outcome of a minimisation.
    struct OptimizeResult {
        std::vector<double> x; // Final iterate (the minimiser on success).
        double value = 0.0; // f(x).
        std::vector<double> gradient; // grad f(x).
        int iterations = 0;
        bool converged = false;
    };

    // Default gradient-norm tolerance. Looser than DEFAULT_TOLERANCE because gradients are
    // obtained by finite differences, which are accurate to roughly eps^(4/5).
    inline constexpr double DEFAULT_OPTIMIZE_TOLERANCE = 1e-8;

    // Minimises f : R^n -> R from x0 using a damped Newton method with numerical derivatives.
    // Converges when ||grad f|| <= tol or the step is below tol * max(1, ||x||).
    OptimizeResult newton_minimize(const MultivariateFunction &f, const std::vector<double> &x0,
                                   double tolerance = DEFAULT_OPTIMIZE_TOLERANCE,
                                   int max_iterations = config::DEFAULT_MAX_ITERATIONS);

} // namespace physicsforge::numerics

#endif // PHYSICSFORGE_NUMERICS_OPTIMIZE_H

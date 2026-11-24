// include/physicsforge/numerics/roots.h
//
// Root finding: Newton-Raphson (scalar and multivariate) and bisection.
//
// Solvers never throw on non-convergence; they return a result whose `converged`
// flag the caller must check (the Python layer raises ConvergenceError by default).
// Invalid arguments (non-finite start, bad bracket, ...) do throw.

#ifndef PHYSICSFORGE_NUMERICS_ROOTS_H
#define PHYSICSFORGE_NUMERICS_ROOTS_H

#include <vector>

#include "physicsforge/config.h"
#include "physicsforge/numerics/functions.h"

namespace physicsforge::numerics {

    // Outcome of a scalar root search.
    struct RootResult {
        double root = 0.0; // Final iterate.
        double residual = 0.0; // f(root).
        int iterations = 0; // Iterations performed.
        bool converged = false;
    };

    // Outcome of a multivariate root search.
    struct SystemRootResult {
        std::vector<double> root; // Final iterate.
        double residual_norm = 0.0; // ||F(root)||_2.
        int iterations = 0;
        bool converged = false;
    };

    // Newton-Raphson with an analytic derivative: x_{k+1} = x_k - f(x_k) / f'(x_k).
    // Quadratic convergence near a simple root. Stops when |dx| <= tol * max(1, |x|) or
    // f(x) == 0; reports non-convergence on a zero/non-finite derivative or iteration cap.
    RootResult newton_raphson(const ScalarFunction &f, const ScalarFunction &derivative, double x0,
                              double tolerance = config::DEFAULT_TOLERANCE,
                              int max_iterations = config::DEFAULT_MAX_ITERATIONS);

    // Newton-Raphson with a numerically differentiated derivative (five-point stencil).
    RootResult newton_raphson(const ScalarFunction &f, double x0, double tolerance = config::DEFAULT_TOLERANCE,
                              int max_iterations = config::DEFAULT_MAX_ITERATIONS);

    // Bisection on a sign-changing bracket [a, b]. Guaranteed (linear) convergence.
    // Throws if f(a) and f(b) have the same strict sign.
    RootResult bisection(const ScalarFunction &f, double a, double b, double tolerance = config::DEFAULT_TOLERANCE,
                         int max_iterations = 200);

    // Multivariate Newton's method for F(x) = 0 with F : R^n -> R^n:
    //   J(x_k) dx = -F(x_k),  x_{k+1} = x_k + dx,
    // with a numerically evaluated Jacobian.
    SystemRootResult newton_system(const VectorFunction &f, const std::vector<double> &x0,
                                   double tolerance = config::DEFAULT_TOLERANCE,
                                   int max_iterations = config::DEFAULT_MAX_ITERATIONS);

} // namespace physicsforge::numerics

#endif // PHYSICSFORGE_NUMERICS_ROOTS_H

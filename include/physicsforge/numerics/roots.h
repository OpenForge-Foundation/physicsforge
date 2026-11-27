/**
 * @file roots.h
 * @brief Root finding: Newton–Raphson (scalar and multivariate) and bisection.
 *
 * Solvers never throw on non-convergence: they return a result whose @c converged flag the caller must
 * check (the Python layer raises @c ConvergenceError by default). Invalid arguments do throw.
 *
 * @ingroup numerics
 */

#ifndef PHYSICSFORGE_NUMERICS_ROOTS_H
#define PHYSICSFORGE_NUMERICS_ROOTS_H

#include <vector>

#include "physicsforge/config.h"
#include "physicsforge/numerics/calculus.h"

namespace physicsforge::numerics {

    /** @brief Outcome of a scalar root search. */
    struct RootResult {
        double root = 0.0; ///< Final iterate.
        double residual = 0.0; ///< \f$ f(\text{root}) \f$.
        int iterations = 0; ///< Iterations performed.
        bool converged = false; ///< Whether the tolerance was met.
    };

    /** @brief Outcome of a multivariate root search. */
    struct SystemRootResult {
        std::vector<double> root; ///< Final iterate.
        double residual_norm = 0.0; ///< \f$ \|F(\text{root})\|_2 \f$.
        int iterations = 0; ///< Iterations performed.
        bool converged = false; ///< Whether the tolerance was met.
    };

    /**
     * @brief Newton–Raphson with an analytic derivative, \f$ x_{k+1} = x_k - f(x_k) / f'(x_k) \f$.
     *
     * Converges quadratically near a simple root. Stops when \f$ |\Delta x| \le \text{tol} \max(1, |x|) \f$
     * or \f$ f(x) = 0 \f$; reports non-convergence on a zero or non-finite derivative or at the iteration cap.
     *
     * @param f Function whose root is sought.
     * @param derivative Its derivative \f$ f' \f$.
     * @param x0 Starting point; must be finite.
     * @param tolerance Relative step tolerance; must be positive.
     * @param max_iterations Iteration cap; must be positive.
     * @return The final iterate and convergence information.
     * @throws std::invalid_argument If an argument is invalid.
     */
    RootResult newton_raphson(const ScalarFunction &f, const ScalarFunction &derivative, double x0,
                              double tolerance = config::DEFAULT_TOLERANCE,
                              int max_iterations = config::DEFAULT_MAX_ITERATIONS);

    /**
     * @brief Newton–Raphson with a numerically differentiated derivative (five-point stencil).
     * @param f Function whose root is sought.
     * @param x0 Starting point; must be finite.
     * @param tolerance Relative step tolerance; must be positive.
     * @param max_iterations Iteration cap; must be positive.
     * @return The final iterate and convergence information.
     * @throws std::invalid_argument If an argument is invalid.
     */
    RootResult newton_raphson(const ScalarFunction &f, double x0, double tolerance = config::DEFAULT_TOLERANCE,
                              int max_iterations = config::DEFAULT_MAX_ITERATIONS);

    /**
     * @brief Bisection on a sign-changing bracket; guaranteed linear convergence.
     * @param f Continuous function on \f$ [a, b] \f$.
     * @param a Lower bracket end.
     * @param b Upper bracket end.
     * @param tolerance Absolute bracket-width tolerance; must be positive.
     * @param max_iterations Iteration cap; must be positive.
     * @return The bracket midpoint and convergence information.
     * @throws std::invalid_argument If \f$ f(a) \f$ and \f$ f(b) \f$ have the same strict sign.
     */
    RootResult bisection(const ScalarFunction &f, double a, double b, double tolerance = config::DEFAULT_TOLERANCE,
                         int max_iterations = 200);

    /**
     * @brief Newton's method for a system \f$ F(\mathbf x) = 0 \f$ with \f$ F : \mathbb R^n \to \mathbb R^n \f$.
     *
     * Each iteration solves \f$ J(\mathbf x_k)\, \Delta\mathbf x = -F(\mathbf x_k) \f$ with a numerical Jacobian.
     *
     * @param f Vector function \f$ F \f$.
     * @param x0 Starting point; non-empty and finite.
     * @param tolerance Tolerance on \f$ \|F\|_2 \f$ and on the relative step; must be positive.
     * @param max_iterations Iteration cap; must be positive.
     * @return The final iterate and convergence information.
     * @throws std::invalid_argument If an argument is invalid or @p f is not square.
     */
    SystemRootResult newton_system(const VectorFunction &f, const std::vector<double> &x0,
                                   double tolerance = config::DEFAULT_TOLERANCE,
                                   int max_iterations = config::DEFAULT_MAX_ITERATIONS);

} // namespace physicsforge::numerics

#endif // PHYSICSFORGE_NUMERICS_ROOTS_H

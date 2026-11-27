/**
 * @file config.h
 * @brief Library-wide numerical defaults.
 *
 * Used as default arguments by physicsforge::numerics and the PDE solvers. They are deliberately
 * conservative: tight enough for double-precision physics, loose enough to terminate quickly on
 * well-conditioned problems.
 */

#ifndef PHYSICSFORGE_CONFIG_H
#define PHYSICSFORGE_CONFIG_H

namespace physicsforge::config {

    /** @brief Default tolerance for iterative solvers (root finding, optimisation). */
    inline constexpr double DEFAULT_TOLERANCE = 1e-12;

    /** @brief Default iteration cap for iterative solvers. */
    inline constexpr int DEFAULT_MAX_ITERATIONS = 100;

    /** @brief Default absolute tolerance for adaptive quadrature. */
    inline constexpr double DEFAULT_INTEGRATION_TOLERANCE = 1e-10;

    /** @brief Maximum recursion depth for adaptive quadrature. */
    inline constexpr int DEFAULT_INTEGRATION_MAX_DEPTH = 50;

    /** @brief Default number of Gaver–Stehfest terms (even; above ~18 loses accuracy in double precision). */
    inline constexpr int DEFAULT_STEHFEST_TERMS = 14;

    /** @brief Largest stable Fourier number \f$ r = \alpha \Delta t / \Delta x^2 \f$ for the explicit heat solver. */
    inline constexpr double FTCS_STABILITY_LIMIT = 0.5;

    /** @brief Smallest pivot magnitude the dense solver accepts before treating a matrix as singular. */
    inline constexpr double SINGULAR_PIVOT_THRESHOLD = 1e-300;

} // namespace physicsforge::config

#endif // PHYSICSFORGE_CONFIG_H

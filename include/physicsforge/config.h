// include/physicsforge/config.h
//
// Library-wide numerical defaults.
//
// These values are used as default arguments by the numerical routines in
// physicsforge::numerics and by the PDE solvers. They are intentionally
// conservative: tight enough for double-precision physics work, loose enough to
// terminate quickly on well-conditioned problems.

#ifndef PHYSICSFORGE_CONFIG_H
#define PHYSICSFORGE_CONFIG_H

#include <cstddef>

namespace physicsforge::config {

    // Default absolute tolerance for iterative solvers (root finding, optimisation).
    inline constexpr double DEFAULT_TOLERANCE = 1e-12;

    // Default iteration cap for iterative solvers.
    inline constexpr int DEFAULT_MAX_ITERATIONS = 100;

    // Default tolerance for adaptive quadrature.
    inline constexpr double DEFAULT_INTEGRATION_TOLERANCE = 1e-10;

    // Maximum recursion depth for adaptive quadrature.
    inline constexpr int DEFAULT_INTEGRATION_MAX_DEPTH = 50;

    // Default number of terms for the Gaver-Stehfest inverse Laplace transform.
    // Must be even. Values above ~18 lose accuracy in double precision.
    inline constexpr int DEFAULT_STEHFEST_TERMS = 14;

    // Upper bound on the explicit (FTCS) heat-equation Fourier number r = alpha*dt/dx^2.
    // The scheme is unstable for r > 0.5.
    inline constexpr double FTCS_STABILITY_LIMIT = 0.5;

    // Smallest pivot magnitude accepted by the dense linear solver before the matrix
    // is treated as singular.
    inline constexpr double SINGULAR_PIVOT_THRESHOLD = 1e-300;

} // namespace physicsforge::config

#endif // PHYSICSFORGE_CONFIG_H

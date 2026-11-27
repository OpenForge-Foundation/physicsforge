/**
 * @file integrate.h
 * @brief Numerical quadrature for sampled data and for callables.
 *
 * @ingroup numerics
 */

#ifndef PHYSICSFORGE_NUMERICS_INTEGRATE_H
#define PHYSICSFORGE_NUMERICS_INTEGRATE_H

#include <span>

#include "physicsforge/config.h"
#include "physicsforge/numerics/calculus.h"

namespace physicsforge::numerics {

    /**
     * @brief Composite trapezoidal rule over uniformly spaced samples.
     * @param samples Function values; at least 2.
     * @param spacing Sample spacing \f$ \Delta x \f$; must be positive.
     * @return Approximation of the integral.
     * @throws std::invalid_argument If there are fewer than 2 samples or @p spacing is not positive.
     */
    double trapezoid(std::span<const double> samples, double spacing);

    /**
     * @brief Composite Simpson rule over uniformly spaced samples.
     *
     * Uses Simpson's 1/3 rule and, when the number of intervals is odd, Simpson's 3/8 rule on the last
     * three intervals, so it is fourth-order for any count \f$ \ge 4 \f$. Two samples fall back to the
     * trapezoid rule and three to a single Simpson panel.
     *
     * @param samples Function values; at least 2.
     * @param spacing Sample spacing \f$ \Delta x \f$; must be positive.
     * @return Approximation of the integral.
     * @throws std::invalid_argument If there are fewer than 2 samples or @p spacing is not positive.
     */
    double simpson(std::span<const double> samples, double spacing);

    /**
     * @brief Adaptive Simpson quadrature of \f$ \int_a^b f(x)\,dx \f$ with Richardson extrapolation.
     * @param f Integrand.
     * @param a Lower limit; must be finite.
     * @param b Upper limit; must be finite (\f$ a > b \f$ yields the negated integral).
     * @param tolerance Absolute error target; must be positive.
     * @param max_depth Maximum recursion depth; must be positive.
     * @return Approximation of the integral.
     * @throws std::invalid_argument If an argument is invalid.
     */
    double integrate(const ScalarFunction &f, double a, double b,
                     double tolerance = config::DEFAULT_INTEGRATION_TOLERANCE,
                     int max_depth = config::DEFAULT_INTEGRATION_MAX_DEPTH);

} // namespace physicsforge::numerics

#endif // PHYSICSFORGE_NUMERICS_INTEGRATE_H

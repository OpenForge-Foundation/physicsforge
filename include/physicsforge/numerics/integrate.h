// include/physicsforge/numerics/integrate.h
//
// Numerical quadrature for sampled data and for callables.

#ifndef PHYSICSFORGE_NUMERICS_INTEGRATE_H
#define PHYSICSFORGE_NUMERICS_INTEGRATE_H

#include <span>

#include "physicsforge/config.h"
#include "physicsforge/numerics/functions.h"

namespace physicsforge::numerics {

    // Composite trapezoidal rule over uniformly spaced samples. Needs >= 2 samples, dx > 0.
    double trapezoid(std::span<const double> samples, double spacing);

    // Composite Simpson rule over uniformly spaced samples. Uses Simpson's 1/3 rule and,
    // when the number of intervals is odd, Simpson's 3/8 rule on the final three intervals,
    // so it is fourth-order accurate for any sample count >= 4 (2 or 3 samples fall back to
    // the trapezoid / single Simpson panel).
    double simpson(std::span<const double> samples, double spacing);

    // Adaptive Simpson quadrature of f over [a, b] to an absolute tolerance, with Richardson
    // extrapolation on accepted panels. a > b yields the negated integral.
    double integrate(const ScalarFunction &f, double a, double b,
                     double tolerance = config::DEFAULT_INTEGRATION_TOLERANCE,
                     int max_depth = config::DEFAULT_INTEGRATION_MAX_DEPTH);

} // namespace physicsforge::numerics

#endif // PHYSICSFORGE_NUMERICS_INTEGRATE_H

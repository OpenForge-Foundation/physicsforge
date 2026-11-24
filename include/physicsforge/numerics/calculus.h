// include/physicsforge/numerics/calculus.h
//
// Numerical differentiation: derivatives, gradients, Hessians, Jacobians and the
// vector-calculus operators (divergence, curl, Laplacian) used by the differential
// forms of Maxwell's equations and the heat equation.
//
// Step sizes are chosen automatically to balance truncation against round-off
// error (h ~ eps^(1/(p+1)) * max(1, |x|) for a method of order p) unless a
// positive step is supplied explicitly.

#ifndef PHYSICSFORGE_NUMERICS_CALCULUS_H
#define PHYSICSFORGE_NUMERICS_CALCULUS_H

#include <vector>

#include "physicsforge/math/vector3.h"
#include "physicsforge/numerics/functions.h"

namespace physicsforge::numerics {

    // First derivative f'(x) by the fourth-order five-point central stencil. step <= 0 selects
    // an automatic step.
    double derivative(const ScalarFunction &f, double x, double step = 0.0);

    // Second derivative f''(x) by the second-order three-point central stencil.
    double second_derivative(const ScalarFunction &f, double x, double step = 0.0);

    // Gradient of f : R^n -> R at x (five-point stencil per coordinate; 4n evaluations).
    std::vector<double> gradient(const MultivariateFunction &f, const std::vector<double> &x);

    // Hessian of f : R^n -> R at x, returned row-major (n*n) and exactly symmetric.
    std::vector<double> hessian(const MultivariateFunction &f, const std::vector<double> &x);

    // Jacobian of F : R^n -> R^m at x, returned row-major (m*n): J[i*n + j] = dF_i/dx_j.
    std::vector<double> jacobian(const VectorFunction &f, const std::vector<double> &x);

    // Laplacian (sum of unmixed second partials) of f : R^n -> R at x.
    double laplacian(const MultivariateFunction &f, const std::vector<double> &x);

    // Divergence of a vector field at point p: div F = dFx/dx + dFy/dy + dFz/dz.
    double divergence(const VectorField &field, const math::Vec3 &point);

    // Curl of a vector field at point p: curl F = nabla x F.
    math::Vec3 curl(const VectorField &field, const math::Vec3 &point);

} // namespace physicsforge::numerics

#endif // PHYSICSFORGE_NUMERICS_CALCULUS_H

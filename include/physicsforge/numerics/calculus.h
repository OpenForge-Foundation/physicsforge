/**
 * @file calculus.h
 * @brief Numerical differentiation and vector calculus, plus the callable types used across numerics.
 *
 * Derivatives, gradients, Hessians, Jacobians and the vector-calculus operators (divergence, curl,
 * Laplacian) used by the differential forms of Maxwell's equations and the heat equation.
 *
 * Step sizes are chosen automatically to balance truncation against round-off error,
 * \f$ h \approx \epsilon^{1/(p+1)} \max(1, |x|) \f$ for a stencil of order \f$ p \f$, unless a positive
 * step is supplied. Because of the absolute floor, pass an explicit step for variables much smaller
 * than 1 (e.g. wavelengths in metres).
 *
 * @ingroup numerics
 */

#ifndef PHYSICSFORGE_NUMERICS_CALCULUS_H
#define PHYSICSFORGE_NUMERICS_CALCULUS_H

#include <functional>
#include <vector>

#include "physicsforge/math/vector3.h"

/** @defgroup numerics Numerics
 *  @brief Calculus, root finding, optimisation, quadrature, ODE integration and transforms.
 */

namespace physicsforge::numerics {

    /**
     * @name Callable types
     * @c std::function is used rather than templates so the routines compile once into the library and
     * bind directly to Python callables; the type-erasure cost is negligible next to the user function.
     * @{
     */
    /** @brief \f$ f : \mathbb R \to \mathbb R \f$. */
    using ScalarFunction = std::function<double(double)>;
    /** @brief \f$ f : \mathbb R^n \to \mathbb R \f$. */
    using MultivariateFunction = std::function<double(const std::vector<double> &)>;
    /** @brief \f$ F : \mathbb R^n \to \mathbb R^m \f$. */
    using VectorFunction = std::function<std::vector<double>(const std::vector<double> &)>;
    /** @brief \f$ F : \mathbb R^3 \to \mathbb R^3 \f$, a vector field such as \f$ \mathbf E(\mathbf r) \f$. */
    using VectorField = std::function<math::Vec3(const math::Vec3 &)>;
    /** @} */

    /**
     * @brief First derivative \f$ f'(x) \f$ by the fourth-order five-point central stencil.
     * @param f Function to differentiate.
     * @param x Evaluation point; must be finite.
     * @param step Step \f$ h \f$; zero or negative selects the automatic step.
     * @return Approximation of \f$ f'(x) \f$.
     * @throws std::invalid_argument If @p x or @p step is non-finite.
     */
    double derivative(const ScalarFunction &f, double x, double step = 0.0);

    /**
     * @brief Second derivative \f$ f''(x) \f$ by the second-order three-point central stencil.
     * @param f Function to differentiate.
     * @param x Evaluation point; must be finite.
     * @param step Step \f$ h \f$; zero or negative selects the automatic step.
     * @return Approximation of \f$ f''(x) \f$.
     * @throws std::invalid_argument If @p x or @p step is non-finite.
     */
    double second_derivative(const ScalarFunction &f, double x, double step = 0.0);

    /**
     * @brief Gradient \f$ \nabla f(\mathbf x) \f$ (five-point stencil per coordinate, \f$ 4n \f$ evaluations).
     * @param f Scalar field \f$ f : \mathbb R^n \to \mathbb R \f$.
     * @param x Evaluation point; non-empty and finite.
     * @return Gradient vector of length \f$ n \f$.
     * @throws std::invalid_argument If @p x is empty or non-finite.
     */
    std::vector<double> gradient(const MultivariateFunction &f, const std::vector<double> &x);

    /**
     * @brief Hessian \f$ H_{ij} = \partial^2 f / \partial x_i \partial x_j \f$.
     * @param f Scalar field \f$ f : \mathbb R^n \to \mathbb R \f$.
     * @param x Evaluation point; non-empty and finite.
     * @return Row-major \f$ n \times n \f$ matrix, exactly symmetric.
     * @throws std::invalid_argument If @p x is empty or non-finite.
     */
    std::vector<double> hessian(const MultivariateFunction &f, const std::vector<double> &x);

    /**
     * @brief Jacobian \f$ J_{ij} = \partial F_i / \partial x_j \f$.
     * @param f Vector function \f$ F : \mathbb R^n \to \mathbb R^m \f$.
     * @param x Evaluation point; non-empty and finite.
     * @return Row-major \f$ m \times n \f$ matrix, @c J[i*n + j].
     * @throws std::invalid_argument If @p x is empty or non-finite, or @p f returns inconsistent sizes.
     */
    std::vector<double> jacobian(const VectorFunction &f, const std::vector<double> &x);

    /**
     * @brief Laplacian \f$ \nabla^2 f = \sum_i \partial^2 f / \partial x_i^2 \f$.
     * @param f Scalar field \f$ f : \mathbb R^n \to \mathbb R \f$.
     * @param x Evaluation point; non-empty and finite.
     * @return Laplacian at @p x.
     * @throws std::invalid_argument If @p x is empty or non-finite.
     */
    double laplacian(const MultivariateFunction &f, const std::vector<double> &x);

    /**
     * @brief Divergence \f$ \nabla \cdot \mathbf F \f$ of a vector field.
     * @param field Vector field \f$ \mathbf F : \mathbb R^3 \to \mathbb R^3 \f$.
     * @param point Evaluation point; must be finite.
     * @return Divergence at @p point.
     * @throws std::invalid_argument If @p point is non-finite.
     */
    double divergence(const VectorField &field, const math::Vec3 &point);

    /**
     * @brief Curl \f$ \nabla \times \mathbf F \f$ of a vector field.
     * @param field Vector field \f$ \mathbf F : \mathbb R^3 \to \mathbb R^3 \f$.
     * @param point Evaluation point; must be finite.
     * @return Curl at @p point.
     * @throws std::invalid_argument If @p point is non-finite.
     */
    math::Vec3 curl(const VectorField &field, const math::Vec3 &point);

} // namespace physicsforge::numerics

#endif // PHYSICSFORGE_NUMERICS_CALCULUS_H

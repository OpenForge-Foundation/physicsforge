/**
 * @file linalg.h
 * @brief Dense linear-algebra kernels for \f$ \mathbb R^n \f$.
 *
 * Small, dependency-free routines sized for the problems physicsforge solves internally (Newton steps on
 * the gradients, Hessians and Jacobians of a handful of variables). Matrices are row-major in a
 * contiguous buffer.
 *
 * @ingroup math
 */

#ifndef PHYSICSFORGE_MATH_LINALG_H
#define PHYSICSFORGE_MATH_LINALG_H

#include <span>
#include <vector>

namespace physicsforge::math {

    /**
     * @brief Inner product \f$ \sum_i a_i b_i \f$.
     * @param a First vector.
     * @param b Second vector; same length as @p a.
     * @return The inner product.
     * @throws std::invalid_argument If the lengths differ.
     */
    double dot(std::span<const double> a, std::span<const double> b);

    /**
     * @brief Euclidean norm \f$ \|a\|_2 \f$, scaled internally to avoid overflow and underflow.
     * @param a Vector.
     * @return The norm.
     */
    double norm(std::span<const double> a);

    /**
     * @brief Solves \f$ A\mathbf x = \mathbf b \f$ by Gaussian elimination with partial pivoting.
     *
     * O(n³) time and O(n²) extra memory.
     *
     * @param matrix \f$ A \f$, \f$ n \times n \f$ row-major (length \f$ n^2 \f$).
     * @param rhs \f$ \mathbf b \f$, length \f$ n \f$.
     * @return Solution \f$ \mathbf x \f$.
     * @throws std::invalid_argument On a shape mismatch.
     * @throws std::domain_error If \f$ A \f$ is numerically singular.
     */
    std::vector<double> solve_linear_system(std::span<const double> matrix, std::span<const double> rhs);

} // namespace physicsforge::math

#endif // PHYSICSFORGE_MATH_LINALG_H

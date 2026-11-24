// include/physicsforge/math/linalg.h
//
// Dense linear-algebra kernels for R^n.
//
// These are small, dependency-free routines sized for the problems physicsforge
// solves internally (Newton steps on gradients/Hessians/Jacobians of a handful of
// variables). Matrices are stored row-major in contiguous std::vector<double>.

#ifndef PHYSICSFORGE_MATH_LINALG_H
#define PHYSICSFORGE_MATH_LINALG_H

#include <cstddef>
#include <span>
#include <vector>

namespace physicsforge::math {

    // Inner product sum_i a_i b_i. Throws std::invalid_argument on length mismatch.
    double dot(std::span<const double> a, std::span<const double> b);

    // Euclidean (L2) norm, computed with scaling to avoid overflow/underflow.
    double norm(std::span<const double> a);

    // Solves the square system A x = b by Gaussian elimination with partial pivoting.
    //
    // A is n x n, row-major, length n*n; b has length n. Throws std::invalid_argument
    // on shape mismatch and std::domain_error if A is (numerically) singular.
    // Complexity: O(n^3) time, O(n^2) extra memory.
    std::vector<double> solve_linear_system(std::span<const double> matrix, std::span<const double> rhs);

} // namespace physicsforge::math

#endif // PHYSICSFORGE_MATH_LINALG_H

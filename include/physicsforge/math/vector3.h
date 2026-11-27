/**
 * @file vector3.h
 * @brief Three-dimensional Cartesian vector type with dot and cross products.
 *
 * Vec3 is a trivially copyable aggregate of three doubles. Everything here is header-only and inline so
 * compilers can keep vectors in registers; an out-of-line call would cost more than the arithmetic.
 *
 * @ingroup math
 */

#ifndef PHYSICSFORGE_MATH_VECTOR3_H
#define PHYSICSFORGE_MATH_VECTOR3_H

#include <cmath>
#include <stdexcept>

/** @defgroup math Math
 *  @brief Vectors and linear algebra.
 */

namespace physicsforge::math {

    /** @brief A vector in \f$ \mathbb R^3 \f$ with Cartesian components. */
    struct Vec3 {
        double x = 0.0; ///< x component.
        double y = 0.0; ///< y component.
        double z = 0.0; ///< z component.

        constexpr Vec3 &operator+=(const Vec3 &o) noexcept {
            x += o.x;
            y += o.y;
            z += o.z;
            return *this;
        }
        constexpr Vec3 &operator-=(const Vec3 &o) noexcept {
            x -= o.x;
            y -= o.y;
            z -= o.z;
            return *this;
        }
        constexpr Vec3 &operator*=(const double s) noexcept {
            x *= s;
            y *= s;
            z *= s;
            return *this;
        }
        constexpr Vec3 &operator/=(const double s) noexcept {
            x /= s;
            y /= s;
            z /= s;
            return *this;
        }

        friend constexpr bool operator==(const Vec3 &, const Vec3 &) noexcept = default;
    };

    /** @name Arithmetic operators
     *  Component-wise addition and subtraction, negation, and scaling by a scalar.
     *  @{ */
    constexpr Vec3 operator+(Vec3 a, const Vec3 &b) noexcept { return a += b; }
    constexpr Vec3 operator-(Vec3 a, const Vec3 &b) noexcept { return a -= b; }
    constexpr Vec3 operator-(const Vec3 &a) noexcept { return {-a.x, -a.y, -a.z}; }
    constexpr Vec3 operator*(Vec3 a, const double s) noexcept { return a *= s; }
    constexpr Vec3 operator*(const double s, Vec3 a) noexcept { return a *= s; }
    constexpr Vec3 operator/(Vec3 a, const double s) noexcept { return a /= s; }
    /** @} */

    /** @brief Dot product, \f$ \mathbf a \cdot \mathbf b = |\mathbf a||\mathbf b| \cos\theta \f$. */
    constexpr double dot(const Vec3 &a, const Vec3 &b) noexcept { return a.x * b.x + a.y * b.y + a.z * b.z; }

    /** @brief Cross product \f$ \mathbf a \times \mathbf b \f$, perpendicular to both, with magnitude
     *         \f$ |\mathbf a||\mathbf b| \sin\theta \f$. */
    constexpr Vec3 cross(const Vec3 &a, const Vec3 &b) noexcept {
        return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
    }

    /** @brief Squared norm \f$ |\mathbf a|^2 \f$; cheaper than norm() when only comparing lengths. */
    constexpr double norm_squared(const Vec3 &a) noexcept { return dot(a, a); }

    /** @brief Euclidean norm \f$ |\mathbf a| \f$, computed with @c hypot to avoid intermediate overflow. */
    inline double norm(const Vec3 &a) noexcept { return std::hypot(a.x, a.y, a.z); }

    /**
     * @brief Unit vector in the direction of @p a.
     * @throws std::invalid_argument If @p a is the zero vector.
     */
    inline Vec3 normalized(const Vec3 &a) {
        const double n = norm(a);
        if (n == 0.0) {
            throw std::invalid_argument("cannot normalize the zero vector");
        }
        return a / n;
    }

    /**
     * @brief Angle between two vectors in \f$ [0, \pi] \f$ [rad].
     *
     * Computed as \f$ \operatorname{atan2}(|\mathbf a \times \mathbf b|, \mathbf a \cdot \mathbf b) \f$, which,
     * unlike @c acos, stays accurate for nearly parallel and antiparallel vectors.
     *
     * @throws std::invalid_argument If either vector is zero.
     */
    inline double angle_between(const Vec3 &a, const Vec3 &b) {
        if (norm_squared(a) == 0.0 || norm_squared(b) == 0.0) {
            throw std::invalid_argument("angle is undefined for a zero vector");
        }
        return std::atan2(norm(cross(a, b)), dot(a, b));
    }

    /**
     * @brief Vector projection of @p a onto @p b, \f$ (\mathbf a \cdot \mathbf b / |\mathbf b|^2)\, \mathbf b \f$.
     * @throws std::invalid_argument If @p b is the zero vector.
     */
    inline Vec3 project(const Vec3 &a, const Vec3 &b) {
        const double denom = norm_squared(b);
        if (denom == 0.0) {
            throw std::invalid_argument("cannot project onto the zero vector");
        }
        return b * (dot(a, b) / denom);
    }

    /** @brief Scalar triple product \f$ \mathbf a \cdot (\mathbf b \times \mathbf c) \f$, the signed volume
     *         of the parallelepiped they span. */
    constexpr double triple_product(const Vec3 &a, const Vec3 &b, const Vec3 &c) noexcept {
        return dot(a, cross(b, c));
    }

    /** @brief Euclidean distance between points @p a and @p b. */
    inline double distance(const Vec3 &a, const Vec3 &b) noexcept { return norm(a - b); }

} // namespace physicsforge::math

#endif // PHYSICSFORGE_MATH_VECTOR3_H

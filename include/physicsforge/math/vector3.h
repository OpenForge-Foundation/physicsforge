// include/physicsforge/math/vector3.h
//
// Three-dimensional Cartesian vector type and its core operations.
//
// Vec3 is a trivially-copyable aggregate of three doubles. All operations are
// header-only and inline so that compilers can keep vectors in registers; this is
// the one deliberately header-only component of the core library, because the
// cost of an out-of-line call would dominate the arithmetic.

#ifndef PHYSICSFORGE_MATH_VECTOR3_H
#define PHYSICSFORGE_MATH_VECTOR3_H

#include <cmath>
#include <stdexcept>

namespace physicsforge::math {

    // A vector in R^3 with Cartesian components (x, y, z).
    struct Vec3 {
        double x = 0.0;
        double y = 0.0;
        double z = 0.0;

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

    constexpr Vec3 operator+(Vec3 a, const Vec3 &b) noexcept { return a += b; }
    constexpr Vec3 operator-(Vec3 a, const Vec3 &b) noexcept { return a -= b; }
    constexpr Vec3 operator-(const Vec3 &a) noexcept { return {-a.x, -a.y, -a.z}; }
    constexpr Vec3 operator*(Vec3 a, const double s) noexcept { return a *= s; }
    constexpr Vec3 operator*(const double s, Vec3 a) noexcept { return a *= s; }
    constexpr Vec3 operator/(Vec3 a, const double s) noexcept { return a /= s; }

    // Dot (scalar) product: a . b = |a||b| cos(theta).
    constexpr double dot(const Vec3 &a, const Vec3 &b) noexcept { return a.x * b.x + a.y * b.y + a.z * b.z; }

    // Cross (vector) product: a x b, perpendicular to both with |a x b| = |a||b| sin(theta).
    constexpr Vec3 cross(const Vec3 &a, const Vec3 &b) noexcept {
        return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
    }

    // Squared Euclidean norm, |a|^2. Cheaper than norm() when only comparisons are needed.
    constexpr double norm_squared(const Vec3 &a) noexcept { return dot(a, a); }

    // Euclidean norm (magnitude), |a|. Uses hypot to avoid intermediate overflow.
    inline double norm(const Vec3 &a) noexcept { return std::hypot(a.x, a.y, a.z); }

    // Unit vector in the direction of a. Throws std::invalid_argument for the zero vector.
    inline Vec3 normalized(const Vec3 &a) {
        const double n = norm(a);
        if (n == 0.0) {
            throw std::invalid_argument("cannot normalize the zero vector");
        }
        return a / n;
    }

    // Angle between a and b in radians, in [0, pi]. Uses atan2(|a x b|, a . b), which is
    // accurate for nearly parallel and nearly antiparallel vectors (unlike acos).
    inline double angle_between(const Vec3 &a, const Vec3 &b) {
        if (norm_squared(a) == 0.0 || norm_squared(b) == 0.0) {
            throw std::invalid_argument("angle is undefined for a zero vector");
        }
        return std::atan2(norm(cross(a, b)), dot(a, b));
    }

    // Vector projection of a onto b: (a . b / |b|^2) b.
    inline Vec3 project(const Vec3 &a, const Vec3 &b) {
        const double denom = norm_squared(b);
        if (denom == 0.0) {
            throw std::invalid_argument("cannot project onto the zero vector");
        }
        return b * (dot(a, b) / denom);
    }

    // Scalar triple product a . (b x c): the signed volume of the parallelepiped.
    constexpr double triple_product(const Vec3 &a, const Vec3 &b, const Vec3 &c) noexcept {
        return dot(a, cross(b, c));
    }

    // Euclidean distance between points a and b.
    inline double distance(const Vec3 &a, const Vec3 &b) noexcept { return norm(a - b); }

} // namespace physicsforge::math

#endif // PHYSICSFORGE_MATH_VECTOR3_H

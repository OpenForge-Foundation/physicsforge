// tests/test_vector.cpp
#include <stdexcept>
#include <vector>

#include "physicsforge/constants.h"
#include "physicsforge/math/linalg.h"
#include "physicsforge/math/vector3.h"
#include "test_framework.h"

void test_vector() {
    using namespace physicsforge::math;

    constexpr Vec3 i{1, 0, 0};
    constexpr Vec3 j{0, 1, 0};
    constexpr Vec3 k{0, 0, 1};

    // Compile-time evaluation of the core operations.
    static_assert(dot(i, j) == 0.0);
    static_assert(cross(i, j) == k);
    static_assert(cross(j, k) == i);
    static_assert(cross(k, i) == j);
    static_assert(cross(j, i) == -k);
    static_assert(triple_product(i, j, k) == 1.0);

    const Vec3 a{1, 2, 3};
    const Vec3 b{4, 5, 6};
    PF_CHECK(dot(a, b) == 32.0);
    PF_CHECK((cross(a, b) == Vec3{-3, 6, -3}));
    PF_CHECK(dot(cross(a, b), a) == 0.0); // a x b is orthogonal to a
    PF_CHECK_CLOSE(norm(Vec3{3, 4, 12}), 13.0);
    PF_CHECK_CLOSE(norm(normalized(a)), 1.0);
    PF_CHECK_CLOSE(angle_between(i, j), physicsforge::constants::PI / 2.0);
    PF_CHECK_CLOSE(angle_between(i, -i), physicsforge::constants::PI);
    PF_CHECK((project(a, i) == Vec3{1, 0, 0}));
    PF_CHECK_CLOSE(distance(a, b), std::sqrt(27.0));
    PF_CHECK_THROWS(normalized(Vec3{}), std::invalid_argument);
    PF_CHECK_THROWS(angle_between(Vec3{}, i), std::invalid_argument);

    // n-dimensional kernels
    const std::vector<double> x{3.0, 4.0};
    PF_CHECK_CLOSE(norm(x), 5.0);
    PF_CHECK_CLOSE(norm(std::vector<double>{1e200, 1e200}), std::sqrt(2.0) * 1e200); // no overflow
    PF_CHECK_THROWS(dot(x, std::vector<double>{1.0}), std::invalid_argument);

    // 3x3 system with a solution of (1, 2, 3); the zero leading pivot forces row exchange.
    const std::vector<double> A{0, 2, 1, 1, 1, 1, 2, 1, 0};
    const std::vector<double> rhs{7, 6, 4};
    const auto sol = solve_linear_system(A, rhs);
    PF_CHECK_CLOSE(sol[0], 1.0);
    PF_CHECK_CLOSE(sol[1], 2.0);
    PF_CHECK_CLOSE(sol[2], 3.0);
    PF_CHECK_THROWS(solve_linear_system(std::vector<double>{1, 2, 2, 4}, std::vector<double>{1, 2}), std::domain_error);
}

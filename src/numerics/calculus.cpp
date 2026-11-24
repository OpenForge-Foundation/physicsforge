// src/numerics/calculus.cpp
#include "physicsforge/numerics/calculus.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>

#include "physicsforge/detail/validation.h"

namespace physicsforge::numerics {

    namespace {
        constexpr double EPS = std::numeric_limits<double>::epsilon();

        // Optimal relative steps: eps^(1/5) for the 4th-order first-derivative stencil,
        // eps^(1/4) for the 2nd-order second-derivative stencil.
        const double FIRST_ORDER_STEP = std::pow(EPS, 1.0 / 5.0);
        const double SECOND_ORDER_STEP = std::pow(EPS, 1.0 / 4.0);

        double scaled_step(const double base, const double x) { return base * std::max(1.0, std::abs(x)); }

        // Five-point stencil given the four samples f(x-2h), f(x-h), f(x+h), f(x+2h).
        double five_point(const double fm2, const double fm1, const double fp1, const double fp2, const double h) {
            return (fm2 - 8.0 * fm1 + 8.0 * fp1 - fp2) / (12.0 * h);
        }

        void require_non_empty(const std::vector<double> &x) {
            if (x.empty()) {
                detail::throw_invalid("x must contain at least one coordinate");
            }
            for (const double v: x) {
                detail::require_finite(v, "x");
            }
        }

        // 3x3 Jacobian of a vector field: J[i][j] = dF_i / dx_j.
        std::array<std::array<double, 3>, 3> field_jacobian(const VectorField &field, const math::Vec3 &p) {
            detail::require_finite(p.x, "point.x");
            detail::require_finite(p.y, "point.y");
            detail::require_finite(p.z, "point.z");
            std::array<std::array<double, 3>, 3> jac{};
            const std::array<double, 3> coords{p.x, p.y, p.z};
            for (std::size_t j = 0; j < 3; ++j) {
                const double h = scaled_step(FIRST_ORDER_STEP, coords[j]);
                auto shifted = [&](const double delta) {
                    math::Vec3 q = p;
                    (j == 0 ? q.x : j == 1 ? q.y : q.z) += delta;
                    return field(q);
                };
                const math::Vec3 fm2 = shifted(-2.0 * h);
                const math::Vec3 fm1 = shifted(-h);
                const math::Vec3 fp1 = shifted(h);
                const math::Vec3 fp2 = shifted(2.0 * h);
                jac[0][j] = five_point(fm2.x, fm1.x, fp1.x, fp2.x, h);
                jac[1][j] = five_point(fm2.y, fm1.y, fp1.y, fp2.y, h);
                jac[2][j] = five_point(fm2.z, fm1.z, fp1.z, fp2.z, h);
            }
            return jac;
        }
    } // namespace

    double derivative(const ScalarFunction &f, const double x, const double step) {
        detail::require_finite(x, "x");
        detail::require_finite(step, "step");
        const double h = step > 0.0 ? step : scaled_step(FIRST_ORDER_STEP, x);
        return five_point(f(x - 2.0 * h), f(x - h), f(x + h), f(x + 2.0 * h), h);
    }

    double second_derivative(const ScalarFunction &f, const double x, const double step) {
        detail::require_finite(x, "x");
        detail::require_finite(step, "step");
        const double h = step > 0.0 ? step : scaled_step(SECOND_ORDER_STEP, x);
        return (f(x + h) - 2.0 * f(x) + f(x - h)) / (h * h);
    }

    std::vector<double> gradient(const MultivariateFunction &f, const std::vector<double> &x) {
        require_non_empty(x);
        const std::size_t n = x.size();
        std::vector<double> grad(n);
        std::vector<double> probe = x;
        for (std::size_t i = 0; i < n; ++i) {
            const double h = scaled_step(FIRST_ORDER_STEP, x[i]);
            probe[i] = x[i] - 2.0 * h;
            const double fm2 = f(probe);
            probe[i] = x[i] - h;
            const double fm1 = f(probe);
            probe[i] = x[i] + h;
            const double fp1 = f(probe);
            probe[i] = x[i] + 2.0 * h;
            const double fp2 = f(probe);
            probe[i] = x[i];
            grad[i] = five_point(fm2, fm1, fp1, fp2, h);
        }
        return grad;
    }

    std::vector<double> hessian(const MultivariateFunction &f, const std::vector<double> &x) {
        require_non_empty(x);
        const std::size_t n = x.size();
        std::vector<double> hess(n * n);
        std::vector<double> h(n);
        for (std::size_t i = 0; i < n; ++i) {
            h[i] = scaled_step(SECOND_ORDER_STEP, x[i]);
        }

        std::vector<double> probe = x;
        const double f0 = f(x);
        for (std::size_t i = 0; i < n; ++i) {
            // Diagonal: standard second difference.
            probe[i] = x[i] + h[i];
            const double fp = f(probe);
            probe[i] = x[i] - h[i];
            const double fm = f(probe);
            probe[i] = x[i];
            hess[i * n + i] = (fp - 2.0 * f0 + fm) / (h[i] * h[i]);

            // Off-diagonal: four-point mixed central difference, mirrored for exact symmetry.
            for (std::size_t j = i + 1; j < n; ++j) {
                auto eval = [&](const double si, const double sj) {
                    probe[i] = x[i] + si * h[i];
                    probe[j] = x[j] + sj * h[j];
                    const double value = f(probe);
                    probe[i] = x[i];
                    probe[j] = x[j];
                    return value;
                };
                const double mixed =
                        (eval(1.0, 1.0) - eval(1.0, -1.0) - eval(-1.0, 1.0) + eval(-1.0, -1.0)) / (4.0 * h[i] * h[j]);
                hess[i * n + j] = mixed;
                hess[j * n + i] = mixed;
            }
        }
        return hess;
    }

    std::vector<double> jacobian(const VectorFunction &f, const std::vector<double> &x) {
        require_non_empty(x);
        const std::size_t n = x.size();
        std::vector<double> probe = x;
        std::vector<double> jac;
        std::size_t m = 0;
        for (std::size_t j = 0; j < n; ++j) {
            const double h = scaled_step(FIRST_ORDER_STEP, x[j]);
            probe[j] = x[j] - 2.0 * h;
            const std::vector<double> fm2 = f(probe);
            probe[j] = x[j] - h;
            const std::vector<double> fm1 = f(probe);
            probe[j] = x[j] + h;
            const std::vector<double> fp1 = f(probe);
            probe[j] = x[j] + 2.0 * h;
            const std::vector<double> fp2 = f(probe);
            probe[j] = x[j];

            if (j == 0) {
                m = fm2.size();
                if (m == 0) {
                    detail::throw_invalid("vector function must return at least one component");
                }
                jac.assign(m * n, 0.0);
            }
            if (fm2.size() != m || fm1.size() != m || fp1.size() != m || fp2.size() != m) {
                detail::throw_invalid("vector function must return the same number of components on every call");
            }
            for (std::size_t i = 0; i < m; ++i) {
                jac[i * n + j] = five_point(fm2[i], fm1[i], fp1[i], fp2[i], h);
            }
        }
        return jac;
    }

    double laplacian(const MultivariateFunction &f, const std::vector<double> &x) {
        require_non_empty(x);
        const double f0 = f(x);
        std::vector<double> probe = x;
        double sum = 0.0;
        for (std::size_t i = 0; i < x.size(); ++i) {
            const double h = scaled_step(SECOND_ORDER_STEP, x[i]);
            probe[i] = x[i] + h;
            const double fp = f(probe);
            probe[i] = x[i] - h;
            const double fm = f(probe);
            probe[i] = x[i];
            sum += (fp - 2.0 * f0 + fm) / (h * h);
        }
        return sum;
    }

    double divergence(const VectorField &field, const math::Vec3 &point) {
        const auto j = field_jacobian(field, point);
        return j[0][0] + j[1][1] + j[2][2];
    }

    math::Vec3 curl(const VectorField &field, const math::Vec3 &point) {
        const auto j = field_jacobian(field, point);
        return {j[2][1] - j[1][2], j[0][2] - j[2][0], j[1][0] - j[0][1]};
    }

} // namespace physicsforge::numerics

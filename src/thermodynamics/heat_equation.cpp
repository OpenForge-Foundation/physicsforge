// src/thermodynamics/heat_equation.cpp
#include "physicsforge/thermodynamics/heat_equation.h"

#include <cmath>
#include <cstddef>
#include <string>
#include <utility>

#include "physicsforge/config.h"
#include "physicsforge/constants.h"
#include "physicsforge/detail/validation.h"

namespace physicsforge::thermodynamics {

    using detail::require_finite;
    using detail::require_non_negative;
    using detail::require_positive;

    namespace {
        void require_grid(const std::span<const double> temperatures) {
            if (temperatures.size() < 3) {
                detail::throw_invalid("initial_temperatures must contain at least 3 grid points");
            }
            for (const double t: temperatures) {
                require_finite(t, "initial_temperatures");
            }
        }

        void require_steps(const int steps) {
            if (steps < 0) {
                detail::throw_invalid("steps", "non-negative", static_cast<double>(steps));
            }
        }

        void require_boundaries(const double left, const double right) {
            require_finite(left, "left_temperature");
            require_finite(right, "right_temperature");
        }
    } // namespace

    double fourier_number(const double diffusivity, const double time_step, const double grid_spacing) {
        require_non_negative(diffusivity, "diffusivity");
        require_non_negative(time_step, "time_step");
        require_positive(grid_spacing, "grid_spacing");
        return diffusivity * time_step / (grid_spacing * grid_spacing);
    }

    std::vector<double> heat_fourier_coefficients(const std::span<const double> initial_temperatures, const int n_terms,
                                                  const double left_temperature, const double right_temperature) {
        require_grid(initial_temperatures);
        require_boundaries(left_temperature, right_temperature);
        if (n_terms < 1) {
            detail::throw_invalid("n_terms", "at least 1", static_cast<double>(n_terms));
        }

        // Work in normalised coordinate xi = x / L in [0, 1]; then b_n = 2 int_0^1 u(xi) sin(n pi xi) dxi,
        // independent of L. The end-point terms vanish because sin(0) = sin(n pi) = 0, so the
        // trapezoidal rule reduces to h * sum over interior points.
        const std::size_t intervals = initial_temperatures.size() - 1;
        const double h = 1.0 / static_cast<double>(intervals);

        std::vector<double> transient(intervals - 1);
        for (std::size_t j = 1; j < intervals; ++j) {
            const double xi = static_cast<double>(j) * h;
            const double steady = left_temperature + (right_temperature - left_temperature) * xi;
            transient[j - 1] = initial_temperatures[j] - steady;
        }

        std::vector<double> coefficients(static_cast<std::size_t>(n_terms));
        for (int n = 1; n <= n_terms; ++n) {
            const double k = static_cast<double>(n) * constants::PI * h;
            double sum = 0.0;
            for (std::size_t j = 1; j < intervals; ++j) {
                sum += transient[j - 1] * std::sin(k * static_cast<double>(j));
            }
            coefficients[static_cast<std::size_t>(n - 1)] = 2.0 * h * sum;
        }
        return coefficients;
    }

    std::vector<double> heat_fourier_solution(const std::span<const double> coefficients, const double length,
                                              const double diffusivity, const double time,
                                              const std::span<const double> positions, const double left_temperature,
                                              const double right_temperature) {
        require_positive(length, "length");
        require_non_negative(diffusivity, "diffusivity");
        require_non_negative(time, "time");
        require_boundaries(left_temperature, right_temperature);
        for (const double b: coefficients) {
            require_finite(b, "coefficients");
        }

        // Pre-compute each mode's wavenumber and temporal decay once: O(n_terms) exps, not O(n * m).
        const std::size_t n_terms = coefficients.size();
        std::vector<double> wavenumber(n_terms);
        std::vector<double> amplitude(n_terms);
        for (std::size_t i = 0; i < n_terms; ++i) {
            const double k = static_cast<double>(i + 1) * constants::PI / length;
            wavenumber[i] = k;
            amplitude[i] = coefficients[i] * std::exp(-diffusivity * k * k * time);
        }

        std::vector<double> result(positions.size());
        for (std::size_t p = 0; p < positions.size(); ++p) {
            const double x = positions[p];
            detail::require_in_range(x, 0.0, length, "position");
            double value = left_temperature + (right_temperature - left_temperature) * (x / length);
            for (std::size_t i = 0; i < n_terms; ++i) {
                value += amplitude[i] * std::sin(wavenumber[i] * x);
            }
            result[p] = value;
        }
        return result;
    }

    std::vector<double> heat_equation_explicit(const std::span<const double> initial_temperatures,
                                               const double diffusivity, const double grid_spacing,
                                               const double time_step, const int steps, const double left_temperature,
                                               const double right_temperature) {
        require_grid(initial_temperatures);
        require_boundaries(left_temperature, right_temperature);
        require_steps(steps);
        const double r = fourier_number(diffusivity, time_step, grid_spacing);
        if (r > config::FTCS_STABILITY_LIMIT) {
            detail::throw_invalid("Fourier number alpha*dt/dx^2 = " + std::to_string(r) +
                                  " exceeds the FTCS stability limit 0.5; reduce time_step or use "
                                  "heat_equation_crank_nicolson");
        }

        const std::size_t n = initial_temperatures.size();
        std::vector<double> current(initial_temperatures.begin(), initial_temperatures.end());
        current.front() = left_temperature;
        current.back() = right_temperature;
        std::vector<double> next = current;

        for (int step = 0; step < steps; ++step) {
            for (std::size_t i = 1; i + 1 < n; ++i) {
                next[i] = current[i] + r * (current[i - 1] - 2.0 * current[i] + current[i + 1]);
            }
            std::swap(current, next);
        }
        return current;
    }

    std::vector<double> heat_equation_crank_nicolson(const std::span<const double> initial_temperatures,
                                                     const double diffusivity, const double grid_spacing,
                                                     const double time_step, const int steps,
                                                     const double left_temperature, const double right_temperature) {
        require_grid(initial_temperatures);
        require_boundaries(left_temperature, right_temperature);
        require_steps(steps);
        const double r = fourier_number(diffusivity, time_step, grid_spacing);

        const std::size_t n = initial_temperatures.size();
        const std::size_t m = n - 2; // interior unknowns
        std::vector<double> u(initial_temperatures.begin(), initial_temperatures.end());
        u.front() = left_temperature;
        u.back() = right_temperature;

        // Implicit operator (I - r/2 D2) is tridiagonal with constant coefficients:
        //   sub = super = -r/2, diag = 1 + r.
        // Its Thomas-algorithm forward-elimination factors depend only on r, so compute them
        // once and reuse them every step: each step is then a single O(m) sweep.
        const double off = -0.5 * r;
        const double diag = 1.0 + r;
        std::vector<double> c_prime(m);
        std::vector<double> inv_denominator(m);
        inv_denominator[0] = 1.0 / diag;
        c_prime[0] = off * inv_denominator[0];
        for (std::size_t i = 1; i < m; ++i) {
            inv_denominator[i] = 1.0 / (diag - off * c_prime[i - 1]);
            c_prime[i] = off * inv_denominator[i];
        }

        std::vector<double> rhs(m);
        const double half_r = 0.5 * r;
        for (int step = 0; step < steps; ++step) {
            // Explicit half: (I + r/2 D2) u^n, with the fixed boundary values contributing from
            // both the explicit and implicit sides.
            for (std::size_t i = 0; i < m; ++i) {
                const std::size_t g = i + 1;
                rhs[i] = half_r * u[g - 1] + (1.0 - r) * u[g] + half_r * u[g + 1];
            }
            rhs[0] += half_r * left_temperature;
            rhs[m - 1] += half_r * right_temperature;

            // Thomas forward sweep (with precomputed factors) and back substitution.
            rhs[0] *= inv_denominator[0];
            for (std::size_t i = 1; i < m; ++i) {
                rhs[i] = (rhs[i] - off * rhs[i - 1]) * inv_denominator[i];
            }
            for (std::size_t i = m - 1; i-- > 0;) {
                rhs[i] -= c_prime[i] * rhs[i + 1];
            }
            for (std::size_t i = 0; i < m; ++i) {
                u[i + 1] = rhs[i];
            }
        }
        return u;
    }

} // namespace physicsforge::thermodynamics

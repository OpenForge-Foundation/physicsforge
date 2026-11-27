// src/numerics/ode.cpp
#include "physicsforge/numerics/ode.h"

#include <cstddef>
#include <string>
#include <vector>

#include "physicsforge/validation.h"

namespace physicsforge::numerics {

    namespace {

        void require_state(const std::vector<double> &state, const char *name) {
            if (state.empty()) {
                detail::throw_invalid(std::string(name) + " must contain at least one value");
            }
            for (const double v: state) {
                detail::require_finite(v, name);
            }
        }

        // Validates the shared time arguments and returns the step size h.
        double step_size(const double t0, const double t_end, const int steps) {
            detail::require_finite(t0, "t0");
            detail::require_finite(t_end, "t_end");
            if (!(t_end > t0)) {
                detail::throw_invalid("t_end must be greater than t0");
            }
            if (steps < 1) {
                detail::throw_invalid("steps", "at least 1", steps);
            }
            return (t_end - t0) / steps;
        }

        // Calls a user right-hand side and checks the result has the expected dimension.
        template<typename Fn>
        std::vector<double> evaluate(const Fn &f, const double t, const std::vector<double> &y, const char *name) {
            std::vector<double> out = f(t, y);
            if (out.size() != y.size()) {
                detail::throw_invalid(std::string(name) + " returned " + std::to_string(out.size()) +
                                      " values; expected " + std::to_string(y.size()));
            }
            return out;
        }

        std::vector<double> time_grid(const double t0, const double t_end, const double h, const int steps) {
            std::vector<double> t(static_cast<std::size_t>(steps) + 1);
            for (int i = 0; i < steps; ++i) {
                t[static_cast<std::size_t>(i)] = t0 + i * h;
            }
            t.back() = t_end; // land exactly on t_end, free of accumulated rounding
            return t;
        }

        // Shared driver for one-step methods: advance(t, h, y) returns y_{n+1}.
        template<typename Advance>
        OdeSolution integrate_first_order(const double t0, const std::vector<double> &y0, const double t_end,
                                          const int steps, Advance advance) {
            require_state(y0, "y0");
            const double h = step_size(t0, t_end, steps);
            const std::size_t n = y0.size();

            OdeSolution solution;
            solution.dimension = n;
            solution.t = time_grid(t0, t_end, h, steps);
            solution.y.reserve((static_cast<std::size_t>(steps) + 1) * n);
            solution.y.insert(solution.y.end(), y0.begin(), y0.end());

            std::vector<double> y = y0;
            for (int i = 0; i < steps; ++i) {
                y = advance(solution.t[static_cast<std::size_t>(i)], h, y);
                solution.y.insert(solution.y.end(), y.begin(), y.end());
            }
            return solution;
        }

    } // namespace

    OdeSolution euler(const OdeFunction &f, const double t0, const std::vector<double> &y0, const double t_end,
                      const int steps) {
        return integrate_first_order(t0, y0, t_end, steps, [&](const double t, const double h, std::vector<double> y) {
            const std::vector<double> k = evaluate(f, t, y, "f");
            for (std::size_t j = 0; j < y.size(); ++j) {
                y[j] += h * k[j];
            }
            return y;
        });
    }

    OdeSolution runge_kutta4(const OdeFunction &f, const double t0, const std::vector<double> &y0, const double t_end,
                             const int steps) {
        const std::size_t n = y0.size();
        std::vector<double> stage(n);
        return integrate_first_order(t0, y0, t_end, steps, [&](const double t, const double h, std::vector<double> y) {
            const auto shifted = [&](const std::vector<double> &k, const double scale) -> const std::vector<double> & {
                for (std::size_t j = 0; j < n; ++j) {
                    stage[j] = y[j] + scale * k[j];
                }
                return stage;
            };
            const std::vector<double> k1 = evaluate(f, t, y, "f");
            const std::vector<double> k2 = evaluate(f, t + 0.5 * h, shifted(k1, 0.5 * h), "f");
            const std::vector<double> k3 = evaluate(f, t + 0.5 * h, shifted(k2, 0.5 * h), "f");
            const std::vector<double> k4 = evaluate(f, t + h, shifted(k3, h), "f");
            for (std::size_t j = 0; j < n; ++j) {
                y[j] += h / 6.0 * (k1[j] + 2.0 * k2[j] + 2.0 * k3[j] + k4[j]);
            }
            return y;
        });
    }

    SecondOrderSolution semi_implicit_euler(const AccelerationFunction &acceleration, const double t0,
                                            const std::vector<double> &x0, const std::vector<double> &v0,
                                            const double t_end, const int steps) {
        require_state(x0, "x0");
        require_state(v0, "v0");
        if (v0.size() != x0.size()) {
            detail::throw_invalid("v0 must have the same length as x0");
        }
        const double h = step_size(t0, t_end, steps);
        const std::size_t n = x0.size();
        const std::size_t rows = static_cast<std::size_t>(steps) + 1;

        SecondOrderSolution solution;
        solution.dimension = n;
        solution.t = time_grid(t0, t_end, h, steps);
        solution.position.reserve(rows * n);
        solution.velocity.reserve(rows * n);
        solution.position.insert(solution.position.end(), x0.begin(), x0.end());
        solution.velocity.insert(solution.velocity.end(), v0.begin(), v0.end());

        std::vector<double> x = x0;
        std::vector<double> v = v0;
        for (int i = 0; i < steps; ++i) {
            const std::vector<double> a =
                    evaluate(acceleration, solution.t[static_cast<std::size_t>(i)], x, "acceleration");
            for (std::size_t j = 0; j < n; ++j) {
                v[j] += h * a[j];
                x[j] += h * v[j]; // the new velocity: this ordering is what makes the method symplectic
            }
            solution.position.insert(solution.position.end(), x.begin(), x.end());
            solution.velocity.insert(solution.velocity.end(), v.begin(), v.end());
        }
        return solution;
    }

} // namespace physicsforge::numerics

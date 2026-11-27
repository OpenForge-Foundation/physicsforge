/**
 * @file ode.h
 * @brief Fixed-step integrators for ordinary differential equations.
 *
 * | Method               | Order | Evaluations/step | Use it for                                                  |
 * |----------------------|-------|------------------|-------------------------------------------------------------|
 * | euler()              | 1     | 1                | Teaching and quick sketches; energy drifts in oscillators.  |
 * | semi_implicit_euler()| 1     | 1                | Mechanics \f$ \ddot x = a(t, x) \f$; symplectic, so energy stays
 * bounded. | | runge_kutta4()       | 4     | 4                | General smooth problems where accuracy matters. |
 *
 * All three take @p steps equal steps of size \f$ h = (t_{end} - t_0) / \text{steps} \f$ and return every
 * intermediate state, so the trajectory can be plotted or compared directly.
 *
 * @ingroup numerics
 */

#ifndef PHYSICSFORGE_NUMERICS_ODE_H
#define PHYSICSFORGE_NUMERICS_ODE_H

#include <cstddef>
#include <functional>
#include <vector>

namespace physicsforge::numerics {

    /** @brief Right-hand side of a first-order system, \f$ \mathbf y' = \mathbf f(t, \mathbf y) \f$. */
    using OdeFunction = std::function<std::vector<double>(double, const std::vector<double> &)>;

    /** @brief Acceleration of a second-order system, \f$ \ddot{\mathbf x} = \mathbf a(t, \mathbf x) \f$. */
    using AccelerationFunction = std::function<std::vector<double>(double, const std::vector<double> &)>;

    /** @brief Trajectory of a first-order system. */
    struct OdeSolution {
        std::vector<double> t; ///< Times \f$ t_0 \ldots t_{end} \f$, length steps + 1.
        std::vector<double> y; ///< States, row-major (steps + 1) × dimension: @c y[i*dimension + j].
        std::size_t dimension = 0; ///< Number of state variables.
    };

    /** @brief Trajectory of a second-order system. */
    struct SecondOrderSolution {
        std::vector<double> t; ///< Times \f$ t_0 \ldots t_{end} \f$, length steps + 1.
        std::vector<double> position; ///< Positions, row-major (steps + 1) × dimension.
        std::vector<double> velocity; ///< Velocities, row-major (steps + 1) × dimension.
        std::size_t dimension = 0; ///< Number of coordinates.
    };

    /**
     * @brief Forward (explicit) Euler, \f$ \mathbf y_{n+1} = \mathbf y_n + h\, \mathbf f(t_n, \mathbf y_n) \f$.
     *
     * First order: halving \f$ h \f$ halves the global error. Applied to an undamped oscillator it gains
     * energy every step, which makes it a useful baseline for seeing why better methods exist.
     *
     * @param f Right-hand side \f$ \mathbf f(t, \mathbf y) \f$; must return a vector the size of @p y0.
     * @param t0 Initial time; must be finite.
     * @param y0 Initial state; non-empty and finite.
     * @param t_end Final time; must be finite and greater than @p t0.
     * @param steps Number of steps; must be positive.
     * @return Times and states at every step.
     * @throws std::invalid_argument If an argument is invalid or @p f returns the wrong size.
     */
    OdeSolution euler(const OdeFunction &f, double t0, const std::vector<double> &y0, double t_end, int steps);

    /**
     * @brief Classical fourth-order Runge–Kutta.
     *
     * \f$ \mathbf y_{n+1} = \mathbf y_n + \tfrac h6 (\mathbf k_1 + 2\mathbf k_2 + 2\mathbf k_3 + \mathbf k_4) \f$
     * with the usual midpoint stages. Halving \f$ h \f$ reduces the global error about sixteenfold.
     *
     * @param f Right-hand side \f$ \mathbf f(t, \mathbf y) \f$; must return a vector the size of @p y0.
     * @param t0 Initial time; must be finite.
     * @param y0 Initial state; non-empty and finite.
     * @param t_end Final time; must be finite and greater than @p t0.
     * @param steps Number of steps; must be positive.
     * @return Times and states at every step.
     * @throws std::invalid_argument If an argument is invalid or @p f returns the wrong size.
     */
    OdeSolution runge_kutta4(const OdeFunction &f, double t0, const std::vector<double> &y0, double t_end, int steps);

    /**
     * @brief Semi-implicit (symplectic) Euler for \f$ \ddot{\mathbf x} = \mathbf a(t, \mathbf x) \f$.
     *
     * Updates velocity first, then position with the *new* velocity:
     * \f$ \mathbf v_{n+1} = \mathbf v_n + h\, \mathbf a(t_n, \mathbf x_n) \f$,
     * \f$ \mathbf x_{n+1} = \mathbf x_n + h\, \mathbf v_{n+1} \f$.
     * Still first order and no more expensive than euler(), but symplectic: for conservative forces the
     * energy error stays bounded instead of growing, so orbits and oscillators remain stable for long runs.
     *
     * @param acceleration Acceleration \f$ \mathbf a(t, \mathbf x) \f$; must return a vector the size of @p x0.
     * @param t0 Initial time; must be finite.
     * @param x0 Initial position; non-empty and finite.
     * @param v0 Initial velocity; same size as @p x0 and finite.
     * @param t_end Final time; must be finite and greater than @p t0.
     * @param steps Number of steps; must be positive.
     * @return Times, positions and velocities at every step.
     * @throws std::invalid_argument If an argument is invalid or @p acceleration returns the wrong size.
     */
    SecondOrderSolution semi_implicit_euler(const AccelerationFunction &acceleration, double t0,
                                            const std::vector<double> &x0, const std::vector<double> &v0, double t_end,
                                            int steps);

} // namespace physicsforge::numerics

#endif // PHYSICSFORGE_NUMERICS_ODE_H

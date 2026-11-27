/**
 * @file heat_equation.h
 * @brief Solvers for the one-dimensional heat equation with fixed-temperature ends.
 *
 * Solves \f$ \partial_t T = \alpha\, \partial_{xx} T \f$ on \f$ 0 \le x \le L \f$ with
 * \f$ T(0, t) = T_{left} \f$ and \f$ T(L, t) = T_{right} \f$, using three complementary methods:
 *
 * - **Fourier series (analytic).** Writing \f$ T = T_{steady} + u \f$ with a linear steady profile,
 *   \f$ u(x, t) = \sum_n b_n \sin(n\pi x / L)\, e^{-\alpha (n\pi/L)^2 t} \f$. Exact up to truncation.
 * - **FTCS (explicit).** Cheap per step; stable only when \f$ r = \alpha \Delta t / \Delta x^2 \le 1/2 \f$.
 * - **Crank–Nicolson (implicit).** Second order in time and space and unconditionally stable;
 *   each step is an O(N) tridiagonal solve.
 *
 * Grids are uniform and include both end points, \f$ x_j = j \Delta x \f$ for \f$ j = 0 \ldots N-1 \f$.
 *
 * @ingroup thermodynamics
 */

#ifndef PHYSICSFORGE_THERMODYNAMICS_HEAT_EQUATION_H
#define PHYSICSFORGE_THERMODYNAMICS_HEAT_EQUATION_H

#include <span>
#include <vector>

namespace physicsforge::thermodynamics {

    /**
     * @brief Grid Fourier number, \f$ r = \alpha \Delta t / \Delta x^2 \f$.
     * @param diffusivity Thermal diffusivity \f$ \alpha \f$ [m²/s]; must be non-negative.
     * @param time_step Time step \f$ \Delta t \f$ [s]; must be non-negative.
     * @param grid_spacing Grid spacing \f$ \Delta x \f$ [m]; must be positive.
     * @return Dimensionless Fourier number \f$ r \f$.
     * @throws std::invalid_argument If any precondition above is violated.
     */
    double fourier_number(double diffusivity, double time_step, double grid_spacing);

    /**
     * @brief Sine-series coefficients \f$ b_1 \ldots b_N \f$ of the transient part of an initial profile.
     *
     * The linear steady profile between the end temperatures is subtracted first; the remaining integral
     * is evaluated with the trapezoidal rule, which equals the discrete sine transform of the samples.
     *
     * @param initial_temperatures Initial profile [K] on a uniform grid including both ends (at least 3 points).
     * @param n_terms Number of coefficients \f$ N \f$; must be positive.
     * @param left_temperature Temperature at \f$ x = 0 \f$ [K].
     * @param right_temperature Temperature at \f$ x = L \f$ [K].
     * @return Coefficients \f$ b_1 \ldots b_N \f$ [K].
     * @throws std::invalid_argument If the grid is too short, @p n_terms is not positive or a value is non-finite.
     */
    std::vector<double> heat_fourier_coefficients(std::span<const double> initial_temperatures, int n_terms,
                                                  double left_temperature, double right_temperature);

    /**
     * @brief Evaluates the Fourier-series solution \f$ T(x, t) \f$ at the requested positions.
     * @param coefficients Sine coefficients \f$ b_n \f$ from heat_fourier_coefficients() [K].
     * @param length Rod length \f$ L \f$ [m]; must be positive.
     * @param diffusivity Thermal diffusivity \f$ \alpha \f$ [m²/s]; must be non-negative.
     * @param time Time \f$ t \f$ [s]; must be non-negative.
     * @param positions Positions in [0, L] at which to evaluate [m].
     * @param left_temperature Temperature at \f$ x = 0 \f$ [K].
     * @param right_temperature Temperature at \f$ x = L \f$ [K].
     * @return Temperatures at @p positions [K].
     * @throws std::invalid_argument If a position lies outside [0, L] or any precondition is violated.
     */
    std::vector<double> heat_fourier_solution(std::span<const double> coefficients, double length, double diffusivity,
                                              double time, std::span<const double> positions, double left_temperature,
                                              double right_temperature);

    /**
     * @brief Advances a temperature profile with explicit FTCS steps.
     *
     * The first and last samples are replaced by the boundary temperatures.
     *
     * @param initial_temperatures Initial profile [K] on a uniform grid (at least 3 points).
     * @param diffusivity Thermal diffusivity \f$ \alpha \f$ [m²/s]; must be non-negative.
     * @param grid_spacing Grid spacing \f$ \Delta x \f$ [m]; must be positive.
     * @param time_step Time step \f$ \Delta t \f$ [s]; must be non-negative.
     * @param steps Number of steps; must be non-negative.
     * @param left_temperature Temperature at \f$ x = 0 \f$ [K].
     * @param right_temperature Temperature at \f$ x = L \f$ [K].
     * @return Profile after @p steps steps [K].
     * @throws std::invalid_argument If \f$ r > 1/2 \f$ (unstable) or any precondition is violated.
     */
    std::vector<double> heat_equation_explicit(std::span<const double> initial_temperatures, double diffusivity,
                                               double grid_spacing, double time_step, int steps,
                                               double left_temperature, double right_temperature);

    /**
     * @brief Advances a temperature profile with Crank–Nicolson steps (unconditionally stable).
     * @param initial_temperatures Initial profile [K] on a uniform grid (at least 3 points).
     * @param diffusivity Thermal diffusivity \f$ \alpha \f$ [m²/s]; must be non-negative.
     * @param grid_spacing Grid spacing \f$ \Delta x \f$ [m]; must be positive.
     * @param time_step Time step \f$ \Delta t \f$ [s]; must be non-negative.
     * @param steps Number of steps; must be non-negative.
     * @param left_temperature Temperature at \f$ x = 0 \f$ [K].
     * @param right_temperature Temperature at \f$ x = L \f$ [K].
     * @return Profile after @p steps steps [K].
     * @throws std::invalid_argument If any precondition is violated.
     */
    std::vector<double> heat_equation_crank_nicolson(std::span<const double> initial_temperatures, double diffusivity,
                                                     double grid_spacing, double time_step, int steps,
                                                     double left_temperature, double right_temperature);

} // namespace physicsforge::thermodynamics

#endif // PHYSICSFORGE_THERMODYNAMICS_HEAT_EQUATION_H

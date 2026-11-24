// include/physicsforge/thermodynamics/heat_equation.h
//
// Solvers for the one-dimensional heat (diffusion) equation
//
//     dT/dt = alpha d^2T/dx^2,    0 <= x <= L,
//
// with fixed-temperature (Dirichlet) ends T(0, t) = T_left, T(L, t) = T_right.
//
// Three complementary methods are provided:
//
//   * Fourier series (analytic). Fourier introduced his series precisely to solve
//     this problem. Writing T = T_steady + u, with T_steady linear between the end
//     temperatures, gives
//         u(x, t) = sum_n b_n sin(n pi x / L) exp(-alpha (n pi / L)^2 t),
//         b_n     = (2 / L) int_0^L u(x, 0) sin(n pi x / L) dx.
//     Exact for all t up to series truncation; high modes decay fastest.
//   * FTCS (explicit finite differences). Cheap per step, conditionally stable:
//     requires the Fourier number r = alpha dt / dx^2 <= 1/2.
//   * Crank-Nicolson (implicit finite differences). Second order in time and space
//     and unconditionally stable; each step is an O(N) tridiagonal solve.
//
// Grids are uniform and include both end points: x_j = j dx, j = 0 .. N-1.

#ifndef PHYSICSFORGE_THERMODYNAMICS_HEAT_EQUATION_H
#define PHYSICSFORGE_THERMODYNAMICS_HEAT_EQUATION_H

#include <span>
#include <vector>

namespace physicsforge::thermodynamics {

    // Grid Fourier number r = alpha dt / dx^2. Requires alpha >= 0, dt >= 0, dx > 0.
    double fourier_number(double diffusivity, double time_step, double grid_spacing);

    // Sine-series coefficients b_1..b_n of the transient part of an initial profile sampled
    // on a uniform grid that includes both ends (needs >= 3 samples). The steady linear
    // profile between left/right temperatures is subtracted first. Integration uses the
    // trapezoidal rule, which equals the discrete sine transform of the samples.
    std::vector<double> heat_fourier_coefficients(std::span<const double> initial_temperatures, int n_terms,
                                                  double left_temperature, double right_temperature);

    // Evaluates the Fourier-series solution T(x, t) at each position in [0, L].
    std::vector<double> heat_fourier_solution(std::span<const double> coefficients, double length, double diffusivity,
                                              double time, std::span<const double> positions, double left_temperature,
                                              double right_temperature);

    // Advances the temperature profile `steps` explicit FTCS steps. The first and last
    // samples are overwritten with the boundary temperatures. Throws if r > 1/2.
    std::vector<double> heat_equation_explicit(std::span<const double> initial_temperatures, double diffusivity,
                                               double grid_spacing, double time_step, int steps,
                                               double left_temperature, double right_temperature);

    // Advances the temperature profile `steps` Crank-Nicolson steps (unconditionally stable).
    std::vector<double> heat_equation_crank_nicolson(std::span<const double> initial_temperatures, double diffusivity,
                                                     double grid_spacing, double time_step, int steps,
                                                     double left_temperature, double right_temperature);

} // namespace physicsforge::thermodynamics

#endif // PHYSICSFORGE_THERMODYNAMICS_HEAT_EQUATION_H

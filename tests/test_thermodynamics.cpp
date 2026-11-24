// tests/test_thermodynamics.cpp
#include <cmath>
#include <stdexcept>
#include <vector>

#include "physicsforge/constants.h"
#include "physicsforge/thermodynamics/heat.h"
#include "physicsforge/thermodynamics/heat_equation.h"
#include "test_framework.h"

namespace {
    // Uniform grid of n points on [0, length], ends included.
    std::vector<double> grid(const std::size_t n, const double length) {
        std::vector<double> x(n);
        for (std::size_t i = 0; i < n; ++i) {
            x[i] = length * static_cast<double>(i) / static_cast<double>(n - 1);
        }
        return x;
    }
} // namespace

void test_thermodynamics() {
    using namespace physicsforge::thermodynamics;
    namespace c = physicsforge::constants;

    // --- Heat laws ---
    PF_CHECK_CLOSE(sensible_heat(1.0, 4186.0, 10.0), 41860.0);
    PF_CHECK_CLOSE(latent_heat(2.0, 334000.0), 668000.0);
    PF_CHECK_CLOSE(fourier_heat_flux(200.0, -50.0), 10000.0); // heat flows down the gradient
    PF_CHECK_CLOSE(conduction_rate(0.8, 10.0, 20.0, 0.2), 800.0);
    // Resistance form agrees with Fourier's law.
    PF_CHECK_CLOSE(heat_rate_from_resistance(20.0, conduction_resistance(0.2, 0.8, 10.0)), 800.0);
    PF_CHECK_CLOSE(convection_resistance(25.0, 2.0), 0.02);
    PF_CHECK_CLOSE(newton_cooling_rate(10.0, 2.0, 80.0, 20.0), 1200.0);
    PF_CHECK_CLOSE(newton_cooling_temperature(90.0, 20.0, 0.1, 0.0), 90.0);
    PF_CHECK_CLOSE(newton_cooling_temperature(90.0, 20.0, 0.1, 10.0), 20.0 + 70.0 * std::exp(-1.0));
    PF_CHECK_CLOSE(radiated_power(1.0, 1.0, 1000.0), c::STEFAN_BOLTZMANN * 1e12);
    PF_CHECK_CLOSE_TOL(net_radiated_power(0.9, 1.0, 300.0, 300.0), 0.0, 0.0, 1e-12);
    PF_CHECK_THROWS(radiated_power(1.5, 1.0, 300.0), std::invalid_argument);
    PF_CHECK_CLOSE(internal_energy_change(100.0, 40.0), 60.0);
    PF_CHECK_CLOSE(carnot_efficiency(300.0, 600.0), 0.5);
    PF_CHECK_CLOSE(carnot_cop_refrigerator(250.0, 300.0), 5.0);
    PF_CHECK_THROWS(carnot_efficiency(600.0, 300.0), std::invalid_argument);
    // One mole at STP occupies ~22.414 L.
    PF_CHECK_CLOSE_TOL(ideal_gas_volume(1.0, 273.15, c::STANDARD_ATMOSPHERE), 0.022414, 1e-4, 0.0);
    PF_CHECK_CLOSE(ideal_gas_pressure(1.0, 300.0, ideal_gas_volume(1.0, 300.0, 1e5)), 1e5);
    PF_CHECK_CLOSE(celsius_to_kelvin(0.0), 273.15);
    PF_CHECK_CLOSE(kelvin_to_celsius(373.15), 100.0);
    PF_CHECK_THROWS(celsius_to_kelvin(-300.0), std::invalid_argument);
    PF_CHECK_CLOSE(thermal_diffusivity(237.0, 2700.0, 897.0), 237.0 / (2700.0 * 897.0));

    // --- Heat equation ---
    // Exact single-mode solution: T(x, 0) = sin(pi x / L) decays as exp(-alpha (pi/L)^2 t).
    constexpr double length = 1.0;
    constexpr double alpha = 1e-2;
    constexpr double t_end = 5.0;
    constexpr std::size_t n = 51;
    const std::vector<double> x = grid(n, length);
    std::vector<double> initial(n);
    for (std::size_t i = 0; i < n; ++i) {
        initial[i] = std::sin(c::PI * x[i] / length);
    }
    const double decay = std::exp(-alpha * c::PI * c::PI * t_end / (length * length));

    // Fourier series reproduces the mode exactly.
    const auto coeffs = heat_fourier_coefficients(initial, 5, 0.0, 0.0);
    PF_CHECK_CLOSE_TOL(coeffs[0], 1.0, 1e-12, 0.0);
    PF_CHECK_CLOSE_TOL(coeffs[1], 0.0, 0.0, 1e-12);
    const auto fourier = heat_fourier_solution(coeffs, length, alpha, t_end, x, 0.0, 0.0);
    PF_CHECK_CLOSE_TOL(fourier[n / 2], decay, 1e-12, 0.0);

    // Finite-difference schemes converge to the same answer.
    const double dx = length / static_cast<double>(n - 1);
    const double dt = 0.01;
    const int steps = static_cast<int>(t_end / dt);
    const auto ftcs = heat_equation_explicit(initial, alpha, dx, dt, steps, 0.0, 0.0);
    const auto cn = heat_equation_crank_nicolson(initial, alpha, dx, dt, steps, 0.0, 0.0);
    PF_CHECK_CLOSE_TOL(ftcs[n / 2], decay, 1e-3, 0.0);
    PF_CHECK_CLOSE_TOL(cn[n / 2], decay, 1e-3, 0.0);
    PF_CHECK(ftcs.front() == 0.0 && ftcs.back() == 0.0);

    // FTCS rejects unstable steps; Crank-Nicolson accepts them and stays bounded.
    PF_CHECK_THROWS(heat_equation_explicit(initial, alpha, dx, 1.0, 10, 0.0, 0.0), std::invalid_argument);
    const auto cn_big = heat_equation_crank_nicolson(initial, alpha, dx, 1.0, 5, 0.0, 0.0);
    PF_CHECK_CLOSE_TOL(cn_big[n / 2], decay, 2e-3, 0.0);

    // Non-zero boundaries relax to the linear steady state.
    const std::vector<double> cold(n, 0.0);
    const auto steady =
            heat_fourier_solution(heat_fourier_coefficients(cold, 50, 100.0, 0.0), length, alpha, 1e6, x, 100.0, 0.0);
    PF_CHECK_CLOSE_TOL(steady[n / 2], 50.0, 1e-9, 0.0);
    const auto steady_cn = heat_equation_crank_nicolson(cold, alpha, dx, 10.0, 2000, 100.0, 0.0);
    PF_CHECK_CLOSE_TOL(steady_cn[n / 2], 50.0, 1e-6, 0.0);

    PF_CHECK_CLOSE(fourier_number(1e-4, 0.5, 0.01), 0.5);
}

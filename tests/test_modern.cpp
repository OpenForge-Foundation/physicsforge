// tests/test_modern.cpp
#include <cmath>
#include <stdexcept>

#include "physicsforge/constants.h"
#include "physicsforge/modern/quantum.h"
#include "physicsforge/modern/relativity.h"
#include "test_framework.h"

void test_modern() {
    using namespace physicsforge::modern;
    namespace c = physicsforge::constants;

    // --- Quantum ---
    PF_CHECK_CLOSE(photon_energy_from_frequency(1.0), c::PLANCK);
    const double green = 550e-9;
    const double green_frequency = c::SPEED_OF_LIGHT / green;
    PF_CHECK_CLOSE(wavelength_from_frequency(green_frequency), green);
    PF_CHECK_CLOSE(photon_energy_from_wavelength(green), photon_energy_from_frequency(green_frequency));
    PF_CHECK_CLOSE_TOL(photon_energy_from_wavelength(green) / c::ELECTRON_VOLT, 2.254, 1e-3, 0.0);
    PF_CHECK_CLOSE(photon_momentum(green) * c::SPEED_OF_LIGHT, photon_energy_from_wavelength(green));
    PF_CHECK_CLOSE(de_broglie_wavelength(c::ELECTRON_MASS, 1e6), c::PLANCK / (c::ELECTRON_MASS * 1e6));
    PF_CHECK(photoelectric_kinetic_energy(1e14, 2.0 * c::ELECTRON_VOLT) == 0.0); // below threshold
    PF_CHECK_CLOSE(photoelectric_kinetic_energy(1e15, 2.0 * c::ELECTRON_VOLT),
                   c::PLANCK * 1e15 - 2.0 * c::ELECTRON_VOLT);
    // The Sun (~5772 K) peaks near 502 nm.
    PF_CHECK_CLOSE_TOL(wien_peak_wavelength(5772.0), 502.0e-9, 1e-3, 0.0);
    // Planck's law peaks at Wien's wavelength.
    const double peak = wien_peak_wavelength(5772.0);
    PF_CHECK(planck_spectral_radiance(peak, 5772.0) > planck_spectral_radiance(0.98 * peak, 5772.0));
    PF_CHECK(planck_spectral_radiance(peak, 5772.0) > planck_spectral_radiance(1.02 * peak, 5772.0));
    PF_CHECK(planck_spectral_radiance(1e-9, 300.0) == 0.0); // deep Wien tail underflows cleanly
    // Integrating B over wavelength and solid angle (x pi) gives sigma T^4.
    {
        constexpr double T = 1000.0;
        double total = 0.0;
        const double lo = 1e-7, hi = 1e-3;
        constexpr int panels = 200000;
        const double log_step = std::log(hi / lo) / panels;
        for (int i = 0; i < panels; ++i) {
            const double l = lo * std::exp((i + 0.5) * log_step);
            total += planck_spectral_radiance(l, T) * l * log_step;
        }
        PF_CHECK_CLOSE_TOL(c::PI * total, c::STEFAN_BOLTZMANN * T * T * T * T, 1e-4, 0.0);
    }
    PF_CHECK_CLOSE_TOL(bohr_energy_level(1.0) / c::ELECTRON_VOLT, -13.6057, 1e-5, 0.0);
    PF_CHECK_CLOSE(bohr_energy_level(2.0), bohr_energy_level(1.0) / 4.0);
    PF_CHECK_THROWS(bohr_energy_level(1.5), std::invalid_argument);
    PF_CHECK_CLOSE(heisenberg_min_momentum_uncertainty(1e-10), c::REDUCED_PLANCK / 2e-10);

    // --- Relativity ---
    PF_CHECK_CLOSE(lorentz_factor(0.0), 1.0);
    PF_CHECK_CLOSE(lorentz_factor(0.6 * c::SPEED_OF_LIGHT), 1.25);
    PF_CHECK_CLOSE(time_dilation(1.0, 0.6 * c::SPEED_OF_LIGHT), 1.25);
    PF_CHECK_CLOSE(length_contraction(1.0, 0.6 * c::SPEED_OF_LIGHT), 0.8);
    PF_CHECK_THROWS(lorentz_factor(c::SPEED_OF_LIGHT), std::invalid_argument);
    PF_CHECK_CLOSE(rest_energy(1.0), c::SPEED_OF_LIGHT * c::SPEED_OF_LIGHT);
    PF_CHECK_CLOSE(relativistic_kinetic_energy(1.0, 0.6 * c::SPEED_OF_LIGHT), 0.25 * rest_energy(1.0));
    // Low-speed limit matches m v^2 / 2 without cancellation error.
    PF_CHECK_CLOSE_TOL(relativistic_kinetic_energy(1.0, 1.0), 0.5, 1e-12, 0.0);
    PF_CHECK_CLOSE(relativistic_momentum(1.0, 0.6 * c::SPEED_OF_LIGHT), 1.25 * 0.6 * c::SPEED_OF_LIGHT);
}

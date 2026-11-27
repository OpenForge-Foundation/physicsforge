// src/quantum.cpp
#include "physicsforge/quantum.h"

#include <cmath>
#include <limits>

#include "physicsforge/constants.h"
#include "physicsforge/validation.h"

namespace physicsforge::quantum {

    using constants::BOLTZMANN;
    using constants::PLANCK;
    using constants::SPEED_OF_LIGHT;
    using detail::require_finite;
    using detail::require_non_negative;
    using detail::require_non_zero;
    using detail::require_positive;

    double photon_energy_from_frequency(const double frequency) {
        require_non_negative(frequency, "frequency");
        return PLANCK * frequency;
    }

    double photon_energy_from_wavelength(const double wavelength) {
        require_positive(wavelength, "wavelength");
        return PLANCK * SPEED_OF_LIGHT / wavelength;
    }

    double photon_momentum(const double wavelength) {
        require_positive(wavelength, "wavelength");
        return PLANCK / wavelength;
    }

    double wavelength_from_frequency(const double frequency) {
        require_positive(frequency, "frequency");
        return SPEED_OF_LIGHT / frequency;
    }

    double de_broglie_wavelength(const double mass, const double velocity) {
        require_positive(mass, "mass");
        require_non_zero(velocity, "velocity");
        return PLANCK / (mass * std::abs(velocity));
    }

    double photoelectric_kinetic_energy(const double frequency, const double work_function) {
        require_non_negative(frequency, "frequency");
        require_non_negative(work_function, "work_function");
        const double energy = PLANCK * frequency - work_function;
        return energy > 0.0 ? energy : 0.0;
    }

    double planck_spectral_radiance(const double wavelength, const double temperature) {
        require_positive(wavelength, "wavelength");
        require_non_negative(temperature, "temperature");
        if (temperature == 0.0) {
            return 0.0;
        }
        const double exponent = PLANCK * SPEED_OF_LIGHT / (wavelength * BOLTZMANN * temperature);
        const double l2 = wavelength * wavelength;
        const double prefactor = 2.0 * PLANCK * SPEED_OF_LIGHT * SPEED_OF_LIGHT / (l2 * l2 * wavelength);
        // expm1 keeps full precision in the Rayleigh-Jeans regime (small exponent);
        // for huge exponents expm1 overflows to +inf and the radiance correctly becomes 0.
        return prefactor / std::expm1(exponent);
    }

    double wien_peak_wavelength(const double temperature) {
        require_positive(temperature, "temperature");
        return constants::WIEN_DISPLACEMENT / temperature;
    }

    double bohr_energy_level(const double principal_quantum_number) {
        detail::require_in_range(principal_quantum_number, 1.0, std::numeric_limits<double>::max(),
                                 "principal_quantum_number");
        if (principal_quantum_number != std::floor(principal_quantum_number)) {
            detail::throw_invalid("principal_quantum_number", "an integer", principal_quantum_number);
        }
        constexpr double e2 = constants::ELEMENTARY_CHARGE * constants::ELEMENTARY_CHARGE;
        constexpr double rydberg_energy =
                constants::ELECTRON_MASS * e2 * e2 /
                (8.0 * constants::VACUUM_PERMITTIVITY * constants::VACUUM_PERMITTIVITY * PLANCK * PLANCK);
        return -rydberg_energy / (principal_quantum_number * principal_quantum_number);
    }

    double heisenberg_min_momentum_uncertainty(const double position_uncertainty) {
        require_positive(position_uncertainty, "position_uncertainty");
        return constants::REDUCED_PLANCK / (2.0 * position_uncertainty);
    }

} // namespace physicsforge::quantum

// src/electromagnetism/maxwell.cpp
#include "physicsforge/electromagnetism/maxwell.h"

#include <cmath>
#include <limits>

#include "physicsforge/constants.h"
#include "physicsforge/detail/validation.h"

namespace physicsforge::electromagnetism {

    using constants::SPEED_OF_LIGHT;
    using constants::VACUUM_PERMEABILITY;
    using constants::VACUUM_PERMITTIVITY;
    using detail::require_finite;
    using detail::require_non_negative;
    using detail::require_positive;

    double gauss_electric_flux(const double enclosed_charge) {
        require_finite(enclosed_charge, "enclosed_charge");
        return enclosed_charge / VACUUM_PERMITTIVITY;
    }

    double gauss_enclosed_charge(const double electric_flux) {
        require_finite(electric_flux, "electric_flux");
        return VACUUM_PERMITTIVITY * electric_flux;
    }

    double magnetic_flux(const double field_strength, const double area, const double angle) {
        require_non_negative(field_strength, "field_strength");
        require_non_negative(area, "area");
        require_finite(angle, "angle");
        return field_strength * area * std::cos(angle);
    }

    double faraday_emf(const double flux_rate_of_change, const double turns) {
        require_finite(flux_rate_of_change, "flux_rate_of_change");
        detail::require_in_range(turns, 1.0, std::numeric_limits<double>::max(), "turns");
        return -turns * flux_rate_of_change;
    }

    double displacement_current(const double electric_flux_rate_of_change) {
        require_finite(electric_flux_rate_of_change, "electric_flux_rate_of_change");
        return VACUUM_PERMITTIVITY * electric_flux_rate_of_change;
    }

    double ampere_maxwell_circulation(const double enclosed_current, const double electric_flux_rate_of_change) {
        require_finite(enclosed_current, "enclosed_current");
        return VACUUM_PERMEABILITY * (enclosed_current + displacement_current(electric_flux_rate_of_change));
    }

    double magnetic_field_straight_wire(const double current, const double distance) {
        require_finite(current, "current");
        require_positive(distance, "distance");
        return VACUUM_PERMEABILITY * current / (constants::TAU * distance);
    }

    double magnetic_field_solenoid(const double turns_per_length, const double current) {
        require_non_negative(turns_per_length, "turns_per_length");
        require_finite(current, "current");
        return VACUUM_PERMEABILITY * turns_per_length * current;
    }

    double magnetic_field_loop_center(const double current, const double radius) {
        require_finite(current, "current");
        require_positive(radius, "radius");
        return VACUUM_PERMEABILITY * current / (2.0 * radius);
    }

    math::Vec3 lorentz_force(const double charge, const math::Vec3 &electric_field, const math::Vec3 &velocity,
                             const math::Vec3 &magnetic_field) {
        require_finite(charge, "charge");
        return (electric_field + math::cross(velocity, magnetic_field)) * charge;
    }

    double em_wave_speed(const double relative_permittivity, const double relative_permeability) {
        require_positive(relative_permittivity, "relative_permittivity");
        require_positive(relative_permeability, "relative_permeability");
        return SPEED_OF_LIGHT / std::sqrt(relative_permittivity * relative_permeability);
    }

    double speed_of_light_from_vacuum_constants() { return 1.0 / std::sqrt(VACUUM_PERMEABILITY * VACUUM_PERMITTIVITY); }

    double em_wave_intensity(const double electric_field_amplitude) {
        require_non_negative(electric_field_amplitude, "electric_field_amplitude");
        return 0.5 * SPEED_OF_LIGHT * VACUUM_PERMITTIVITY * electric_field_amplitude * electric_field_amplitude;
    }

} // namespace physicsforge::electromagnetism

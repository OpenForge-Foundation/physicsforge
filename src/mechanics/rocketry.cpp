// src/mechanics/rocketry.cpp
#include "physicsforge/mechanics/rocketry.h"

#include <cmath>

#include "physicsforge/detail/validation.h"

namespace physicsforge::mechanics {

    using detail::require_finite;
    using detail::require_non_negative;
    using detail::require_positive;

    double rocket_delta_v(const double exhaust_velocity, const double initial_mass, const double final_mass) {
        require_non_negative(exhaust_velocity, "exhaust_velocity");
        require_positive(final_mass, "final_mass");
        require_positive(initial_mass, "initial_mass");
        if (initial_mass < final_mass) {
            detail::throw_invalid("initial_mass must be greater than or equal to final_mass");
        }
        return exhaust_velocity * std::log(initial_mass / final_mass);
    }

    double rocket_mass_ratio(const double delta_v, const double exhaust_velocity) {
        require_non_negative(delta_v, "delta_v");
        require_positive(exhaust_velocity, "exhaust_velocity");
        return std::exp(delta_v / exhaust_velocity);
    }

    double rocket_propellant_mass(const double delta_v, const double exhaust_velocity, const double final_mass) {
        require_non_negative(delta_v, "delta_v");
        require_positive(exhaust_velocity, "exhaust_velocity");
        require_positive(final_mass, "final_mass");
        return final_mass * std::expm1(delta_v / exhaust_velocity);
    }

    double exhaust_velocity_from_isp(const double specific_impulse, const double standard_gravity) {
        require_non_negative(specific_impulse, "specific_impulse");
        require_positive(standard_gravity, "standard_gravity");
        return specific_impulse * standard_gravity;
    }

    double specific_impulse(const double exhaust_velocity, const double standard_gravity) {
        require_non_negative(exhaust_velocity, "exhaust_velocity");
        require_positive(standard_gravity, "standard_gravity");
        return exhaust_velocity / standard_gravity;
    }

    double rocket_thrust(const double mass_flow_rate, const double exhaust_velocity, const double exit_pressure,
                         const double ambient_pressure, const double exit_area) {
        require_non_negative(mass_flow_rate, "mass_flow_rate");
        require_finite(exhaust_velocity, "exhaust_velocity");
        require_non_negative(exit_pressure, "exit_pressure");
        require_non_negative(ambient_pressure, "ambient_pressure");
        require_non_negative(exit_area, "exit_area");
        return mass_flow_rate * exhaust_velocity + (exit_pressure - ambient_pressure) * exit_area;
    }

    double burn_time(const double propellant_mass, const double mass_flow_rate) {
        require_non_negative(propellant_mass, "propellant_mass");
        require_positive(mass_flow_rate, "mass_flow_rate");
        return propellant_mass / mass_flow_rate;
    }

} // namespace physicsforge::mechanics

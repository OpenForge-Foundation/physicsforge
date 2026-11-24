// src/mechanics/kinematics.cpp
#include "physicsforge/mechanics/kinematics.h"

#include <cmath>

#include "physicsforge/detail/validation.h"

namespace physicsforge::mechanics {

    using detail::require_finite;
    using detail::require_non_negative;
    using detail::require_positive;

    double displacement(const double initial_velocity, const double acceleration, const double time) {
        require_finite(initial_velocity, "initial_velocity");
        require_finite(acceleration, "acceleration");
        require_non_negative(time, "time");
        return time * (initial_velocity + 0.5 * acceleration * time);
    }

    double final_velocity(const double initial_velocity, const double acceleration, const double time) {
        require_finite(initial_velocity, "initial_velocity");
        require_finite(acceleration, "acceleration");
        require_non_negative(time, "time");
        return initial_velocity + acceleration * time;
    }

    double final_velocity_from_displacement(const double initial_velocity, const double acceleration,
                                            const double displacement) {
        require_finite(initial_velocity, "initial_velocity");
        require_finite(acceleration, "acceleration");
        require_finite(displacement, "displacement");
        const double v_squared = initial_velocity * initial_velocity + 2.0 * acceleration * displacement;
        if (v_squared < 0.0) {
            detail::throw_invalid("v0^2 + 2*a*s must be non-negative; the displacement is never reached");
        }
        return std::sqrt(v_squared);
    }

    double average_velocity(const double displacement, const double time) {
        require_finite(displacement, "displacement");
        require_positive(time, "time");
        return displacement / time;
    }

    double average_acceleration(const double velocity_change, const double time) {
        require_finite(velocity_change, "velocity_change");
        require_positive(time, "time");
        return velocity_change / time;
    }

    double projectile_range(const double initial_speed, const double launch_angle, const double gravity) {
        require_non_negative(initial_speed, "initial_speed");
        require_finite(launch_angle, "launch_angle");
        require_positive(gravity, "gravity");
        return initial_speed * initial_speed * std::sin(2.0 * launch_angle) / gravity;
    }

    double projectile_time_of_flight(const double initial_speed, const double launch_angle, const double gravity) {
        require_non_negative(initial_speed, "initial_speed");
        require_finite(launch_angle, "launch_angle");
        require_positive(gravity, "gravity");
        return 2.0 * initial_speed * std::sin(launch_angle) / gravity;
    }

    double projectile_max_height(const double initial_speed, const double launch_angle, const double gravity) {
        require_non_negative(initial_speed, "initial_speed");
        require_finite(launch_angle, "launch_angle");
        require_positive(gravity, "gravity");
        const double vertical = initial_speed * std::sin(launch_angle);
        return vertical * vertical / (2.0 * gravity);
    }

} // namespace physicsforge::mechanics

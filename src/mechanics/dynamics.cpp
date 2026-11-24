// src/mechanics/dynamics.cpp
#include "physicsforge/mechanics/dynamics.h"

#include <cmath>

#include "physicsforge/detail/validation.h"

namespace physicsforge::mechanics {

    using detail::require_finite;
    using detail::require_non_negative;
    using detail::require_positive;

    double net_force(const double mass, const double acceleration) {
        require_non_negative(mass, "mass");
        require_finite(acceleration, "acceleration");
        return mass * acceleration;
    }

    double acceleration_from_force(const double force, const double mass) {
        require_finite(force, "force");
        require_positive(mass, "mass");
        return force / mass;
    }

    double weight(const double mass, const double gravity) {
        require_non_negative(mass, "mass");
        require_non_negative(gravity, "gravity");
        return mass * gravity;
    }

    double momentum(const double mass, const double velocity) {
        require_non_negative(mass, "mass");
        require_finite(velocity, "velocity");
        return mass * velocity;
    }

    double impulse(const double force, const double duration) {
        require_finite(force, "force");
        require_non_negative(duration, "duration");
        return force * duration;
    }

    double kinetic_energy(const double mass, const double velocity) {
        require_non_negative(mass, "mass");
        require_finite(velocity, "velocity");
        return 0.5 * mass * velocity * velocity;
    }

    double potential_energy(const double mass, const double height, const double gravity) {
        require_non_negative(mass, "mass");
        require_finite(height, "height");
        require_non_negative(gravity, "gravity");
        return mass * gravity * height;
    }

    double work(const double force, const double distance, const double angle) {
        require_finite(force, "force");
        require_non_negative(distance, "distance");
        require_finite(angle, "angle");
        return force * distance * std::cos(angle);
    }

    double power(const double work, const double time) {
        require_finite(work, "work");
        require_positive(time, "time");
        return work / time;
    }

    double spring_force(const double spring_constant, const double extension) {
        require_non_negative(spring_constant, "spring_constant");
        require_finite(extension, "extension");
        return -spring_constant * extension;
    }

    double spring_potential_energy(const double spring_constant, const double extension) {
        require_non_negative(spring_constant, "spring_constant");
        require_finite(extension, "extension");
        return 0.5 * spring_constant * extension * extension;
    }

    double friction_force(const double coefficient, const double normal_force) {
        require_non_negative(coefficient, "coefficient");
        require_non_negative(normal_force, "normal_force");
        return coefficient * normal_force;
    }

} // namespace physicsforge::mechanics

// src/mechanics/rotation.cpp
#include "physicsforge/mechanics/rotation.h"

#include <cmath>

#include "physicsforge/detail/validation.h"

namespace physicsforge::mechanics {

    using detail::require_finite;
    using detail::require_non_negative;
    using detail::require_positive;

    namespace {
        // Shared implementation for I = k m s^2 shapes.
        double inertia(const double factor, const double mass, const double size, const char *size_name) {
            require_non_negative(mass, "mass");
            require_non_negative(size, size_name);
            return factor * mass * size * size;
        }
    } // namespace

    double torque(const double lever_arm, const double force, const double angle) {
        require_non_negative(lever_arm, "lever_arm");
        require_finite(force, "force");
        require_finite(angle, "angle");
        return lever_arm * force * std::sin(angle);
    }

    math::Vec3 torque_vector(const math::Vec3 &lever_arm, const math::Vec3 &force) {
        return math::cross(lever_arm, force);
    }

    double angular_acceleration(const double torque, const double moment_of_inertia) {
        require_finite(torque, "torque");
        require_positive(moment_of_inertia, "moment_of_inertia");
        return torque / moment_of_inertia;
    }

    double angular_momentum(const double moment_of_inertia, const double angular_velocity) {
        require_non_negative(moment_of_inertia, "moment_of_inertia");
        require_finite(angular_velocity, "angular_velocity");
        return moment_of_inertia * angular_velocity;
    }

    math::Vec3 angular_momentum_vector(const math::Vec3 &position, const math::Vec3 &linear_momentum) {
        return math::cross(position, linear_momentum);
    }

    double rotational_kinetic_energy(const double moment_of_inertia, const double angular_velocity) {
        require_non_negative(moment_of_inertia, "moment_of_inertia");
        require_finite(angular_velocity, "angular_velocity");
        return 0.5 * moment_of_inertia * angular_velocity * angular_velocity;
    }

    double tangential_velocity(const double angular_velocity, const double radius) {
        require_finite(angular_velocity, "angular_velocity");
        require_non_negative(radius, "radius");
        return angular_velocity * radius;
    }

    double centripetal_acceleration(const double speed, const double radius) {
        require_finite(speed, "speed");
        require_positive(radius, "radius");
        return speed * speed / radius;
    }

    double centripetal_force(const double mass, const double speed, const double radius) {
        require_non_negative(mass, "mass");
        return mass * centripetal_acceleration(speed, radius);
    }

    double moment_of_inertia_point_mass(const double mass, const double radius) {
        return inertia(1.0, mass, radius, "radius");
    }
    double moment_of_inertia_solid_sphere(const double mass, const double radius) {
        return inertia(2.0 / 5.0, mass, radius, "radius");
    }
    double moment_of_inertia_hollow_sphere(const double mass, const double radius) {
        return inertia(2.0 / 3.0, mass, radius, "radius");
    }
    double moment_of_inertia_solid_cylinder(const double mass, const double radius) {
        return inertia(0.5, mass, radius, "radius");
    }
    double moment_of_inertia_hoop(const double mass, const double radius) {
        return inertia(1.0, mass, radius, "radius");
    }
    double moment_of_inertia_rod_center(const double mass, const double length) {
        return inertia(1.0 / 12.0, mass, length, "length");
    }
    double moment_of_inertia_rod_end(const double mass, const double length) {
        return inertia(1.0 / 3.0, mass, length, "length");
    }

    double parallel_axis(const double moment_of_inertia_cm, const double mass, const double distance) {
        require_non_negative(moment_of_inertia_cm, "moment_of_inertia_cm");
        require_non_negative(mass, "mass");
        require_finite(distance, "distance");
        return moment_of_inertia_cm + mass * distance * distance;
    }

} // namespace physicsforge::mechanics

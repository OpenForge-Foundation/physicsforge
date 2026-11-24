// src/mechanics/gravitation.cpp
#include "physicsforge/mechanics/gravitation.h"

#include <cmath>

#include "physicsforge/constants.h"
#include "physicsforge/detail/validation.h"

namespace physicsforge::mechanics {

    using constants::GRAVITATIONAL_CONSTANT;
    using detail::require_finite;
    using detail::require_non_negative;
    using detail::require_positive;

    double gravitational_force(const double mass1, const double mass2, const double distance) {
        require_non_negative(mass1, "mass1");
        require_non_negative(mass2, "mass2");
        require_positive(distance, "distance");
        return GRAVITATIONAL_CONSTANT * mass1 * mass2 / (distance * distance);
    }

    math::Vec3 gravitational_force_vector(const double mass1, const double mass2, const math::Vec3 &position1,
                                          const math::Vec3 &position2) {
        require_non_negative(mass1, "mass1");
        require_non_negative(mass2, "mass2");
        const math::Vec3 separation = position1 - position2;
        const double r = math::norm(separation);
        require_positive(r, "distance between bodies");
        return separation * (-GRAVITATIONAL_CONSTANT * mass1 * mass2 / (r * r * r));
    }

    double gravitational_field_strength(const double mass, const double distance) {
        require_non_negative(mass, "mass");
        require_positive(distance, "distance");
        return GRAVITATIONAL_CONSTANT * mass / (distance * distance);
    }

    double gravitational_potential_energy(const double mass1, const double mass2, const double distance) {
        require_non_negative(mass1, "mass1");
        require_non_negative(mass2, "mass2");
        require_positive(distance, "distance");
        return -GRAVITATIONAL_CONSTANT * mass1 * mass2 / distance;
    }

    double escape_velocity(const double central_mass, const double distance) {
        require_non_negative(central_mass, "central_mass");
        require_positive(distance, "distance");
        return std::sqrt(2.0 * GRAVITATIONAL_CONSTANT * central_mass / distance);
    }

    double orbital_velocity(const double central_mass, const double distance) {
        require_non_negative(central_mass, "central_mass");
        require_positive(distance, "distance");
        return std::sqrt(GRAVITATIONAL_CONSTANT * central_mass / distance);
    }

    double orbital_period(const double semi_major_axis, const double central_mass) {
        require_positive(semi_major_axis, "semi_major_axis");
        require_positive(central_mass, "central_mass");
        return constants::TAU *
               std::sqrt(semi_major_axis * semi_major_axis * semi_major_axis / (GRAVITATIONAL_CONSTANT * central_mass));
    }

} // namespace physicsforge::mechanics

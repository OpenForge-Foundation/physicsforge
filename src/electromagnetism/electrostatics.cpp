// src/electromagnetism/electrostatics.cpp
#include "physicsforge/electromagnetism/electrostatics.h"

#include <limits>

#include "physicsforge/constants.h"
#include "physicsforge/detail/validation.h"

namespace physicsforge::electromagnetism {

    using constants::COULOMB_CONSTANT;
    using detail::require_finite;
    using detail::require_non_negative;
    using detail::require_positive;

    double coulomb_force(const double charge1, const double charge2, const double distance) {
        require_finite(charge1, "charge1");
        require_finite(charge2, "charge2");
        require_positive(distance, "distance");
        return COULOMB_CONSTANT * charge1 * charge2 / (distance * distance);
    }

    math::Vec3 coulomb_force_vector(const double charge1, const double charge2, const math::Vec3 &position1,
                                    const math::Vec3 &position2) {
        require_finite(charge1, "charge1");
        require_finite(charge2, "charge2");
        const math::Vec3 separation = position1 - position2;
        const double r = math::norm(separation);
        require_positive(r, "distance between charges");
        return separation * (COULOMB_CONSTANT * charge1 * charge2 / (r * r * r));
    }

    double electric_field_point_charge(const double charge, const double distance) {
        require_finite(charge, "charge");
        require_positive(distance, "distance");
        return COULOMB_CONSTANT * charge / (distance * distance);
    }

    math::Vec3 electric_field_vector(const double charge, const math::Vec3 &source, const math::Vec3 &point) {
        require_finite(charge, "charge");
        const math::Vec3 separation = point - source;
        const double r = math::norm(separation);
        require_positive(r, "distance from charge");
        return separation * (COULOMB_CONSTANT * charge / (r * r * r));
    }

    double electric_potential_point_charge(const double charge, const double distance) {
        require_finite(charge, "charge");
        require_positive(distance, "distance");
        return COULOMB_CONSTANT * charge / distance;
    }

    double electric_potential_energy(const double charge1, const double charge2, const double distance) {
        require_finite(charge1, "charge1");
        require_finite(charge2, "charge2");
        require_positive(distance, "distance");
        return COULOMB_CONSTANT * charge1 * charge2 / distance;
    }

    double parallel_plate_capacitance(const double relative_permittivity, const double area, const double separation) {
        detail::require_in_range(relative_permittivity, 1.0, std::numeric_limits<double>::max(),
                                 "relative_permittivity");
        require_non_negative(area, "area");
        require_positive(separation, "separation");
        return relative_permittivity * constants::VACUUM_PERMITTIVITY * area / separation;
    }

    double capacitor_energy(const double capacitance, const double voltage) {
        require_non_negative(capacitance, "capacitance");
        require_finite(voltage, "voltage");
        return 0.5 * capacitance * voltage * voltage;
    }

} // namespace physicsforge::electromagnetism

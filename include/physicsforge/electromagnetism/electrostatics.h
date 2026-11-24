// include/physicsforge/electromagnetism/electrostatics.h
//
// Coulomb's law, point-charge fields and potentials, and capacitors.
// SI units: coulombs (C), metres (m), newtons (N), volts (V), farads (F).

#ifndef PHYSICSFORGE_ELECTROMAGNETISM_ELECTROSTATICS_H
#define PHYSICSFORGE_ELECTROMAGNETISM_ELECTROSTATICS_H

#include "physicsforge/math/vector3.h"

namespace physicsforge::electromagnetism {

    // Coulomb's law (signed): F = k_e q1 q2 / r^2. Positive = repulsive. Requires r > 0.
    double coulomb_force(double charge1, double charge2, double distance);

    // Force vector on charge 1 due to charge 2: F_12 = k_e q1 q2 (r1 - r2) / |r1 - r2|^3.
    math::Vec3 coulomb_force_vector(double charge1, double charge2, const math::Vec3 &position1,
                                    const math::Vec3 &position2);

    // Radial electric field of a point charge (signed): E = k_e q / r^2. Requires r > 0.
    double electric_field_point_charge(double charge, double distance);

    // Electric field vector at `point` from a point charge at `source`.
    math::Vec3 electric_field_vector(double charge, const math::Vec3 &source, const math::Vec3 &point);

    // Electric potential of a point charge (zero at infinity): V = k_e q / r. Requires r > 0.
    double electric_potential_point_charge(double charge, double distance);

    // Electric potential energy of two point charges: U = k_e q1 q2 / r. Requires r > 0.
    double electric_potential_energy(double charge1, double charge2, double distance);

    // Ideal parallel-plate capacitance: C = eps_r eps_0 A / d. Requires eps_r >= 1, A >= 0, d > 0.
    double parallel_plate_capacitance(double relative_permittivity, double area, double separation);

    // Energy stored in a capacitor: U = C V^2 / 2. Requires C >= 0.
    double capacitor_energy(double capacitance, double voltage);

} // namespace physicsforge::electromagnetism

#endif // PHYSICSFORGE_ELECTROMAGNETISM_ELECTROSTATICS_H

// include/physicsforge/mechanics/gravitation.h
//
// Newton's law of universal gravitation and two-body orbital mechanics.
//
// Orbital helpers assume a test body whose mass is negligible compared with the
// central mass M (restricted two-body problem). SI units throughout.

#ifndef PHYSICSFORGE_MECHANICS_GRAVITATION_H
#define PHYSICSFORGE_MECHANICS_GRAVITATION_H

#include "physicsforge/math/vector3.h"

namespace physicsforge::mechanics {

    // Magnitude of the mutual attraction: F = G m1 m2 / r^2. Requires m1, m2 >= 0, r > 0.
    double gravitational_force(double mass1, double mass2, double distance);

    // Force vector on body 1 due to body 2: F_12 = -G m1 m2 (r1 - r2) / |r1 - r2|^3.
    // Throws if the bodies coincide.
    math::Vec3 gravitational_force_vector(double mass1, double mass2, const math::Vec3 &position1,
                                          const math::Vec3 &position2);

    // Gravitational field strength of a point/spherical mass: g = G M / r^2. Requires M >= 0, r > 0.
    double gravitational_field_strength(double mass, double distance);

    // Two-body gravitational potential energy: U = -G m1 m2 / r. Requires m1, m2 >= 0, r > 0.
    double gravitational_potential_energy(double mass1, double mass2, double distance);

    // Escape velocity from radius r: v_esc = sqrt(2 G M / r). Requires M >= 0, r > 0.
    double escape_velocity(double central_mass, double distance);

    // Circular orbital speed at radius r: v = sqrt(G M / r). Requires M >= 0, r > 0.
    double orbital_velocity(double central_mass, double distance);

    // Orbital period from Kepler's third law: T = 2 pi sqrt(a^3 / (G M)). Requires a > 0, M > 0.
    double orbital_period(double semi_major_axis, double central_mass);

} // namespace physicsforge::mechanics

#endif // PHYSICSFORGE_MECHANICS_GRAVITATION_H

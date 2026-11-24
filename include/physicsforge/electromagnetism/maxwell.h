// include/physicsforge/electromagnetism/maxwell.h
//
// Maxwell's equations (integral form) and their standard closed-form consequences.
//
//   Gauss (electric):   oint E . dA = Q_enc / eps_0
//   Gauss (magnetic):   oint B . dA = 0
//   Faraday:            oint E . dl = -dPhi_B/dt
//   Ampere-Maxwell:     oint B . dl = mu_0 (I_enc + eps_0 dPhi_E/dt)
//
// The differential forms (div E = rho/eps_0, curl B = mu_0 J + mu_0 eps_0 dE/dt, ...)
// can be evaluated numerically for arbitrary fields with
// physicsforge::numerics::divergence and physicsforge::numerics::curl.

#ifndef PHYSICSFORGE_ELECTROMAGNETISM_MAXWELL_H
#define PHYSICSFORGE_ELECTROMAGNETISM_MAXWELL_H

#include "physicsforge/math/vector3.h"

namespace physicsforge::electromagnetism {

    // Gauss's law: electric flux through a closed surface, Phi_E = Q_enc / eps_0.
    double gauss_electric_flux(double enclosed_charge);

    // Gauss's law inverted: enclosed charge from flux, Q_enc = eps_0 Phi_E.
    double gauss_enclosed_charge(double electric_flux);

    // Magnetic flux through a flat surface in a uniform field: Phi_B = B A cos(theta),
    // theta between B and the surface normal. Requires B >= 0, A >= 0.
    double magnetic_flux(double field_strength, double area, double angle);

    // Faraday's law (with Lenz sign): EMF = -N dPhi_B/dt. Requires N >= 1.
    double faraday_emf(double flux_rate_of_change, double turns);

    // Maxwell's displacement current: I_d = eps_0 dPhi_E/dt.
    double displacement_current(double electric_flux_rate_of_change);

    // Ampere-Maxwell law: circulation of B, oint B . dl = mu_0 (I_enc + eps_0 dPhi_E/dt).
    double ampere_maxwell_circulation(double enclosed_current, double electric_flux_rate_of_change);

    // Field of an infinite straight wire (Ampere's law): B = mu_0 I / (2 pi r). Requires r > 0.
    double magnetic_field_straight_wire(double current, double distance);

    // Field inside an ideal long solenoid: B = mu_0 n I, n = turns per metre. Requires n >= 0.
    double magnetic_field_solenoid(double turns_per_length, double current);

    // Field at the centre of a circular current loop: B = mu_0 I / (2 R). Requires R > 0.
    double magnetic_field_loop_center(double current, double radius);

    // Lorentz force on a point charge: F = q (E + v x B).
    math::Vec3 lorentz_force(double charge, const math::Vec3 &electric_field, const math::Vec3 &velocity,
                             const math::Vec3 &magnetic_field);

    // Phase velocity of an EM wave in a linear medium: v = c / sqrt(eps_r mu_r). Requires eps_r, mu_r > 0.
    double em_wave_speed(double relative_permittivity, double relative_permeability);

    // Speed of light derived from the vacuum constants: c = 1 / sqrt(mu_0 eps_0).
    double speed_of_light_from_vacuum_constants();

    // Time-averaged intensity of a plane wave in vacuum: I = c eps_0 E0^2 / 2. Requires E0 >= 0.
    double em_wave_intensity(double electric_field_amplitude);

} // namespace physicsforge::electromagnetism

#endif // PHYSICSFORGE_ELECTROMAGNETISM_MAXWELL_H

/**
 * @file maxwell.h
 * @brief Maxwell's equations in integral form and their standard closed-form consequences.
 *
 * | Law             | Integral form                                                                     |
 * |-----------------|-----------------------------------------------------------------------------------|
 * | Gauss           | \f$ \oint \mathbf E \cdot d\mathbf A = Q_{enc} / \varepsilon_0 \f$                 |
 * | Gauss (magnetic)| \f$ \oint \mathbf B \cdot d\mathbf A = 0 \f$                                       |
 * | Faraday         | \f$ \oint \mathbf E \cdot d\boldsymbol\ell = -d\Phi_B/dt \f$                       |
 * | Ampère–Maxwell  | \f$ \oint \mathbf B \cdot d\boldsymbol\ell = \mu_0 (I_{enc} + \varepsilon_0\, d\Phi_E/dt) \f$ |
 *
 * The differential forms can be evaluated numerically for arbitrary fields with
 * physicsforge::numerics::divergence and physicsforge::numerics::curl.
 *
 * @ingroup electromagnetism
 */

#ifndef PHYSICSFORGE_ELECTROMAGNETISM_MAXWELL_H
#define PHYSICSFORGE_ELECTROMAGNETISM_MAXWELL_H

#include "physicsforge/math/vector3.h"

namespace physicsforge::electromagnetism {

    /**
     * @brief Gauss's law: electric flux through a closed surface, \f$ \Phi_E = Q_{enc} / \varepsilon_0 \f$.
     * @param enclosed_charge Charge enclosed by the surface \f$ Q_{enc} \f$ [C].
     * @return Electric flux \f$ \Phi_E \f$ [V·m].
     */
    double gauss_electric_flux(double enclosed_charge);

    /**
     * @brief Gauss's law inverted: enclosed charge from flux, \f$ Q_{enc} = \varepsilon_0 \Phi_E \f$.
     * @param electric_flux Electric flux through a closed surface \f$ \Phi_E \f$ [V·m].
     * @return Enclosed charge [C].
     */
    double gauss_enclosed_charge(double electric_flux);

    /**
     * @brief Magnetic flux through a flat surface in a uniform field, \f$ \Phi_B = B A \cos\theta \f$.
     * @param field_strength Field magnitude \f$ B \f$ [T]; must be non-negative.
     * @param area Surface area \f$ A \f$ [m²]; must be non-negative.
     * @param angle Angle between the field and the surface normal \f$ \theta \f$ [rad].
     * @return Magnetic flux \f$ \Phi_B \f$ [Wb].
     * @throws std::invalid_argument If @p field_strength or @p area is negative.
     */
    double magnetic_flux(double field_strength, double area, double angle);

    /**
     * @brief Faraday's law with Lenz's sign, \f$ \mathcal E = -N\, d\Phi_B/dt \f$.
     * @param flux_rate_of_change Rate of change of flux through one turn \f$ d\Phi_B/dt \f$ [Wb/s].
     * @param turns Number of turns \f$ N \f$; must be at least 1.
     * @return Induced EMF [V].
     * @throws std::invalid_argument If @p turns is less than 1.
     */
    double faraday_emf(double flux_rate_of_change, double turns);

    /**
     * @brief Maxwell's displacement current, \f$ I_d = \varepsilon_0\, d\Phi_E/dt \f$.
     * @param electric_flux_rate_of_change Rate of change of electric flux \f$ d\Phi_E/dt \f$ [V·m/s].
     * @return Displacement current [A].
     */
    double displacement_current(double electric_flux_rate_of_change);

    /**
     * @brief Ampère–Maxwell law: circulation of \f$ \mathbf B \f$ around a closed loop,
     *        \f$ \oint \mathbf B \cdot d\boldsymbol\ell = \mu_0 (I_{enc} + \varepsilon_0\, d\Phi_E/dt) \f$.
     * @param enclosed_current Conduction current through the loop \f$ I_{enc} \f$ [A].
     * @param electric_flux_rate_of_change Rate of change of electric flux through the loop [V·m/s].
     * @return Circulation [T·m].
     */
    double ampere_maxwell_circulation(double enclosed_current, double electric_flux_rate_of_change);

    /**
     * @brief Field of an infinite straight wire (Ampère's law), \f$ B = \mu_0 I / (2\pi r) \f$.
     * @param current Wire current \f$ I \f$ [A].
     * @param distance Perpendicular distance from the wire \f$ r \f$ [m]; must be positive.
     * @return Field magnitude [T].
     * @throws std::invalid_argument If @p distance is not positive.
     */
    double magnetic_field_straight_wire(double current, double distance);

    /**
     * @brief Field inside an ideal long solenoid, \f$ B = \mu_0 n I \f$.
     * @param turns_per_length Winding density \f$ n \f$ [turns/m]; must be non-negative.
     * @param current Current \f$ I \f$ [A].
     * @return Field magnitude [T].
     * @throws std::invalid_argument If @p turns_per_length is negative.
     */
    double magnetic_field_solenoid(double turns_per_length, double current);

    /**
     * @brief Field at the centre of a circular current loop, \f$ B = \mu_0 I / (2R) \f$.
     * @param current Loop current \f$ I \f$ [A].
     * @param radius Loop radius \f$ R \f$ [m]; must be positive.
     * @return Field magnitude [T].
     * @throws std::invalid_argument If @p radius is not positive.
     */
    double magnetic_field_loop_center(double current, double radius);

    /**
     * @brief Lorentz force on a point charge, \f$ \mathbf F = q (\mathbf E + \mathbf v \times \mathbf B) \f$.
     * @param charge Charge \f$ q \f$ [C].
     * @param electric_field Electric field \f$ \mathbf E \f$ [V/m].
     * @param velocity Particle velocity \f$ \mathbf v \f$ [m/s].
     * @param magnetic_field Magnetic field \f$ \mathbf B \f$ [T].
     * @return Force vector [N].
     */
    math::Vec3 lorentz_force(double charge, const math::Vec3 &electric_field, const math::Vec3 &velocity,
                             const math::Vec3 &magnetic_field);

    /**
     * @brief Phase velocity of an EM wave in a linear medium, \f$ v = c / \sqrt{\varepsilon_r \mu_r} \f$.
     * @param relative_permittivity Relative permittivity \f$ \varepsilon_r \f$; must be positive.
     * @param relative_permeability Relative permeability \f$ \mu_r \f$; must be positive.
     * @return Wave speed [m/s].
     * @throws std::invalid_argument If either input is not positive.
     */
    double em_wave_speed(double relative_permittivity, double relative_permeability);

    /**
     * @brief Speed of light derived from the vacuum constants, \f$ c = 1 / \sqrt{\mu_0 \varepsilon_0} \f$.
     * @return Speed of light [m/s].
     */
    double speed_of_light_from_vacuum_constants();

    /**
     * @brief Time-averaged intensity of a plane wave in vacuum, \f$ I = \tfrac12 c \varepsilon_0 E_0^2 \f$.
     * @param electric_field_amplitude Peak electric field \f$ E_0 \f$ [V/m]; must be non-negative.
     * @return Intensity [W/m²].
     * @throws std::invalid_argument If @p electric_field_amplitude is negative.
     */
    double em_wave_intensity(double electric_field_amplitude);

} // namespace physicsforge::electromagnetism

#endif // PHYSICSFORGE_ELECTROMAGNETISM_MAXWELL_H

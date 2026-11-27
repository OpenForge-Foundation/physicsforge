/**
 * @file gravitation.h
 * @brief Newton's law of universal gravitation and two-body orbital mechanics.
 *
 * Orbital helpers assume a test body whose mass is negligible compared with the central mass
 * (restricted two-body problem). SI units throughout.
 *
 * @ingroup mechanics
 */

#ifndef PHYSICSFORGE_MECHANICS_GRAVITATION_H
#define PHYSICSFORGE_MECHANICS_GRAVITATION_H

#include "physicsforge/math/vector3.h"

namespace physicsforge::mechanics {

    /**
     * @brief Magnitude of the mutual gravitational attraction, \f$ F = G m_1 m_2 / r^2 \f$.
     * @param mass1 Mass of the first body \f$ m_1 \f$ [kg]; must be non-negative.
     * @param mass2 Mass of the second body \f$ m_2 \f$ [kg]; must be non-negative.
     * @param distance Centre-to-centre separation \f$ r \f$ [m]; must be positive.
     * @return Force magnitude \f$ F \f$ [N].
     * @throws std::invalid_argument If a mass is negative or @p distance is not positive.
     */
    double gravitational_force(double mass1, double mass2, double distance);

    /**
     * @brief Gravitational force on body 1 due to body 2,
     *        \f$ \mathbf F_{12} = -G m_1 m_2 (\mathbf r_1 - \mathbf r_2) / |\mathbf r_1 - \mathbf r_2|^3 \f$.
     * @param mass1 Mass of body 1 [kg]; must be non-negative.
     * @param mass2 Mass of body 2 [kg]; must be non-negative.
     * @param position1 Position of body 1 [m].
     * @param position2 Position of body 2 [m].
     * @return Force vector acting on body 1 [N].
     * @throws std::invalid_argument If a mass is negative or the positions coincide.
     */
    math::Vec3 gravitational_force_vector(double mass1, double mass2, const math::Vec3 &position1,
                                          const math::Vec3 &position2);

    /**
     * @brief Field strength of a point or spherical mass, \f$ g = G M / r^2 \f$.
     * @param mass Source mass \f$ M \f$ [kg]; must be non-negative.
     * @param distance Distance from the centre \f$ r \f$ [m]; must be positive.
     * @return Field strength \f$ g \f$ [m/s²].
     * @throws std::invalid_argument If @p mass is negative or @p distance is not positive.
     */
    double gravitational_field_strength(double mass, double distance);

    /**
     * @brief Two-body gravitational potential energy, \f$ U = -G m_1 m_2 / r \f$.
     * @param mass1 Mass of the first body [kg]; must be non-negative.
     * @param mass2 Mass of the second body [kg]; must be non-negative.
     * @param distance Separation \f$ r \f$ [m]; must be positive.
     * @return Potential energy \f$ U \f$ [J] (zero at infinite separation).
     * @throws std::invalid_argument If a mass is negative or @p distance is not positive.
     */
    double gravitational_potential_energy(double mass1, double mass2, double distance);

    /**
     * @brief Escape velocity, \f$ v_{esc} = \sqrt{2 G M / r} \f$.
     * @param central_mass Central mass \f$ M \f$ [kg]; must be non-negative.
     * @param distance Distance from the centre \f$ r \f$ [m]; must be positive.
     * @return Escape speed [m/s].
     * @throws std::invalid_argument If @p central_mass is negative or @p distance is not positive.
     */
    double escape_velocity(double central_mass, double distance);

    /**
     * @brief Circular orbital speed, \f$ v = \sqrt{G M / r} \f$.
     * @param central_mass Central mass \f$ M \f$ [kg]; must be non-negative.
     * @param distance Orbital radius \f$ r \f$ [m]; must be positive.
     * @return Orbital speed [m/s].
     * @throws std::invalid_argument If @p central_mass is negative or @p distance is not positive.
     */
    double orbital_velocity(double central_mass, double distance);

    /**
     * @brief Orbital period from Kepler's third law, \f$ T = 2\pi \sqrt{a^3 / (G M)} \f$.
     * @param semi_major_axis Semi-major axis \f$ a \f$ [m]; must be positive.
     * @param central_mass Central mass \f$ M \f$ [kg]; must be positive.
     * @return Period \f$ T \f$ [s].
     * @throws std::invalid_argument If either input is not positive.
     */
    double orbital_period(double semi_major_axis, double central_mass);

} // namespace physicsforge::mechanics

#endif // PHYSICSFORGE_MECHANICS_GRAVITATION_H

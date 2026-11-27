/**
 * @file electrostatics.h
 * @brief Coulomb's law, point-charge fields and potentials, and capacitors.
 *
 * SI units: C, m, N, V, F. Signed results follow the usual convention (positive force = repulsive).
 *
 * @ingroup electromagnetism
 */

#ifndef PHYSICSFORGE_ELECTROMAGNETISM_ELECTROSTATICS_H
#define PHYSICSFORGE_ELECTROMAGNETISM_ELECTROSTATICS_H

#include "physicsforge/math/vector3.h"

namespace physicsforge::electromagnetism {

    /**
     * @brief Coulomb's law (signed), \f$ F = k_e q_1 q_2 / r^2 \f$.
     * @param charge1 First charge \f$ q_1 \f$ [C].
     * @param charge2 Second charge \f$ q_2 \f$ [C].
     * @param distance Separation \f$ r \f$ [m]; must be positive.
     * @return Force [N]; positive is repulsive, negative attractive.
     * @throws std::invalid_argument If @p distance is not positive.
     */
    double coulomb_force(double charge1, double charge2, double distance);

    /**
     * @brief Coulomb force on charge 1 due to charge 2,
     *        \f$ \mathbf F_{12} = k_e q_1 q_2 (\mathbf r_1 - \mathbf r_2) / |\mathbf r_1 - \mathbf r_2|^3 \f$.
     * @param charge1 Charge 1 [C].
     * @param charge2 Charge 2 [C].
     * @param position1 Position of charge 1 [m].
     * @param position2 Position of charge 2 [m].
     * @return Force vector on charge 1 [N].
     * @throws std::invalid_argument If the positions coincide.
     */
    math::Vec3 coulomb_force_vector(double charge1, double charge2, const math::Vec3 &position1,
                                    const math::Vec3 &position2);

    /**
     * @brief Radial electric field of a point charge (signed), \f$ E = k_e q / r^2 \f$.
     * @param charge Source charge \f$ q \f$ [C].
     * @param distance Distance from the charge \f$ r \f$ [m]; must be positive.
     * @return Radial field component [V/m]; positive points away from the charge.
     * @throws std::invalid_argument If @p distance is not positive.
     */
    double electric_field_point_charge(double charge, double distance);

    /**
     * @brief Electric field vector of a point charge, \f$ \mathbf E = k_e q \hat{\mathbf r} / r^2 \f$.
     * @param charge Source charge \f$ q \f$ [C].
     * @param source Position of the charge [m].
     * @param point Field point [m].
     * @return Field at @p point [V/m].
     * @throws std::invalid_argument If @p point coincides with @p source.
     */
    math::Vec3 electric_field_vector(double charge, const math::Vec3 &source, const math::Vec3 &point);

    /**
     * @brief Electric potential of a point charge, \f$ V = k_e q / r \f$ (zero at infinity).
     * @param charge Source charge \f$ q \f$ [C].
     * @param distance Distance from the charge \f$ r \f$ [m]; must be positive.
     * @return Potential [V].
     * @throws std::invalid_argument If @p distance is not positive.
     */
    double electric_potential_point_charge(double charge, double distance);

    /**
     * @brief Electric potential energy of two point charges, \f$ U = k_e q_1 q_2 / r \f$.
     * @param charge1 First charge [C].
     * @param charge2 Second charge [C].
     * @param distance Separation \f$ r \f$ [m]; must be positive.
     * @return Potential energy [J].
     * @throws std::invalid_argument If @p distance is not positive.
     */
    double electric_potential_energy(double charge1, double charge2, double distance);

    /**
     * @brief Ideal parallel-plate capacitance, \f$ C = \varepsilon_r \varepsilon_0 A / d \f$.
     * @param relative_permittivity Dielectric constant \f$ \varepsilon_r \f$; must be at least 1.
     * @param area Plate area \f$ A \f$ [m²]; must be non-negative.
     * @param separation Plate separation \f$ d \f$ [m]; must be positive.
     * @return Capacitance [F].
     * @throws std::invalid_argument If any precondition above is violated.
     */
    double parallel_plate_capacitance(double relative_permittivity, double area, double separation);

    /**
     * @brief Energy stored in a capacitor, \f$ U = \tfrac12 C V^2 \f$.
     * @param capacitance Capacitance \f$ C \f$ [F]; must be non-negative.
     * @param voltage Voltage across the capacitor \f$ V \f$ [V].
     * @return Stored energy [J].
     * @throws std::invalid_argument If @p capacitance is negative.
     */
    double capacitor_energy(double capacitance, double voltage);

} // namespace physicsforge::electromagnetism

#endif // PHYSICSFORGE_ELECTROMAGNETISM_ELECTROSTATICS_H

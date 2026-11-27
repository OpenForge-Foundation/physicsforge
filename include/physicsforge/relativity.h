/**
 * @file relativity.h
 * @brief Special relativity: Lorentz factor, time dilation, length contraction, energy and momentum.
 *
 * Every speed must satisfy \f$ |v| < c \f$.
 *
 * @ingroup relativity
 */

#ifndef PHYSICSFORGE_RELATIVITY_H
#define PHYSICSFORGE_RELATIVITY_H

/** @defgroup relativity Relativity
 *  @brief Special relativity.
 */

namespace physicsforge::relativity {

    /**
     * @brief Lorentz factor, \f$ \gamma = 1 / \sqrt{1 - v^2/c^2} \f$.
     * @param velocity Speed \f$ v \f$ [m/s]; must satisfy \f$ |v| < c \f$.
     * @return Lorentz factor \f$ \gamma \ge 1 \f$.
     * @throws std::invalid_argument If \f$ |v| \ge c \f$.
     */
    double lorentz_factor(double velocity);

    /**
     * @brief Dilated (coordinate) time interval, \f$ \Delta t = \gamma \Delta t_0 \f$.
     * @param proper_time Proper time interval \f$ \Delta t_0 \f$ [s]; must be non-negative.
     * @param velocity Relative speed \f$ v \f$ [m/s]; must satisfy \f$ |v| < c \f$.
     * @return Dilated interval [s].
     * @throws std::invalid_argument If @p proper_time is negative or \f$ |v| \ge c \f$.
     */
    double time_dilation(double proper_time, double velocity);

    /**
     * @brief Contracted length, \f$ L = L_0 / \gamma \f$.
     * @param proper_length Proper (rest-frame) length \f$ L_0 \f$ [m]; must be non-negative.
     * @param velocity Relative speed \f$ v \f$ [m/s]; must satisfy \f$ |v| < c \f$.
     * @return Contracted length [m].
     * @throws std::invalid_argument If @p proper_length is negative or \f$ |v| \ge c \f$.
     */
    double length_contraction(double proper_length, double velocity);

    /**
     * @brief Mass–energy equivalence, \f$ E = m c^2 \f$.
     * @param mass Rest mass \f$ m \f$ [kg]; must be non-negative.
     * @return Rest energy [J].
     * @throws std::invalid_argument If @p mass is negative.
     */
    double rest_energy(double mass);

    /**
     * @brief Relativistic momentum, \f$ p = \gamma m v \f$.
     * @param mass Rest mass \f$ m \f$ [kg]; must be non-negative.
     * @param velocity Velocity \f$ v \f$ [m/s]; must satisfy \f$ |v| < c \f$.
     * @return Momentum [kg·m/s].
     * @throws std::invalid_argument If @p mass is negative or \f$ |v| \ge c \f$.
     */
    double relativistic_momentum(double mass, double velocity);

    /**
     * @brief Relativistic kinetic energy, \f$ K = (\gamma - 1) m c^2 \f$.
     *
     * \f$ \gamma - 1 \f$ is computed without cancellation, so the Newtonian limit \f$ \tfrac12 m v^2 \f$
     * is reproduced to full precision at low speed.
     *
     * @param mass Rest mass \f$ m \f$ [kg]; must be non-negative.
     * @param velocity Velocity \f$ v \f$ [m/s]; must satisfy \f$ |v| < c \f$.
     * @return Kinetic energy [J].
     * @throws std::invalid_argument If @p mass is negative or \f$ |v| \ge c \f$.
     */
    double relativistic_kinetic_energy(double mass, double velocity);

} // namespace physicsforge::relativity

#endif // PHYSICSFORGE_RELATIVITY_H

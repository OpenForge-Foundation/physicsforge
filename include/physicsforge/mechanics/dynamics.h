/**
 * @file dynamics.h
 * @brief Newtonian dynamics: forces, momentum, work and energy.
 *
 * SI units throughout (kg, m, s, N, J, W); angles in radians.
 *
 * @ingroup mechanics
 */

#ifndef PHYSICSFORGE_MECHANICS_DYNAMICS_H
#define PHYSICSFORGE_MECHANICS_DYNAMICS_H

namespace physicsforge::mechanics {

    /**
     * @brief Newton's second law, \f$ F = m a \f$.
     * @param mass Mass \f$ m \f$ [kg]; must be non-negative.
     * @param acceleration Acceleration \f$ a \f$ [m/s²].
     * @return Net force \f$ F \f$ [N].
     * @throws std::invalid_argument If @p mass is negative or any input is non-finite.
     */
    double net_force(double mass, double acceleration);

    /**
     * @brief Acceleration from Newton's second law, \f$ a = F / m \f$.
     * @param force Net force \f$ F \f$ [N].
     * @param mass Mass \f$ m \f$ [kg]; must be positive.
     * @return Acceleration \f$ a \f$ [m/s²].
     * @throws std::invalid_argument If @p mass is not positive.
     */
    double acceleration_from_force(double force, double mass);

    /**
     * @brief Weight near a surface, \f$ W = m g \f$.
     * @param mass Mass \f$ m \f$ [kg]; must be non-negative.
     * @param gravity Gravitational acceleration \f$ g \f$ [m/s²]; must be non-negative.
     * @return Weight \f$ W \f$ [N].
     * @throws std::invalid_argument If @p mass or @p gravity is negative.
     */
    double weight(double mass, double gravity);

    /**
     * @brief Linear momentum, \f$ p = m v \f$.
     * @param mass Mass \f$ m \f$ [kg]; must be non-negative.
     * @param velocity Velocity \f$ v \f$ [m/s].
     * @return Momentum \f$ p \f$ [kg·m/s].
     * @throws std::invalid_argument If @p mass is negative.
     */
    double momentum(double mass, double velocity);

    /**
     * @brief Impulse of a constant force, \f$ J = F \Delta t \f$.
     * @param force Force \f$ F \f$ [N].
     * @param duration Duration \f$ \Delta t \f$ [s]; must be non-negative.
     * @return Impulse \f$ J \f$ [N·s].
     * @throws std::invalid_argument If @p duration is negative.
     */
    double impulse(double force, double duration);

    /**
     * @brief Translational kinetic energy, \f$ K = \tfrac12 m v^2 \f$.
     * @param mass Mass \f$ m \f$ [kg]; must be non-negative.
     * @param velocity Velocity \f$ v \f$ [m/s].
     * @return Kinetic energy \f$ K \f$ [J].
     * @throws std::invalid_argument If @p mass is negative.
     */
    double kinetic_energy(double mass, double velocity);

    /**
     * @brief Gravitational potential energy in a uniform field, \f$ U = m g h \f$.
     * @param mass Mass \f$ m \f$ [kg]; must be non-negative.
     * @param height Height above the reference level \f$ h \f$ [m].
     * @param gravity Gravitational acceleration \f$ g \f$ [m/s²]; must be non-negative.
     * @return Potential energy \f$ U \f$ [J].
     * @throws std::invalid_argument If @p mass or @p gravity is negative.
     */
    double potential_energy(double mass, double height, double gravity);

    /**
     * @brief Work done by a constant force, \f$ W = F d \cos\theta \f$.
     * @param force Force magnitude \f$ F \f$ [N].
     * @param distance Distance moved \f$ d \f$ [m]; must be non-negative.
     * @param angle Angle between force and displacement \f$ \theta \f$ [rad].
     * @return Work \f$ W \f$ [J].
     * @throws std::invalid_argument If @p distance is negative.
     */
    double work(double force, double distance, double angle);

    /**
     * @brief Average power, \f$ P = W / t \f$.
     * @param work Work done \f$ W \f$ [J].
     * @param time Time interval \f$ t \f$ [s]; must be positive.
     * @return Power \f$ P \f$ [W].
     * @throws std::invalid_argument If @p time is not positive.
     */
    double power(double work, double time);

    /**
     * @brief Hooke's law restoring force, \f$ F = -k x \f$.
     * @param spring_constant Spring constant \f$ k \f$ [N/m]; must be non-negative.
     * @param extension Extension from equilibrium \f$ x \f$ [m].
     * @return Restoring force \f$ F \f$ [N].
     * @throws std::invalid_argument If @p spring_constant is negative.
     */
    double spring_force(double spring_constant, double extension);

    /**
     * @brief Elastic potential energy of an ideal spring, \f$ U = \tfrac12 k x^2 \f$.
     * @param spring_constant Spring constant \f$ k \f$ [N/m]; must be non-negative.
     * @param extension Extension from equilibrium \f$ x \f$ [m].
     * @return Stored energy \f$ U \f$ [J].
     * @throws std::invalid_argument If @p spring_constant is negative.
     */
    double spring_potential_energy(double spring_constant, double extension);

    /**
     * @brief Kinetic (Coulomb) friction magnitude, \f$ f = \mu N \f$.
     * @param coefficient Coefficient of friction \f$ \mu \f$; must be non-negative.
     * @param normal_force Normal force \f$ N \f$ [N]; must be non-negative.
     * @return Friction force \f$ f \f$ [N].
     * @throws std::invalid_argument If either input is negative.
     */
    double friction_force(double coefficient, double normal_force);

} // namespace physicsforge::mechanics

#endif // PHYSICSFORGE_MECHANICS_DYNAMICS_H

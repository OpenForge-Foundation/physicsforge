/**
 * @file kinematics.h
 * @brief One-dimensional constant-acceleration kinematics and ideal projectile motion.
 *
 * SI units throughout (m, s, m/s, m/s², rad). Projectile functions assume level ground and no drag.
 *
 * @ingroup mechanics
 */

#ifndef PHYSICSFORGE_MECHANICS_KINEMATICS_H
#define PHYSICSFORGE_MECHANICS_KINEMATICS_H

/** @defgroup mechanics Mechanics
 *  @brief Kinematics, Newtonian dynamics, gravitation, rotation and rocketry.
 */

namespace physicsforge::mechanics {

    /**
     * @brief Displacement under constant acceleration, \f$ s = v_0 t + \tfrac12 a t^2 \f$.
     * @param initial_velocity Initial velocity \f$ v_0 \f$ [m/s].
     * @param acceleration Constant acceleration \f$ a \f$ [m/s²].
     * @param time Elapsed time \f$ t \f$ [s]; must be non-negative.
     * @return Displacement \f$ s \f$ [m].
     * @throws std::invalid_argument If @p time is negative or any input is non-finite.
     */
    double displacement(double initial_velocity, double acceleration, double time);

    /**
     * @brief Velocity under constant acceleration, \f$ v = v_0 + a t \f$.
     * @param initial_velocity Initial velocity \f$ v_0 \f$ [m/s].
     * @param acceleration Constant acceleration \f$ a \f$ [m/s²].
     * @param time Elapsed time \f$ t \f$ [s]; must be non-negative.
     * @return Final velocity \f$ v \f$ [m/s].
     * @throws std::invalid_argument If @p time is negative or any input is non-finite.
     */
    double final_velocity(double initial_velocity, double acceleration, double time);

    /**
     * @brief Speed after a displacement under constant acceleration, \f$ |v| = \sqrt{v_0^2 + 2 a s} \f$.
     * @param initial_velocity Initial velocity \f$ v_0 \f$ [m/s].
     * @param acceleration Constant acceleration \f$ a \f$ [m/s²].
     * @param displacement Displacement \f$ s \f$ [m].
     * @return Final speed \f$ |v| \f$ [m/s].
     * @throws std::invalid_argument If \f$ v_0^2 + 2 a s < 0 \f$ (the displacement is never reached).
     */
    double final_velocity_from_displacement(double initial_velocity, double acceleration, double displacement);

    /**
     * @brief Average velocity, \f$ \bar v = \Delta x / \Delta t \f$.
     * @param displacement Displacement \f$ \Delta x \f$ [m].
     * @param time Time interval \f$ \Delta t \f$ [s]; must be positive.
     * @return Average velocity [m/s].
     * @throws std::invalid_argument If @p time is not positive.
     */
    double average_velocity(double displacement, double time);

    /**
     * @brief Average acceleration, \f$ \bar a = \Delta v / \Delta t \f$.
     * @param velocity_change Change in velocity \f$ \Delta v \f$ [m/s].
     * @param time Time interval \f$ \Delta t \f$ [s]; must be positive.
     * @return Average acceleration [m/s²].
     * @throws std::invalid_argument If @p time is not positive.
     */
    double average_acceleration(double velocity_change, double time);

    /**
     * @brief Horizontal range of a projectile, \f$ R = v_0^2 \sin(2\theta) / g \f$.
     * @param initial_speed Launch speed \f$ v_0 \f$ [m/s]; must be non-negative.
     * @param launch_angle Elevation angle \f$ \theta \f$ [rad].
     * @param gravity Gravitational acceleration \f$ g \f$ [m/s²]; must be positive.
     * @return Range \f$ R \f$ [m].
     * @throws std::invalid_argument If @p initial_speed is negative or @p gravity is not positive.
     */
    double projectile_range(double initial_speed, double launch_angle, double gravity);

    /**
     * @brief Time of flight of a projectile, \f$ T = 2 v_0 \sin\theta / g \f$.
     * @param initial_speed Launch speed \f$ v_0 \f$ [m/s]; must be non-negative.
     * @param launch_angle Elevation angle \f$ \theta \f$ [rad].
     * @param gravity Gravitational acceleration \f$ g \f$ [m/s²]; must be positive.
     * @return Time of flight \f$ T \f$ [s].
     * @throws std::invalid_argument If @p initial_speed is negative or @p gravity is not positive.
     */
    double projectile_time_of_flight(double initial_speed, double launch_angle, double gravity);

    /**
     * @brief Apex height above the launch point, \f$ H = (v_0 \sin\theta)^2 / (2 g) \f$.
     * @param initial_speed Launch speed \f$ v_0 \f$ [m/s]; must be non-negative.
     * @param launch_angle Elevation angle \f$ \theta \f$ [rad].
     * @param gravity Gravitational acceleration \f$ g \f$ [m/s²]; must be positive.
     * @return Maximum height \f$ H \f$ [m].
     * @throws std::invalid_argument If @p initial_speed is negative or @p gravity is not positive.
     */
    double projectile_max_height(double initial_speed, double launch_angle, double gravity);

} // namespace physicsforge::mechanics

#endif // PHYSICSFORGE_MECHANICS_KINEMATICS_H

/**
 * @file rotation.h
 * @brief Rigid-body rotation: torque, angular momentum, moments of inertia and circular motion.
 *
 * SI units; angles in radians; angular velocity in rad/s.
 *
 * @ingroup mechanics
 */

#ifndef PHYSICSFORGE_MECHANICS_ROTATION_H
#define PHYSICSFORGE_MECHANICS_ROTATION_H

#include "physicsforge/math/vector3.h"

namespace physicsforge::mechanics {

    /**
     * @brief Torque magnitude, \f$ \tau = r F \sin\theta \f$.
     * @param lever_arm Distance from the axis to the point of application \f$ r \f$ [m]; must be non-negative.
     * @param force Force magnitude \f$ F \f$ [N].
     * @param angle Angle between the lever arm and the force \f$ \theta \f$ [rad].
     * @return Torque \f$ \tau \f$ [N·m].
     * @throws std::invalid_argument If @p lever_arm is negative.
     */
    double torque(double lever_arm, double force, double angle);

    /**
     * @brief Torque vector, \f$ \boldsymbol\tau = \mathbf r \times \mathbf F \f$.
     * @param lever_arm Position of the point of application relative to the axis [m].
     * @param force Applied force [N].
     * @return Torque vector [N·m].
     */
    math::Vec3 torque_vector(const math::Vec3 &lever_arm, const math::Vec3 &force);

    /**
     * @brief Newton's second law for rotation, \f$ \alpha = \tau / I \f$.
     * @param torque Net torque \f$ \tau \f$ [N·m].
     * @param moment_of_inertia Moment of inertia \f$ I \f$ [kg·m²]; must be positive.
     * @return Angular acceleration \f$ \alpha \f$ [rad/s²].
     * @throws std::invalid_argument If @p moment_of_inertia is not positive.
     */
    double angular_acceleration(double torque, double moment_of_inertia);

    /**
     * @brief Angular momentum about a fixed axis, \f$ L = I \omega \f$.
     * @param moment_of_inertia Moment of inertia \f$ I \f$ [kg·m²]; must be non-negative.
     * @param angular_velocity Angular velocity \f$ \omega \f$ [rad/s].
     * @return Angular momentum \f$ L \f$ [kg·m²/s].
     * @throws std::invalid_argument If @p moment_of_inertia is negative.
     */
    double angular_momentum(double moment_of_inertia, double angular_velocity);

    /**
     * @brief Angular momentum of a particle, \f$ \mathbf L = \mathbf r \times \mathbf p \f$.
     * @param position Position relative to the origin [m].
     * @param linear_momentum Linear momentum [kg·m/s].
     * @return Angular momentum vector [kg·m²/s].
     */
    math::Vec3 angular_momentum_vector(const math::Vec3 &position, const math::Vec3 &linear_momentum);

    /**
     * @brief Rotational kinetic energy, \f$ K = \tfrac12 I \omega^2 \f$.
     * @param moment_of_inertia Moment of inertia \f$ I \f$ [kg·m²]; must be non-negative.
     * @param angular_velocity Angular velocity \f$ \omega \f$ [rad/s].
     * @return Kinetic energy \f$ K \f$ [J].
     * @throws std::invalid_argument If @p moment_of_inertia is negative.
     */
    double rotational_kinetic_energy(double moment_of_inertia, double angular_velocity);

    /**
     * @brief Tangential speed of a point on a rotating body, \f$ v = \omega r \f$.
     * @param angular_velocity Angular velocity \f$ \omega \f$ [rad/s].
     * @param radius Distance from the axis \f$ r \f$ [m]; must be non-negative.
     * @return Tangential velocity \f$ v \f$ [m/s].
     * @throws std::invalid_argument If @p radius is negative.
     */
    double tangential_velocity(double angular_velocity, double radius);

    /**
     * @brief Centripetal acceleration, \f$ a_c = v^2 / r \f$.
     * @param speed Speed \f$ v \f$ [m/s].
     * @param radius Radius of the circular path \f$ r \f$ [m]; must be positive.
     * @return Centripetal acceleration [m/s²].
     * @throws std::invalid_argument If @p radius is not positive.
     */
    double centripetal_acceleration(double speed, double radius);

    /**
     * @brief Centripetal force, \f$ F_c = m v^2 / r \f$.
     * @param mass Mass \f$ m \f$ [kg]; must be non-negative.
     * @param speed Speed \f$ v \f$ [m/s].
     * @param radius Radius of the circular path \f$ r \f$ [m]; must be positive.
     * @return Centripetal force [N].
     * @throws std::invalid_argument If @p mass is negative or @p radius is not positive.
     */
    double centripetal_force(double mass, double speed, double radius);

    /**
     * @name Moments of inertia
     * Standard rigid bodies about the stated axis. Every function takes the body's mass [kg] and a
     * characteristic size [m], both non-negative, returns \f$ I \f$ [kg·m²] and throws
     * std::invalid_argument for negative input.
     * @{
     */

    /** @brief Point mass at distance \f$ r \f$: \f$ I = m r^2 \f$. */
    double moment_of_inertia_point_mass(double mass, double radius);
    /** @brief Solid sphere about a diameter: \f$ I = \tfrac25 m r^2 \f$. */
    double moment_of_inertia_solid_sphere(double mass, double radius);
    /** @brief Thin spherical shell about a diameter: \f$ I = \tfrac23 m r^2 \f$. */
    double moment_of_inertia_hollow_sphere(double mass, double radius);
    /** @brief Solid cylinder or disc about its symmetry axis: \f$ I = \tfrac12 m r^2 \f$. */
    double moment_of_inertia_solid_cylinder(double mass, double radius);
    /** @brief Thin hoop or ring about its symmetry axis: \f$ I = m r^2 \f$. */
    double moment_of_inertia_hoop(double mass, double radius);
    /** @brief Thin rod about its centre, perpendicular to its length: \f$ I = \tfrac1{12} m L^2 \f$. */
    double moment_of_inertia_rod_center(double mass, double length);
    /** @brief Thin rod about one end, perpendicular to its length: \f$ I = \tfrac13 m L^2 \f$. */
    double moment_of_inertia_rod_end(double mass, double length);

    /** @} */

    /**
     * @brief Parallel-axis theorem, \f$ I = I_{cm} + m d^2 \f$.
     * @param moment_of_inertia_cm Moment of inertia about a parallel axis through the centre of mass [kg·m²];
     *        must be non-negative.
     * @param mass Body mass \f$ m \f$ [kg]; must be non-negative.
     * @param distance Distance between the two axes \f$ d \f$ [m].
     * @return Moment of inertia about the shifted axis [kg·m²].
     * @throws std::invalid_argument If @p moment_of_inertia_cm or @p mass is negative.
     */
    double parallel_axis(double moment_of_inertia_cm, double mass, double distance);

} // namespace physicsforge::mechanics

#endif // PHYSICSFORGE_MECHANICS_ROTATION_H

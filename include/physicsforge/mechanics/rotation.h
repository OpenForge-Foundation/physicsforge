// include/physicsforge/mechanics/rotation.h
//
// Rigid-body rotation: torque, angular momentum, moments of inertia and circular
// motion. SI units; angles in radians; angular velocity in rad/s.

#ifndef PHYSICSFORGE_MECHANICS_ROTATION_H
#define PHYSICSFORGE_MECHANICS_ROTATION_H

#include "physicsforge/math/vector3.h"

namespace physicsforge::mechanics {

    // Torque magnitude: tau = r F sin(theta), theta measured between r and F. Requires r >= 0.
    double torque(double lever_arm, double force, double angle);

    // Torque vector: tau = r x F.
    math::Vec3 torque_vector(const math::Vec3 &lever_arm, const math::Vec3 &force);

    // Angular acceleration from Newton's second law for rotation: alpha = tau / I. Requires I > 0.
    double angular_acceleration(double torque, double moment_of_inertia);

    // Angular momentum about a fixed axis: L = I omega. Requires I >= 0.
    double angular_momentum(double moment_of_inertia, double angular_velocity);

    // Angular momentum of a particle: L = r x p.
    math::Vec3 angular_momentum_vector(const math::Vec3 &position, const math::Vec3 &linear_momentum);

    // Rotational kinetic energy: K = I omega^2 / 2. Requires I >= 0.
    double rotational_kinetic_energy(double moment_of_inertia, double angular_velocity);

    // Tangential speed of a point at radius r: v = omega r. Requires r >= 0.
    double tangential_velocity(double angular_velocity, double radius);

    // Centripetal acceleration: a_c = v^2 / r. Requires r > 0.
    double centripetal_acceleration(double speed, double radius);

    // Centripetal force: F_c = m v^2 / r. Requires m >= 0, r > 0.
    double centripetal_force(double mass, double speed, double radius);

    // --- Moments of inertia (about the stated axis). All require m >= 0 and size >= 0. ---

    // Point mass at distance r: I = m r^2.
    double moment_of_inertia_point_mass(double mass, double radius);
    // Solid sphere about a diameter: I = 2/5 m r^2.
    double moment_of_inertia_solid_sphere(double mass, double radius);
    // Thin spherical shell about a diameter: I = 2/3 m r^2.
    double moment_of_inertia_hollow_sphere(double mass, double radius);
    // Solid cylinder/disc about its symmetry axis: I = 1/2 m r^2.
    double moment_of_inertia_solid_cylinder(double mass, double radius);
    // Thin hoop/ring about its symmetry axis: I = m r^2.
    double moment_of_inertia_hoop(double mass, double radius);
    // Thin rod about its centre, perpendicular to its length: I = 1/12 m L^2.
    double moment_of_inertia_rod_center(double mass, double length);
    // Thin rod about one end, perpendicular to its length: I = 1/3 m L^2.
    double moment_of_inertia_rod_end(double mass, double length);

    // Parallel-axis theorem: I = I_cm + m d^2. Requires I_cm >= 0, m >= 0.
    double parallel_axis(double moment_of_inertia_cm, double mass, double distance);

} // namespace physicsforge::mechanics

#endif // PHYSICSFORGE_MECHANICS_ROTATION_H

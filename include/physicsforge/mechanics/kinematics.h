// include/physicsforge/mechanics/kinematics.h
//
// One-dimensional constant-acceleration kinematics and ideal projectile motion.
//
// All quantities are SI: metres, seconds, m/s, m/s^2, radians. Functions throw
// std::invalid_argument when a physical precondition is violated.

#ifndef PHYSICSFORGE_MECHANICS_KINEMATICS_H
#define PHYSICSFORGE_MECHANICS_KINEMATICS_H

namespace physicsforge::mechanics {

    // Displacement under constant acceleration: s = v0 t + a t^2 / 2. Requires t >= 0.
    double displacement(double initial_velocity, double acceleration, double time);

    // Velocity under constant acceleration: v = v0 + a t. Requires t >= 0.
    double final_velocity(double initial_velocity, double acceleration, double time);

    // Speed after displacement s under constant acceleration: |v| = sqrt(v0^2 + 2 a s).
    // Throws if v0^2 + 2 a s < 0 (the displacement is never reached).
    double final_velocity_from_displacement(double initial_velocity, double acceleration, double displacement);

    // Average velocity: v_avg = delta_x / delta_t. Requires delta_t > 0.
    double average_velocity(double displacement, double time);

    // Average acceleration: a_avg = delta_v / delta_t. Requires delta_t > 0.
    double average_acceleration(double velocity_change, double time);

    // Horizontal range of a projectile launched from and landing on level ground (no drag):
    // R = v0^2 sin(2 theta) / g. Requires v0 >= 0, g > 0.
    double projectile_range(double initial_speed, double launch_angle, double gravity);

    // Time of flight on level ground: T = 2 v0 sin(theta) / g. Requires v0 >= 0, g > 0.
    double projectile_time_of_flight(double initial_speed, double launch_angle, double gravity);

    // Apex height above the launch point: H = (v0 sin(theta))^2 / (2 g). Requires v0 >= 0, g > 0.
    double projectile_max_height(double initial_speed, double launch_angle, double gravity);

} // namespace physicsforge::mechanics

#endif // PHYSICSFORGE_MECHANICS_KINEMATICS_H

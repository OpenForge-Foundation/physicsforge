// include/physicsforge/mechanics/dynamics.h
//
// Newtonian dynamics: forces, momentum, work and energy.
//
// SI units throughout (kg, m, s, N, J, W). Angles are in radians.

#ifndef PHYSICSFORGE_MECHANICS_DYNAMICS_H
#define PHYSICSFORGE_MECHANICS_DYNAMICS_H

namespace physicsforge::mechanics {

    // Newton's second law: F = m a. Requires m >= 0.
    double net_force(double mass, double acceleration);

    // Acceleration from Newton's second law: a = F / m. Requires m > 0.
    double acceleration_from_force(double force, double mass);

    // Weight (gravitational force near a surface): W = m g. Requires m >= 0, g >= 0.
    double weight(double mass, double gravity);

    // Linear momentum: p = m v. Requires m >= 0.
    double momentum(double mass, double velocity);

    // Impulse of a constant force: J = F dt. Requires dt >= 0.
    double impulse(double force, double duration);

    // Translational kinetic energy: K = m v^2 / 2. Requires m >= 0.
    double kinetic_energy(double mass, double velocity);

    // Uniform-field gravitational potential energy: U = m g h. Requires m >= 0, g >= 0.
    double potential_energy(double mass, double height, double gravity);

    // Work done by a constant force: W = F d cos(theta). Requires d >= 0.
    double work(double force, double distance, double angle);

    // Average power: P = W / t. Requires t > 0.
    double power(double work, double time);

    // Hooke's law restoring force: F = -k x. Requires k >= 0.
    double spring_force(double spring_constant, double extension);

    // Elastic potential energy of an ideal spring: U = k x^2 / 2. Requires k >= 0.
    double spring_potential_energy(double spring_constant, double extension);

    // Kinetic (Coulomb) friction magnitude: f = mu N. Requires mu >= 0, N >= 0.
    double friction_force(double coefficient, double normal_force);

} // namespace physicsforge::mechanics

#endif // PHYSICSFORGE_MECHANICS_DYNAMICS_H

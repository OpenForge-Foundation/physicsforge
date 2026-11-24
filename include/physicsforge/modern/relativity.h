// include/physicsforge/modern/relativity.h
//
// Special relativity: Lorentz factor, time dilation, length contraction and
// relativistic energy/momentum. Speeds must satisfy |v| < c.

#ifndef PHYSICSFORGE_MODERN_RELATIVITY_H
#define PHYSICSFORGE_MODERN_RELATIVITY_H

namespace physicsforge::modern {

    // Lorentz factor: gamma = 1 / sqrt(1 - v^2 / c^2). Requires |v| < c.
    double lorentz_factor(double velocity);

    // Dilated (coordinate) time interval: dt = gamma dt0. Requires dt0 >= 0, |v| < c.
    double time_dilation(double proper_time, double velocity);

    // Contracted length: L = L0 / gamma. Requires L0 >= 0, |v| < c.
    double length_contraction(double proper_length, double velocity);

    // Mass-energy equivalence (rest energy): E = m c^2. Requires m >= 0.
    double rest_energy(double mass);

    // Relativistic momentum: p = gamma m v. Requires m >= 0, |v| < c.
    double relativistic_momentum(double mass, double velocity);

    // Relativistic kinetic energy: K = (gamma - 1) m c^2. Requires m >= 0, |v| < c.
    double relativistic_kinetic_energy(double mass, double velocity);

} // namespace physicsforge::modern

#endif // PHYSICSFORGE_MODERN_RELATIVITY_H

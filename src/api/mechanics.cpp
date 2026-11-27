// src/api/mechanics.cpp - C API shims for physicsforge::mechanics
#include "physicsforge/api.h"
#include "physicsforge/mechanics/dynamics.h"
#include "physicsforge/mechanics/gravitation.h"
#include "physicsforge/mechanics/kinematics.h"
#include "physicsforge/mechanics/rocketry.h"
#include "physicsforge/mechanics/rotation.h"
#include "physicsforge/validation.h"

namespace m = physicsforge::mechanics;
using physicsforge::detail::guarded;

extern "C" double pf_mechanics_displacement(double v0, double a, double t) {
    return guarded([&] { return m::displacement(v0, a, t); });
}
extern "C" double pf_mechanics_final_velocity(double v0, double a, double t) {
    return guarded([&] { return m::final_velocity(v0, a, t); });
}
extern "C" double pf_mechanics_projectile_range(double v0, double angle, double g) {
    return guarded([&] { return m::projectile_range(v0, angle, g); });
}
extern "C" double pf_mechanics_net_force(double mass, double a) {
    return guarded([&] { return m::net_force(mass, a); });
}
extern "C" double pf_mechanics_momentum(double mass, double v) {
    return guarded([&] { return m::momentum(mass, v); });
}
extern "C" double pf_mechanics_kinetic_energy(double mass, double v) {
    return guarded([&] { return m::kinetic_energy(mass, v); });
}
extern "C" double pf_mechanics_work(double force, double distance, double angle) {
    return guarded([&] { return m::work(force, distance, angle); });
}
extern "C" double pf_mechanics_gravitational_force(double m1, double m2, double r) {
    return guarded([&] { return m::gravitational_force(m1, m2, r); });
}
extern "C" double pf_mechanics_escape_velocity(double mass, double r) {
    return guarded([&] { return m::escape_velocity(mass, r); });
}
extern "C" double pf_mechanics_orbital_period(double a, double mass) {
    return guarded([&] { return m::orbital_period(a, mass); });
}
extern "C" double pf_mechanics_torque(double r, double f, double angle) {
    return guarded([&] { return m::torque(r, f, angle); });
}
extern "C" double pf_mechanics_angular_momentum(double inertia, double omega) {
    return guarded([&] { return m::angular_momentum(inertia, omega); });
}
extern "C" double pf_mechanics_rocket_delta_v(double ve, double m0, double mf) {
    return guarded([&] { return m::rocket_delta_v(ve, m0, mf); });
}
extern "C" double pf_mechanics_rocket_mass_ratio(double dv, double ve) {
    return guarded([&] { return m::rocket_mass_ratio(dv, ve); });
}

/**
 * @file api.h
 * @brief Stable C ABI for the scalar physics laws.
 *
 * Every function is a thin, exception-safe shim over the C++ function it names:
 * - on success it returns exactly what the C++ function returns;
 * - on invalid input (where the C++ function would throw) it returns NaN, so test with @c isnan().
 *
 * Naming is @c pf_<domain>_<function>. Units, preconditions and formulas are documented on the C++
 * functions linked from each entry. This header is valid C99 and C++.
 *
 * @defgroup c_api C API
 * @brief Exception-free C entry points.
 * @{
 */

#ifndef PHYSICSFORGE_API_H
#define PHYSICSFORGE_API_H

#include "physicsforge/version.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ---- Mechanics: kinematics ---------------------------------------------------------------- */
/** @brief C wrapper for physicsforge::mechanics::displacement(); NaN on invalid input. */
PF_API double pf_mechanics_displacement(double initial_velocity, double acceleration, double time);
/** @brief C wrapper for physicsforge::mechanics::final_velocity(); NaN on invalid input. */
PF_API double pf_mechanics_final_velocity(double initial_velocity, double acceleration, double time);
/** @brief C wrapper for physicsforge::mechanics::projectile_range(); NaN on invalid input. */
PF_API double pf_mechanics_projectile_range(double initial_speed, double launch_angle, double gravity);

/* ---- Mechanics: dynamics ------------------------------------------------------------------ */
/** @brief C wrapper for physicsforge::mechanics::net_force(); NaN on invalid input. */
PF_API double pf_mechanics_net_force(double mass, double acceleration);
/** @brief C wrapper for physicsforge::mechanics::momentum(); NaN on invalid input. */
PF_API double pf_mechanics_momentum(double mass, double velocity);
/** @brief C wrapper for physicsforge::mechanics::kinetic_energy(); NaN on invalid input. */
PF_API double pf_mechanics_kinetic_energy(double mass, double velocity);
/** @brief C wrapper for physicsforge::mechanics::work(); NaN on invalid input. */
PF_API double pf_mechanics_work(double force, double distance, double angle);

/* ---- Mechanics: gravitation --------------------------------------------------------------- */
/** @brief C wrapper for physicsforge::mechanics::gravitational_force(); NaN on invalid input. */
PF_API double pf_mechanics_gravitational_force(double mass1, double mass2, double distance);
/** @brief C wrapper for physicsforge::mechanics::escape_velocity(); NaN on invalid input. */
PF_API double pf_mechanics_escape_velocity(double central_mass, double distance);
/** @brief C wrapper for physicsforge::mechanics::orbital_period(); NaN on invalid input. */
PF_API double pf_mechanics_orbital_period(double semi_major_axis, double central_mass);

/* ---- Mechanics: rotation ------------------------------------------------------------------ */
/** @brief C wrapper for physicsforge::mechanics::torque(); NaN on invalid input. */
PF_API double pf_mechanics_torque(double lever_arm, double force, double angle);
/** @brief C wrapper for physicsforge::mechanics::angular_momentum(); NaN on invalid input. */
PF_API double pf_mechanics_angular_momentum(double moment_of_inertia, double angular_velocity);

/* ---- Mechanics: rocketry ------------------------------------------------------------------ */
/** @brief C wrapper for physicsforge::mechanics::rocket_delta_v(); NaN on invalid input. */
PF_API double pf_mechanics_rocket_delta_v(double exhaust_velocity, double initial_mass, double final_mass);
/** @brief C wrapper for physicsforge::mechanics::rocket_mass_ratio(); NaN on invalid input. */
PF_API double pf_mechanics_rocket_mass_ratio(double delta_v, double exhaust_velocity);

/* ---- Electromagnetism --------------------------------------------------------------------- */
/** @brief C wrapper for physicsforge::electromagnetism::ohms_law_voltage(); NaN on invalid input. */
PF_API double pf_em_ohms_law_voltage(double current, double resistance);
/** @brief C wrapper for physicsforge::electromagnetism::ohms_law_current(); NaN on invalid input. */
PF_API double pf_em_ohms_law_current(double voltage, double resistance);
/** @brief C wrapper for physicsforge::electromagnetism::voltage_divider(); NaN on invalid input. */
PF_API double pf_em_voltage_divider(double input_voltage, double r1, double r2);
/** @brief C wrapper for physicsforge::electromagnetism::coulomb_force(); NaN on invalid input. */
PF_API double pf_em_coulomb_force(double charge1, double charge2, double distance);
/** @brief C wrapper for physicsforge::electromagnetism::gauss_electric_flux(); NaN on invalid input. */
PF_API double pf_em_gauss_electric_flux(double enclosed_charge);
/** @brief C wrapper for physicsforge::electromagnetism::faraday_emf(); NaN on invalid input. */
PF_API double pf_em_faraday_emf(double flux_rate_of_change, double turns);
/** @brief C wrapper for physicsforge::electromagnetism::magnetic_field_straight_wire(); NaN on invalid input. */
PF_API double pf_em_magnetic_field_straight_wire(double current, double distance);

/* ---- Thermodynamics ----------------------------------------------------------------------- */
/** @brief C wrapper for physicsforge::thermodynamics::conduction_rate(); NaN on invalid input. */
PF_API double pf_thermo_conduction_rate(double thermal_conductivity, double area, double temperature_difference,
                                        double thickness);
/** @brief C wrapper for physicsforge::thermodynamics::newton_cooling_temperature(); NaN on invalid input. */
PF_API double pf_thermo_newton_cooling_temperature(double initial_temperature, double ambient_temperature,
                                                   double cooling_constant, double time);
/** @brief C wrapper for physicsforge::thermodynamics::radiated_power(); NaN on invalid input. */
PF_API double pf_thermo_radiated_power(double emissivity, double area, double temperature);
/** @brief C wrapper for physicsforge::thermodynamics::carnot_efficiency(); NaN on invalid input. */
PF_API double pf_thermo_carnot_efficiency(double cold_temperature, double hot_temperature);

/* ---- Quantum ------------------------------------------------------------------------------ */
/** @brief C wrapper for physicsforge::quantum::photon_energy_from_frequency(); NaN on invalid input. */
PF_API double pf_quantum_photon_energy_from_frequency(double frequency);
/** @brief C wrapper for physicsforge::quantum::planck_spectral_radiance(); NaN on invalid input. */
PF_API double pf_quantum_planck_spectral_radiance(double wavelength, double temperature);

/* ---- Relativity --------------------------------------------------------------------------- */
/** @brief C wrapper for physicsforge::relativity::lorentz_factor(); NaN on invalid input. */
PF_API double pf_relativity_lorentz_factor(double velocity);
/** @brief C wrapper for physicsforge::relativity::rest_energy(); NaN on invalid input. */
PF_API double pf_relativity_rest_energy(double mass);

#ifdef __cplusplus
}
#endif

/** @} */

#endif /* PHYSICSFORGE_API_H */

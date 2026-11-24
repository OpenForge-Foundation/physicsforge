/* include/physicsforge/api.h
 *
 * Stable C ABI for the scalar physics laws in physicsforge.
 *
 * Every function is a thin, exception-safe shim over the C++ core:
 *   - On success it returns the same value as the corresponding C++ function.
 *   - On invalid input (the C++ function would throw) it returns NaN; test with isnan().
 *
 * Naming: pf_<domain>_<function>. Units are SI and angles are radians, exactly as in
 * the C++ headers, which carry the full documentation of each formula.
 *
 * This header is valid C99 and C++.
 */

#ifndef PHYSICSFORGE_API_H
#define PHYSICSFORGE_API_H

#include "physicsforge/export.h"
#include "physicsforge/version.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ---- Mechanics: kinematics ---------------------------------------------------------------- */
PF_API double pf_mechanics_displacement(double initial_velocity, double acceleration, double time);
PF_API double pf_mechanics_final_velocity(double initial_velocity, double acceleration, double time);
PF_API double pf_mechanics_projectile_range(double initial_speed, double launch_angle, double gravity);

/* ---- Mechanics: dynamics ------------------------------------------------------------------ */
PF_API double pf_mechanics_net_force(double mass, double acceleration);
PF_API double pf_mechanics_momentum(double mass, double velocity);
PF_API double pf_mechanics_kinetic_energy(double mass, double velocity);
PF_API double pf_mechanics_work(double force, double distance, double angle);

/* ---- Mechanics: gravitation --------------------------------------------------------------- */
PF_API double pf_mechanics_gravitational_force(double mass1, double mass2, double distance);
PF_API double pf_mechanics_escape_velocity(double central_mass, double distance);
PF_API double pf_mechanics_orbital_period(double semi_major_axis, double central_mass);

/* ---- Mechanics: rotation ------------------------------------------------------------------ */
PF_API double pf_mechanics_torque(double lever_arm, double force, double angle);
PF_API double pf_mechanics_angular_momentum(double moment_of_inertia, double angular_velocity);

/* ---- Mechanics: rocketry ------------------------------------------------------------------ */
PF_API double pf_mechanics_rocket_delta_v(double exhaust_velocity, double initial_mass, double final_mass);
PF_API double pf_mechanics_rocket_mass_ratio(double delta_v, double exhaust_velocity);

/* ---- Electromagnetism --------------------------------------------------------------------- */
PF_API double pf_em_ohms_law_voltage(double current, double resistance);
PF_API double pf_em_ohms_law_current(double voltage, double resistance);
PF_API double pf_em_voltage_divider(double input_voltage, double r1, double r2);
PF_API double pf_em_coulomb_force(double charge1, double charge2, double distance);
PF_API double pf_em_gauss_electric_flux(double enclosed_charge);
PF_API double pf_em_faraday_emf(double flux_rate_of_change, double turns);
PF_API double pf_em_magnetic_field_straight_wire(double current, double distance);

/* ---- Thermodynamics ----------------------------------------------------------------------- */
PF_API double pf_thermo_conduction_rate(double thermal_conductivity, double area, double temperature_difference,
                                        double thickness);
PF_API double pf_thermo_newton_cooling_temperature(double initial_temperature, double ambient_temperature,
                                                   double cooling_constant, double time);
PF_API double pf_thermo_radiated_power(double emissivity, double area, double temperature);
PF_API double pf_thermo_carnot_efficiency(double cold_temperature, double hot_temperature);

/* ---- Modern physics ----------------------------------------------------------------------- */
PF_API double pf_modern_photon_energy_from_frequency(double frequency);
PF_API double pf_modern_planck_spectral_radiance(double wavelength, double temperature);
PF_API double pf_modern_lorentz_factor(double velocity);
PF_API double pf_modern_rest_energy(double mass);

#ifdef __cplusplus
}
#endif

#endif /* PHYSICSFORGE_API_H */

// include/physicsforge/thermodynamics/heat.h
//
// Heat transfer (conduction, convection, radiation), calorimetry and the laws of
// thermodynamics. Temperatures are absolute (kelvin) unless a function explicitly
// works with temperature *differences*, which are identical in K and degC.

#ifndef PHYSICSFORGE_THERMODYNAMICS_HEAT_H
#define PHYSICSFORGE_THERMODYNAMICS_HEAT_H

namespace physicsforge::thermodynamics {

    // --- Calorimetry -----------------------------------------------------------------------------

    // Sensible heat: Q = m c dT. Requires m >= 0, c >= 0.
    double sensible_heat(double mass, double specific_heat, double temperature_change);

    // Latent heat of a phase change: Q = m L. Requires m >= 0, L >= 0.
    double latent_heat(double mass, double specific_latent_heat);

    // --- Conduction (Fourier's law) --------------------------------------------------------------

    // Fourier's law, 1-D heat flux: q = -k dT/dx [W m^-2]. Requires k >= 0.
    double fourier_heat_flux(double thermal_conductivity, double temperature_gradient);

    // Steady conduction rate through a plane wall: Qdot = k A (T_hot - T_cold) / L [W].
    // Requires k >= 0, A >= 0, L > 0.
    double conduction_rate(double thermal_conductivity, double area, double temperature_difference, double thickness);

    // Conductive thermal resistance of a plane wall: R = L / (k A) [K W^-1]. Requires L >= 0, k > 0, A > 0.
    double conduction_resistance(double thickness, double thermal_conductivity, double area);

    // Convective thermal resistance: R = 1 / (h A) [K W^-1]. Requires h > 0, A > 0.
    double convection_resistance(double heat_transfer_coefficient, double area);

    // Heat rate across a thermal resistance network: Qdot = dT / R_total. Requires R > 0.
    double heat_rate_from_resistance(double temperature_difference, double thermal_resistance);

    // Thermal diffusivity: alpha = k / (rho c_p) [m^2 s^-1]. Requires k >= 0, rho > 0, c_p > 0.
    double thermal_diffusivity(double thermal_conductivity, double density, double specific_heat);

    // --- Convection (Newton's law of cooling) ----------------------------------------------------

    // Convective heat rate: Qdot = h A (T_s - T_inf). Requires h >= 0, A >= 0.
    double newton_cooling_rate(double heat_transfer_coefficient, double area, double surface_temperature,
                               double ambient_temperature);

    // Lumped-capacitance temperature history: T(t) = T_env + (T0 - T_env) exp(-k t).
    // Requires k >= 0, t >= 0.
    double newton_cooling_temperature(double initial_temperature, double ambient_temperature, double cooling_constant,
                                      double time);

    // --- Radiation (Stefan-Boltzmann law) --------------------------------------------------------

    // Emitted radiant power of a grey body: P = eps sigma A T^4. Requires 0 <= eps <= 1, A >= 0, T >= 0 K.
    double radiated_power(double emissivity, double area, double temperature);

    // Net radiative exchange with large surroundings: P = eps sigma A (T^4 - T_surr^4).
    double net_radiated_power(double emissivity, double area, double temperature, double surroundings_temperature);

    // --- Laws of thermodynamics and ideal gases --------------------------------------------------

    // First law: dU = Q - W, with W the work done *by* the system.
    double internal_energy_change(double heat_added, double work_done_by_system);

    // Carnot (maximum) efficiency of a heat engine: eta = 1 - T_c / T_h. Requires 0 < T_c <= T_h.
    double carnot_efficiency(double cold_temperature, double hot_temperature);

    // Carnot coefficient of performance of a refrigerator: COP = T_c / (T_h - T_c). Requires 0 < T_c < T_h.
    double carnot_cop_refrigerator(double cold_temperature, double hot_temperature);

    // Ideal gas law solved for pressure: P = n R T / V. Requires n >= 0, T >= 0, V > 0.
    double ideal_gas_pressure(double amount, double temperature, double volume);

    // Ideal gas law solved for volume: V = n R T / P. Requires n >= 0, T >= 0, P > 0.
    double ideal_gas_volume(double amount, double temperature, double pressure);

    // Linear thermal expansion: dL = alpha L0 dT. Requires L0 >= 0.
    double linear_thermal_expansion(double expansion_coefficient, double initial_length, double temperature_change);

    // Unit conversions between Celsius and kelvin (exact offset 273.15).
    double celsius_to_kelvin(double celsius);
    double kelvin_to_celsius(double kelvin);

} // namespace physicsforge::thermodynamics

#endif // PHYSICSFORGE_THERMODYNAMICS_HEAT_H

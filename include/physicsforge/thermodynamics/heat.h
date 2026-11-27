/**
 * @file heat.h
 * @brief Heat transfer (conduction, convection, radiation), calorimetry and the laws of thermodynamics.
 *
 * Temperatures are absolute (K) unless a parameter is a temperature *difference*, which is identical
 * in K and °C.
 *
 * @ingroup thermodynamics
 */

#ifndef PHYSICSFORGE_THERMODYNAMICS_HEAT_H
#define PHYSICSFORGE_THERMODYNAMICS_HEAT_H

/** @defgroup thermodynamics Thermodynamics
 *  @brief Heat laws and the one-dimensional heat equation.
 */

namespace physicsforge::thermodynamics {

    /** @name Calorimetry
     *  @{ */

    /**
     * @brief Sensible heat, \f$ Q = m c \Delta T \f$.
     * @param mass Mass \f$ m \f$ [kg]; must be non-negative.
     * @param specific_heat Specific heat capacity \f$ c \f$ [J/(kg·K)]; must be non-negative.
     * @param temperature_change Temperature change \f$ \Delta T \f$ [K].
     * @return Heat \f$ Q \f$ [J].
     * @throws std::invalid_argument If @p mass or @p specific_heat is negative.
     */
    double sensible_heat(double mass, double specific_heat, double temperature_change);

    /**
     * @brief Latent heat of a phase change, \f$ Q = m L \f$.
     * @param mass Mass \f$ m \f$ [kg]; must be non-negative.
     * @param specific_latent_heat Specific latent heat \f$ L \f$ [J/kg]; must be non-negative.
     * @return Heat \f$ Q \f$ [J].
     * @throws std::invalid_argument If either input is negative.
     */
    double latent_heat(double mass, double specific_latent_heat);

    /** @} */

    /** @name Conduction (Fourier's law)
     *  @{ */

    /**
     * @brief Fourier's law of conduction in one dimension, \f$ q = -k\, dT/dx \f$.
     * @param thermal_conductivity Thermal conductivity \f$ k \f$ [W/(m·K)]; must be non-negative.
     * @param temperature_gradient Temperature gradient \f$ dT/dx \f$ [K/m].
     * @return Heat flux \f$ q \f$ [W/m²].
     * @throws std::invalid_argument If @p thermal_conductivity is negative.
     */
    double fourier_heat_flux(double thermal_conductivity, double temperature_gradient);

    /**
     * @brief Steady conduction rate through a plane wall, \f$ \dot Q = k A \Delta T / L \f$.
     * @param thermal_conductivity Thermal conductivity \f$ k \f$ [W/(m·K)]; must be non-negative.
     * @param area Wall area \f$ A \f$ [m²]; must be non-negative.
     * @param temperature_difference Hot-side minus cold-side temperature \f$ \Delta T \f$ [K].
     * @param thickness Wall thickness \f$ L \f$ [m]; must be positive.
     * @return Heat rate \f$ \dot Q \f$ [W].
     * @throws std::invalid_argument If any precondition above is violated.
     */
    double conduction_rate(double thermal_conductivity, double area, double temperature_difference, double thickness);

    /**
     * @brief Conductive thermal resistance of a plane wall, \f$ R = L / (k A) \f$.
     * @param thickness Wall thickness \f$ L \f$ [m]; must be non-negative.
     * @param thermal_conductivity Thermal conductivity \f$ k \f$ [W/(m·K)]; must be positive.
     * @param area Wall area \f$ A \f$ [m²]; must be positive.
     * @return Thermal resistance [K/W].
     * @throws std::invalid_argument If any precondition above is violated.
     */
    double conduction_resistance(double thickness, double thermal_conductivity, double area);

    /**
     * @brief Convective thermal resistance, \f$ R = 1 / (h A) \f$.
     * @param heat_transfer_coefficient Convection coefficient \f$ h \f$ [W/(m²·K)]; must be positive.
     * @param area Surface area \f$ A \f$ [m²]; must be positive.
     * @return Thermal resistance [K/W].
     * @throws std::invalid_argument If either input is not positive.
     */
    double convection_resistance(double heat_transfer_coefficient, double area);

    /**
     * @brief Heat rate across a thermal resistance network, \f$ \dot Q = \Delta T / R_{total} \f$.
     * @param temperature_difference Overall temperature difference \f$ \Delta T \f$ [K].
     * @param thermal_resistance Total thermal resistance \f$ R_{total} \f$ [K/W]; must be positive.
     * @return Heat rate [W].
     * @throws std::invalid_argument If @p thermal_resistance is not positive.
     */
    double heat_rate_from_resistance(double temperature_difference, double thermal_resistance);

    /**
     * @brief Thermal diffusivity, \f$ \alpha = k / (\rho c_p) \f$.
     * @param thermal_conductivity Thermal conductivity \f$ k \f$ [W/(m·K)]; must be non-negative.
     * @param density Density \f$ \rho \f$ [kg/m³]; must be positive.
     * @param specific_heat Specific heat capacity \f$ c_p \f$ [J/(kg·K)]; must be positive.
     * @return Thermal diffusivity \f$ \alpha \f$ [m²/s].
     * @throws std::invalid_argument If any precondition above is violated.
     */
    double thermal_diffusivity(double thermal_conductivity, double density, double specific_heat);

    /** @} */

    /** @name Convection (Newton's law of cooling)
     *  @{ */

    /**
     * @brief Convective heat rate, \f$ \dot Q = h A (T_s - T_\infty) \f$.
     * @param heat_transfer_coefficient Convection coefficient \f$ h \f$ [W/(m²·K)]; must be non-negative.
     * @param area Surface area \f$ A \f$ [m²]; must be non-negative.
     * @param surface_temperature Surface temperature \f$ T_s \f$ [K].
     * @param ambient_temperature Fluid temperature far from the surface \f$ T_\infty \f$ [K].
     * @return Heat rate leaving the surface [W].
     * @throws std::invalid_argument If @p heat_transfer_coefficient or @p area is negative.
     */
    double newton_cooling_rate(double heat_transfer_coefficient, double area, double surface_temperature,
                               double ambient_temperature);

    /**
     * @brief Lumped-capacitance temperature history,
     *        \f$ T(t) = T_{env} + (T_0 - T_{env}) e^{-k t} \f$.
     * @param initial_temperature Initial body temperature \f$ T_0 \f$ [K].
     * @param ambient_temperature Environment temperature \f$ T_{env} \f$ [K].
     * @param cooling_constant Cooling constant \f$ k = hA/(mc) \f$ [1/s]; must be non-negative.
     * @param time Elapsed time \f$ t \f$ [s]; must be non-negative.
     * @return Body temperature at @p time [K].
     * @throws std::invalid_argument If @p cooling_constant or @p time is negative.
     */
    double newton_cooling_temperature(double initial_temperature, double ambient_temperature, double cooling_constant,
                                      double time);

    /** @} */

    /** @name Radiation (Stefan–Boltzmann law)
     *  @{ */

    /**
     * @brief Radiant power emitted by a grey body, \f$ P = \varepsilon \sigma A T^4 \f$.
     * @param emissivity Emissivity \f$ \varepsilon \f$; must lie in [0, 1].
     * @param area Radiating area \f$ A \f$ [m²]; must be non-negative.
     * @param temperature Absolute surface temperature \f$ T \f$ [K]; must be non-negative.
     * @return Emitted power [W].
     * @throws std::invalid_argument If any precondition above is violated.
     */
    double radiated_power(double emissivity, double area, double temperature);

    /**
     * @brief Net radiative exchange with large surroundings, \f$ P = \varepsilon \sigma A (T^4 - T_{surr}^4) \f$.
     * @param emissivity Emissivity \f$ \varepsilon \f$; must lie in [0, 1].
     * @param area Radiating area \f$ A \f$ [m²]; must be non-negative.
     * @param temperature Absolute surface temperature \f$ T \f$ [K]; must be non-negative.
     * @param surroundings_temperature Absolute temperature of the surroundings [K]; must be non-negative.
     * @return Net power leaving the surface [W].
     * @throws std::invalid_argument If any precondition above is violated.
     */
    double net_radiated_power(double emissivity, double area, double temperature, double surroundings_temperature);

    /** @} */

    /** @name Laws of thermodynamics and ideal gases
     *  @{ */

    /**
     * @brief First law of thermodynamics, \f$ \Delta U = Q - W \f$.
     * @param heat_added Heat added to the system \f$ Q \f$ [J].
     * @param work_done_by_system Work done *by* the system \f$ W \f$ [J].
     * @return Change in internal energy [J].
     */
    double internal_energy_change(double heat_added, double work_done_by_system);

    /**
     * @brief Carnot (maximum) efficiency of a heat engine, \f$ \eta = 1 - T_c / T_h \f$.
     * @param cold_temperature Cold reservoir temperature \f$ T_c \f$ [K]; must be positive.
     * @param hot_temperature Hot reservoir temperature \f$ T_h \f$ [K]; must satisfy \f$ T_h \ge T_c \f$.
     * @return Efficiency in [0, 1).
     * @throws std::invalid_argument If \f$ T_c \le 0 \f$ or \f$ T_c > T_h \f$.
     */
    double carnot_efficiency(double cold_temperature, double hot_temperature);

    /**
     * @brief Carnot coefficient of performance of a refrigerator, \f$ COP = T_c / (T_h - T_c) \f$.
     * @param cold_temperature Cold reservoir temperature \f$ T_c \f$ [K]; must be positive.
     * @param hot_temperature Hot reservoir temperature \f$ T_h \f$ [K]; must satisfy \f$ T_h > T_c \f$.
     * @return Coefficient of performance.
     * @throws std::invalid_argument If \f$ T_c \le 0 \f$ or \f$ T_c \ge T_h \f$.
     */
    double carnot_cop_refrigerator(double cold_temperature, double hot_temperature);

    /**
     * @brief Ideal gas law solved for pressure, \f$ P = n R T / V \f$.
     * @param amount Amount of substance \f$ n \f$ [mol]; must be non-negative.
     * @param temperature Absolute temperature \f$ T \f$ [K]; must be non-negative.
     * @param volume Volume \f$ V \f$ [m³]; must be positive.
     * @return Pressure [Pa].
     * @throws std::invalid_argument If any precondition above is violated.
     */
    double ideal_gas_pressure(double amount, double temperature, double volume);

    /**
     * @brief Ideal gas law solved for volume, \f$ V = n R T / P \f$.
     * @param amount Amount of substance \f$ n \f$ [mol]; must be non-negative.
     * @param temperature Absolute temperature \f$ T \f$ [K]; must be non-negative.
     * @param pressure Pressure \f$ P \f$ [Pa]; must be positive.
     * @return Volume [m³].
     * @throws std::invalid_argument If any precondition above is violated.
     */
    double ideal_gas_volume(double amount, double temperature, double pressure);

    /**
     * @brief Linear thermal expansion, \f$ \Delta L = \alpha L_0 \Delta T \f$.
     * @param expansion_coefficient Linear expansion coefficient \f$ \alpha \f$ [1/K].
     * @param initial_length Initial length \f$ L_0 \f$ [m]; must be non-negative.
     * @param temperature_change Temperature change \f$ \Delta T \f$ [K].
     * @return Change in length [m].
     * @throws std::invalid_argument If @p initial_length is negative.
     */
    double linear_thermal_expansion(double expansion_coefficient, double initial_length, double temperature_change);

    /**
     * @brief Converts Celsius to kelvin, \f$ T_K = T_C + 273.15 \f$.
     * @param celsius Temperature [°C]; must be at or above absolute zero.
     * @return Temperature [K].
     * @throws std::invalid_argument If the result would be below 0 K.
     */
    double celsius_to_kelvin(double celsius);

    /**
     * @brief Converts kelvin to Celsius, \f$ T_C = T_K - 273.15 \f$.
     * @param kelvin Temperature [K]; must be non-negative.
     * @return Temperature [°C].
     * @throws std::invalid_argument If @p kelvin is negative.
     */
    double kelvin_to_celsius(double kelvin);

    /** @} */

} // namespace physicsforge::thermodynamics

#endif // PHYSICSFORGE_THERMODYNAMICS_HEAT_H

/**
 * @file circuits.h
 * @brief Lumped-element DC circuits and first-order RC transients.
 *
 * SI units: V, A, Ω, F, W, s.
 *
 * @ingroup electromagnetism
 */

#ifndef PHYSICSFORGE_ELECTROMAGNETISM_CIRCUITS_H
#define PHYSICSFORGE_ELECTROMAGNETISM_CIRCUITS_H

#include <span>

/** @defgroup electromagnetism Electromagnetism
 *  @brief Circuits, electrostatics and Maxwell's equations.
 */

namespace physicsforge::electromagnetism {

    /**
     * @brief Ohm's law solved for voltage, \f$ V = I R \f$.
     * @param current Current \f$ I \f$ [A].
     * @param resistance Resistance \f$ R \f$ [Ω]; must be non-negative.
     * @return Voltage \f$ V \f$ [V].
     * @throws std::invalid_argument If @p resistance is negative.
     */
    double ohms_law_voltage(double current, double resistance);

    /**
     * @brief Ohm's law solved for current, \f$ I = V / R \f$.
     * @param voltage Voltage \f$ V \f$ [V].
     * @param resistance Resistance \f$ R \f$ [Ω]; must be positive.
     * @return Current \f$ I \f$ [A].
     * @throws std::invalid_argument If @p resistance is not positive.
     */
    double ohms_law_current(double voltage, double resistance);

    /**
     * @brief Ohm's law solved for resistance, \f$ R = V / I \f$.
     * @param voltage Voltage \f$ V \f$ [V].
     * @param current Current \f$ I \f$ [A]; must be non-zero.
     * @return Resistance \f$ R \f$ [Ω].
     * @throws std::invalid_argument If @p current is zero.
     */
    double ohms_law_resistance(double voltage, double current);

    /**
     * @brief Electrical power, \f$ P = V I \f$.
     * @param voltage Voltage \f$ V \f$ [V].
     * @param current Current \f$ I \f$ [A].
     * @return Power \f$ P \f$ [W].
     */
    double electrical_power(double voltage, double current);

    /**
     * @brief Joule heating from current, \f$ P = I^2 R \f$.
     * @param current Current \f$ I \f$ [A].
     * @param resistance Resistance \f$ R \f$ [Ω]; must be non-negative.
     * @return Dissipated power [W].
     * @throws std::invalid_argument If @p resistance is negative.
     */
    double power_from_current(double current, double resistance);

    /**
     * @brief Joule heating from voltage, \f$ P = V^2 / R \f$.
     * @param voltage Voltage \f$ V \f$ [V].
     * @param resistance Resistance \f$ R \f$ [Ω]; must be positive.
     * @return Dissipated power [W].
     * @throws std::invalid_argument If @p resistance is not positive.
     */
    double power_from_voltage(double voltage, double resistance);

    /**
     * @brief Unloaded voltage divider with the output taken across @p r2,
     *        \f$ V_{out} = V_{in} R_2 / (R_1 + R_2) \f$.
     * @param input_voltage Input voltage \f$ V_{in} \f$ [V].
     * @param r1 Upper resistor \f$ R_1 \f$ [Ω]; must be non-negative.
     * @param r2 Lower resistor \f$ R_2 \f$ [Ω]; must be non-negative.
     * @return Output voltage \f$ V_{out} \f$ [V].
     * @throws std::invalid_argument If a resistance is negative or both are zero.
     */
    double voltage_divider(double input_voltage, double r1, double r2);

    /**
     * @brief Two-branch current divider, \f$ I_b = I_{total} R_{other} / (R_b + R_{other}) \f$.
     * @param total_current Total current entering the junction [A].
     * @param r_branch Resistance of the branch of interest \f$ R_b \f$ [Ω]; must be non-negative.
     * @param r_other Resistance of the other branch \f$ R_{other} \f$ [Ω]; must be non-negative.
     * @return Current through @p r_branch [A].
     * @throws std::invalid_argument If a resistance is negative or both are zero.
     */
    double current_divider(double total_current, double r_branch, double r_other);

    /**
     * @brief Equivalent series resistance, \f$ R = \sum_i R_i \f$.
     * @param resistances Resistances [Ω]; non-empty, each non-negative.
     * @return Equivalent resistance [Ω].
     * @throws std::invalid_argument If the list is empty or a value is negative.
     */
    double series_resistance(std::span<const double> resistances);

    /**
     * @brief Equivalent parallel resistance, \f$ 1/R = \sum_i 1/R_i \f$.
     * @param resistances Resistances [Ω]; non-empty, each positive.
     * @return Equivalent resistance [Ω].
     * @throws std::invalid_argument If the list is empty or a value is not positive.
     */
    double parallel_resistance(std::span<const double> resistances);

    /**
     * @brief Equivalent series capacitance, \f$ 1/C = \sum_i 1/C_i \f$.
     * @param capacitances Capacitances [F]; non-empty, each positive.
     * @return Equivalent capacitance [F].
     * @throws std::invalid_argument If the list is empty or a value is not positive.
     */
    double series_capacitance(std::span<const double> capacitances);

    /**
     * @brief Equivalent parallel capacitance, \f$ C = \sum_i C_i \f$.
     * @param capacitances Capacitances [F]; non-empty, each non-negative.
     * @return Equivalent capacitance [F].
     * @throws std::invalid_argument If the list is empty or a value is negative.
     */
    double parallel_capacitance(std::span<const double> capacitances);

    /**
     * @brief Resistance of a uniform conductor, \f$ R = \rho L / A \f$.
     * @param resistivity Resistivity \f$ \rho \f$ [Ω·m]; must be non-negative.
     * @param length Conductor length \f$ L \f$ [m]; must be non-negative.
     * @param area Cross-sectional area \f$ A \f$ [m²]; must be positive.
     * @return Resistance [Ω].
     * @throws std::invalid_argument If @p resistivity or @p length is negative or @p area is not positive.
     */
    double resistance_from_resistivity(double resistivity, double length, double area);

    /**
     * @brief RC time constant, \f$ \tau = R C \f$.
     * @param resistance Resistance \f$ R \f$ [Ω]; must be non-negative.
     * @param capacitance Capacitance \f$ C \f$ [F]; must be non-negative.
     * @return Time constant \f$ \tau \f$ [s].
     * @throws std::invalid_argument If either input is negative.
     */
    double rc_time_constant(double resistance, double capacitance);

    /**
     * @brief Voltage across a capacitor charging from 0 V, \f$ V(t) = V_s (1 - e^{-t/RC}) \f$.
     * @param source_voltage Source voltage \f$ V_s \f$ [V].
     * @param time Time since the switch closed \f$ t \f$ [s]; must be non-negative.
     * @param resistance Series resistance \f$ R \f$ [Ω]; must be positive.
     * @param capacitance Capacitance \f$ C \f$ [F]; must be positive.
     * @return Capacitor voltage [V].
     * @throws std::invalid_argument If @p time is negative or @p resistance / @p capacitance is not positive.
     */
    double capacitor_charging_voltage(double source_voltage, double time, double resistance, double capacitance);

    /**
     * @brief Voltage across a discharging capacitor, \f$ V(t) = V_0 e^{-t/RC} \f$.
     * @param initial_voltage Initial voltage \f$ V_0 \f$ [V].
     * @param time Time since discharge began \f$ t \f$ [s]; must be non-negative.
     * @param resistance Discharge resistance \f$ R \f$ [Ω]; must be positive.
     * @param capacitance Capacitance \f$ C \f$ [F]; must be positive.
     * @return Capacitor voltage [V].
     * @throws std::invalid_argument If @p time is negative or @p resistance / @p capacitance is not positive.
     */
    double capacitor_discharging_voltage(double initial_voltage, double time, double resistance, double capacitance);

} // namespace physicsforge::electromagnetism

#endif // PHYSICSFORGE_ELECTROMAGNETISM_CIRCUITS_H

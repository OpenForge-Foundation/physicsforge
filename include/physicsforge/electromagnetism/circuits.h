// include/physicsforge/electromagnetism/circuits.h
//
// Lumped-element DC and first-order RC circuit analysis.
//
// SI units: volts (V), amperes (A), ohms (Ohm), farads (F), watts (W), seconds (s).

#ifndef PHYSICSFORGE_ELECTROMAGNETISM_CIRCUITS_H
#define PHYSICSFORGE_ELECTROMAGNETISM_CIRCUITS_H

#include <span>

namespace physicsforge::electromagnetism {

    // Ohm's law, solved for voltage: V = I R. Requires R >= 0.
    double ohms_law_voltage(double current, double resistance);

    // Ohm's law, solved for current: I = V / R. Requires R > 0.
    double ohms_law_current(double voltage, double resistance);

    // Ohm's law, solved for resistance: R = V / I. Requires I != 0.
    double ohms_law_resistance(double voltage, double current);

    // Electrical power: P = V I.
    double electrical_power(double voltage, double current);

    // Joule heating from current: P = I^2 R. Requires R >= 0.
    double power_from_current(double current, double resistance);

    // Joule heating from voltage: P = V^2 / R. Requires R > 0.
    double power_from_voltage(double voltage, double resistance);

    // Unloaded voltage divider, output taken across r2: V_out = V_in r2 / (r1 + r2).
    // Requires r1, r2 >= 0 and r1 + r2 > 0.
    double voltage_divider(double input_voltage, double r1, double r2);

    // Two-branch current divider, current through r_branch: I_b = I_total r_other / (r_branch + r_other).
    // Requires r_branch, r_other >= 0 and r_branch + r_other > 0.
    double current_divider(double total_current, double r_branch, double r_other);

    // Equivalent resistance in series: R = sum R_i. Requires every R_i >= 0 and a non-empty list.
    double series_resistance(std::span<const double> resistances);

    // Equivalent resistance in parallel: 1/R = sum 1/R_i. Requires every R_i > 0 and a non-empty list.
    double parallel_resistance(std::span<const double> resistances);

    // Equivalent capacitance in series: 1/C = sum 1/C_i. Requires every C_i > 0.
    double series_capacitance(std::span<const double> capacitances);

    // Equivalent capacitance in parallel: C = sum C_i. Requires every C_i >= 0.
    double parallel_capacitance(std::span<const double> capacitances);

    // Resistance of a uniform conductor: R = rho L / A. Requires rho >= 0, L >= 0, A > 0.
    double resistance_from_resistivity(double resistivity, double length, double area);

    // RC time constant: tau = R C. Requires R >= 0, C >= 0.
    double rc_time_constant(double resistance, double capacitance);

    // Capacitor voltage while charging from 0 V through R: V(t) = V_s (1 - exp(-t / RC)).
    // Requires t >= 0, R > 0, C > 0.
    double capacitor_charging_voltage(double source_voltage, double time, double resistance, double capacitance);

    // Capacitor voltage while discharging through R: V(t) = V_0 exp(-t / RC).
    // Requires t >= 0, R > 0, C > 0.
    double capacitor_discharging_voltage(double initial_voltage, double time, double resistance, double capacitance);

} // namespace physicsforge::electromagnetism

#endif // PHYSICSFORGE_ELECTROMAGNETISM_CIRCUITS_H

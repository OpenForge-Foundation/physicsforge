// src/electromagnetism/circuits.cpp
#include "physicsforge/electromagnetism/circuits.h"

#include <cmath>

#include "physicsforge/detail/validation.h"

namespace physicsforge::electromagnetism {

    using detail::require_finite;
    using detail::require_non_negative;
    using detail::require_non_zero;
    using detail::require_positive;

    namespace {
        void require_non_empty(const std::span<const double> values, const char *name) {
            if (values.empty()) {
                detail::throw_invalid(std::string(name) + " must contain at least one element");
            }
        }

        // Sum of values (each must be >= 0).
        double sum_non_negative(const std::span<const double> values, const char *name) {
            require_non_empty(values, name);
            double total = 0.0;
            for (const double v: values) {
                require_non_negative(v, name);
                total += v;
            }
            return total;
        }

        // Reciprocal of the sum of reciprocals (each must be > 0).
        double reciprocal_sum(const std::span<const double> values, const char *name) {
            require_non_empty(values, name);
            double total = 0.0;
            for (const double v: values) {
                require_positive(v, name);
                total += 1.0 / v;
            }
            return 1.0 / total;
        }

        void require_divider_pair(const double a, const double b, const char *name_a, const char *name_b) {
            require_non_negative(a, name_a);
            require_non_negative(b, name_b);
            if (a + b <= 0.0) {
                detail::throw_invalid(std::string(name_a) + " + " + name_b + " must be positive");
            }
        }
    } // namespace

    double ohms_law_voltage(const double current, const double resistance) {
        require_finite(current, "current");
        require_non_negative(resistance, "resistance");
        return current * resistance;
    }

    double ohms_law_current(const double voltage, const double resistance) {
        require_finite(voltage, "voltage");
        require_positive(resistance, "resistance");
        return voltage / resistance;
    }

    double ohms_law_resistance(const double voltage, const double current) {
        require_finite(voltage, "voltage");
        require_non_zero(current, "current");
        return voltage / current;
    }

    double electrical_power(const double voltage, const double current) {
        require_finite(voltage, "voltage");
        require_finite(current, "current");
        return voltage * current;
    }

    double power_from_current(const double current, const double resistance) {
        require_finite(current, "current");
        require_non_negative(resistance, "resistance");
        return current * current * resistance;
    }

    double power_from_voltage(const double voltage, const double resistance) {
        require_finite(voltage, "voltage");
        require_positive(resistance, "resistance");
        return voltage * voltage / resistance;
    }

    double voltage_divider(const double input_voltage, const double r1, const double r2) {
        require_finite(input_voltage, "input_voltage");
        require_divider_pair(r1, r2, "r1", "r2");
        return input_voltage * r2 / (r1 + r2);
    }

    double current_divider(const double total_current, const double r_branch, const double r_other) {
        require_finite(total_current, "total_current");
        require_divider_pair(r_branch, r_other, "r_branch", "r_other");
        return total_current * r_other / (r_branch + r_other);
    }

    double series_resistance(const std::span<const double> resistances) {
        return sum_non_negative(resistances, "resistances");
    }

    double parallel_resistance(const std::span<const double> resistances) {
        return reciprocal_sum(resistances, "resistances");
    }

    double series_capacitance(const std::span<const double> capacitances) {
        return reciprocal_sum(capacitances, "capacitances");
    }

    double parallel_capacitance(const std::span<const double> capacitances) {
        return sum_non_negative(capacitances, "capacitances");
    }

    double resistance_from_resistivity(const double resistivity, const double length, const double area) {
        require_non_negative(resistivity, "resistivity");
        require_non_negative(length, "length");
        require_positive(area, "area");
        return resistivity * length / area;
    }

    double rc_time_constant(const double resistance, const double capacitance) {
        require_non_negative(resistance, "resistance");
        require_non_negative(capacitance, "capacitance");
        return resistance * capacitance;
    }

    double capacitor_charging_voltage(const double source_voltage, const double time, const double resistance,
                                      const double capacitance) {
        require_finite(source_voltage, "source_voltage");
        require_non_negative(time, "time");
        require_positive(resistance, "resistance");
        require_positive(capacitance, "capacitance");
        // -expm1(-x) == 1 - exp(-x), accurate for small t / RC.
        return -source_voltage * std::expm1(-time / (resistance * capacitance));
    }

    double capacitor_discharging_voltage(const double initial_voltage, const double time, const double resistance,
                                         const double capacitance) {
        require_finite(initial_voltage, "initial_voltage");
        require_non_negative(time, "time");
        require_positive(resistance, "resistance");
        require_positive(capacitance, "capacitance");
        return initial_voltage * std::exp(-time / (resistance * capacitance));
    }

} // namespace physicsforge::electromagnetism

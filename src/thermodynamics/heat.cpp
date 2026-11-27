// src/thermodynamics/heat.cpp
#include "physicsforge/thermodynamics/heat.h"

#include <cmath>

#include "physicsforge/constants.h"
#include "physicsforge/validation.h"

namespace physicsforge::thermodynamics {

    using constants::GAS_CONSTANT;
    using constants::STEFAN_BOLTZMANN;
    using detail::require_finite;
    using detail::require_in_range;
    using detail::require_non_negative;
    using detail::require_positive;

    namespace {
        constexpr double CELSIUS_OFFSET = 273.15;

        double fourth_power(const double x) {
            const double sq = x * x;
            return sq * sq;
        }
    } // namespace

    double sensible_heat(const double mass, const double specific_heat, const double temperature_change) {
        require_non_negative(mass, "mass");
        require_non_negative(specific_heat, "specific_heat");
        require_finite(temperature_change, "temperature_change");
        return mass * specific_heat * temperature_change;
    }

    double latent_heat(const double mass, const double specific_latent_heat) {
        require_non_negative(mass, "mass");
        require_non_negative(specific_latent_heat, "specific_latent_heat");
        return mass * specific_latent_heat;
    }

    double fourier_heat_flux(const double thermal_conductivity, const double temperature_gradient) {
        require_non_negative(thermal_conductivity, "thermal_conductivity");
        require_finite(temperature_gradient, "temperature_gradient");
        return -thermal_conductivity * temperature_gradient;
    }

    double conduction_rate(const double thermal_conductivity, const double area, const double temperature_difference,
                           const double thickness) {
        require_non_negative(thermal_conductivity, "thermal_conductivity");
        require_non_negative(area, "area");
        require_finite(temperature_difference, "temperature_difference");
        require_positive(thickness, "thickness");
        return thermal_conductivity * area * temperature_difference / thickness;
    }

    double conduction_resistance(const double thickness, const double thermal_conductivity, const double area) {
        require_non_negative(thickness, "thickness");
        require_positive(thermal_conductivity, "thermal_conductivity");
        require_positive(area, "area");
        return thickness / (thermal_conductivity * area);
    }

    double convection_resistance(const double heat_transfer_coefficient, const double area) {
        require_positive(heat_transfer_coefficient, "heat_transfer_coefficient");
        require_positive(area, "area");
        return 1.0 / (heat_transfer_coefficient * area);
    }

    double heat_rate_from_resistance(const double temperature_difference, const double thermal_resistance) {
        require_finite(temperature_difference, "temperature_difference");
        require_positive(thermal_resistance, "thermal_resistance");
        return temperature_difference / thermal_resistance;
    }

    double thermal_diffusivity(const double thermal_conductivity, const double density, const double specific_heat) {
        require_non_negative(thermal_conductivity, "thermal_conductivity");
        require_positive(density, "density");
        require_positive(specific_heat, "specific_heat");
        return thermal_conductivity / (density * specific_heat);
    }

    double newton_cooling_rate(const double heat_transfer_coefficient, const double area,
                               const double surface_temperature, const double ambient_temperature) {
        require_non_negative(heat_transfer_coefficient, "heat_transfer_coefficient");
        require_non_negative(area, "area");
        require_finite(surface_temperature, "surface_temperature");
        require_finite(ambient_temperature, "ambient_temperature");
        return heat_transfer_coefficient * area * (surface_temperature - ambient_temperature);
    }

    double newton_cooling_temperature(const double initial_temperature, const double ambient_temperature,
                                      const double cooling_constant, const double time) {
        require_finite(initial_temperature, "initial_temperature");
        require_finite(ambient_temperature, "ambient_temperature");
        require_non_negative(cooling_constant, "cooling_constant");
        require_non_negative(time, "time");
        return ambient_temperature + (initial_temperature - ambient_temperature) * std::exp(-cooling_constant * time);
    }

    double radiated_power(const double emissivity, const double area, const double temperature) {
        require_in_range(emissivity, 0.0, 1.0, "emissivity");
        require_non_negative(area, "area");
        require_non_negative(temperature, "temperature");
        return emissivity * STEFAN_BOLTZMANN * area * fourth_power(temperature);
    }

    double net_radiated_power(const double emissivity, const double area, const double temperature,
                              const double surroundings_temperature) {
        require_in_range(emissivity, 0.0, 1.0, "emissivity");
        require_non_negative(area, "area");
        require_non_negative(temperature, "temperature");
        require_non_negative(surroundings_temperature, "surroundings_temperature");
        return emissivity * STEFAN_BOLTZMANN * area *
               (fourth_power(temperature) - fourth_power(surroundings_temperature));
    }

    double internal_energy_change(const double heat_added, const double work_done_by_system) {
        require_finite(heat_added, "heat_added");
        require_finite(work_done_by_system, "work_done_by_system");
        return heat_added - work_done_by_system;
    }

    double carnot_efficiency(const double cold_temperature, const double hot_temperature) {
        require_positive(cold_temperature, "cold_temperature");
        require_positive(hot_temperature, "hot_temperature");
        if (cold_temperature > hot_temperature) {
            detail::throw_invalid("cold_temperature must not exceed hot_temperature");
        }
        return 1.0 - cold_temperature / hot_temperature;
    }

    double carnot_cop_refrigerator(const double cold_temperature, const double hot_temperature) {
        require_positive(cold_temperature, "cold_temperature");
        require_positive(hot_temperature, "hot_temperature");
        if (cold_temperature >= hot_temperature) {
            detail::throw_invalid("cold_temperature must be strictly less than hot_temperature");
        }
        return cold_temperature / (hot_temperature - cold_temperature);
    }

    double ideal_gas_pressure(const double amount, const double temperature, const double volume) {
        require_non_negative(amount, "amount");
        require_non_negative(temperature, "temperature");
        require_positive(volume, "volume");
        return amount * GAS_CONSTANT * temperature / volume;
    }

    double ideal_gas_volume(const double amount, const double temperature, const double pressure) {
        require_non_negative(amount, "amount");
        require_non_negative(temperature, "temperature");
        require_positive(pressure, "pressure");
        return amount * GAS_CONSTANT * temperature / pressure;
    }

    double linear_thermal_expansion(const double expansion_coefficient, const double initial_length,
                                    const double temperature_change) {
        require_finite(expansion_coefficient, "expansion_coefficient");
        require_non_negative(initial_length, "initial_length");
        require_finite(temperature_change, "temperature_change");
        return expansion_coefficient * initial_length * temperature_change;
    }

    double celsius_to_kelvin(const double celsius) {
        require_finite(celsius, "celsius");
        const double kelvin = celsius + CELSIUS_OFFSET;
        if (kelvin < 0.0) {
            detail::throw_invalid("celsius", "at or above absolute zero (-273.15)", celsius);
        }
        return kelvin;
    }

    double kelvin_to_celsius(const double kelvin) {
        require_non_negative(kelvin, "kelvin");
        return kelvin - CELSIUS_OFFSET;
    }

} // namespace physicsforge::thermodynamics

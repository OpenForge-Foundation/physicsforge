// src/bindings/thermodynamics.cpp - exposes physicsforge::thermodynamics
#include <vector>

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include "binding_utils.h"
#include "physicsforge/thermodynamics/heat.h"
#include "physicsforge/thermodynamics/heat_equation.h"

namespace py = pybind11;

void bind_thermodynamics(py::module_ &m) {
    namespace th = physicsforge::thermodynamics;
    using pfpy::as_span;
    using pfpy::DoubleArray;
    using pfpy::to_numpy;

    // Heat laws
    m.def("sensible_heat", py::vectorize(&th::sensible_heat), py::arg("mass"), py::arg("specific_heat"),
          py::arg("temperature_change"));
    m.def("latent_heat", py::vectorize(&th::latent_heat), py::arg("mass"), py::arg("specific_latent_heat"));
    m.def("fourier_heat_flux", py::vectorize(&th::fourier_heat_flux), py::arg("thermal_conductivity"),
          py::arg("temperature_gradient"));
    m.def("conduction_rate", py::vectorize(&th::conduction_rate), py::arg("thermal_conductivity"), py::arg("area"),
          py::arg("temperature_difference"), py::arg("thickness"));
    m.def("conduction_resistance", py::vectorize(&th::conduction_resistance), py::arg("thickness"),
          py::arg("thermal_conductivity"), py::arg("area"));
    m.def("convection_resistance", py::vectorize(&th::convection_resistance), py::arg("heat_transfer_coefficient"),
          py::arg("area"));
    m.def("heat_rate_from_resistance", py::vectorize(&th::heat_rate_from_resistance), py::arg("temperature_difference"),
          py::arg("thermal_resistance"));
    m.def("thermal_diffusivity", py::vectorize(&th::thermal_diffusivity), py::arg("thermal_conductivity"),
          py::arg("density"), py::arg("specific_heat"));
    m.def("newton_cooling_rate", py::vectorize(&th::newton_cooling_rate), py::arg("heat_transfer_coefficient"),
          py::arg("area"), py::arg("surface_temperature"), py::arg("ambient_temperature"));
    m.def("newton_cooling_temperature", py::vectorize(&th::newton_cooling_temperature), py::arg("initial_temperature"),
          py::arg("ambient_temperature"), py::arg("cooling_constant"), py::arg("time"));
    m.def("radiated_power", py::vectorize(&th::radiated_power), py::arg("emissivity"), py::arg("area"),
          py::arg("temperature"));
    m.def("net_radiated_power", py::vectorize(&th::net_radiated_power), py::arg("emissivity"), py::arg("area"),
          py::arg("temperature"), py::arg("surroundings_temperature"));
    m.def("internal_energy_change", py::vectorize(&th::internal_energy_change), py::arg("heat_added"),
          py::arg("work_done_by_system"));
    m.def("carnot_efficiency", py::vectorize(&th::carnot_efficiency), py::arg("cold_temperature"),
          py::arg("hot_temperature"));
    m.def("carnot_cop_refrigerator", py::vectorize(&th::carnot_cop_refrigerator), py::arg("cold_temperature"),
          py::arg("hot_temperature"));
    m.def("ideal_gas_pressure", py::vectorize(&th::ideal_gas_pressure), py::arg("amount"), py::arg("temperature"),
          py::arg("volume"));
    m.def("ideal_gas_volume", py::vectorize(&th::ideal_gas_volume), py::arg("amount"), py::arg("temperature"),
          py::arg("pressure"));
    m.def("linear_thermal_expansion", py::vectorize(&th::linear_thermal_expansion), py::arg("expansion_coefficient"),
          py::arg("initial_length"), py::arg("temperature_change"));
    m.def("celsius_to_kelvin", py::vectorize(&th::celsius_to_kelvin), py::arg("celsius"));
    m.def("kelvin_to_celsius", py::vectorize(&th::kelvin_to_celsius), py::arg("kelvin"));

    // Heat equation (heavy loops run with the GIL released)
    m.def("fourier_number", py::vectorize(&th::fourier_number), py::arg("diffusivity"), py::arg("time_step"),
          py::arg("grid_spacing"));
    m.def(
            "heat_fourier_coefficients",
            [](const DoubleArray &initial, const int n_terms, const double left, const double right) {
                const auto span = as_span(initial, "initial_temperatures");
                std::vector<double> result;
                {
                    py::gil_scoped_release release;
                    result = th::heat_fourier_coefficients(span, n_terms, left, right);
                }
                return to_numpy(std::move(result));
            },
            py::arg("initial_temperatures"), py::arg("n_terms"), py::arg("left_temperature"),
            py::arg("right_temperature"));
    m.def(
            "heat_fourier_solution",
            [](const DoubleArray &coefficients, const double length, const double diffusivity, const double time,
               const DoubleArray &positions, const double left, const double right) {
                const auto c = as_span(coefficients, "coefficients");
                const auto x = as_span(positions, "positions");
                std::vector<double> result;
                {
                    py::gil_scoped_release release;
                    result = th::heat_fourier_solution(c, length, diffusivity, time, x, left, right);
                }
                return to_numpy(std::move(result));
            },
            py::arg("coefficients"), py::arg("length"), py::arg("diffusivity"), py::arg("time"), py::arg("positions"),
            py::arg("left_temperature"), py::arg("right_temperature"));
    m.def(
            "heat_equation_explicit",
            [](const DoubleArray &initial, const double diffusivity, const double dx, const double dt, const int steps,
               const double left, const double right) {
                const auto span = as_span(initial, "initial_temperatures");
                std::vector<double> result;
                {
                    py::gil_scoped_release release;
                    result = th::heat_equation_explicit(span, diffusivity, dx, dt, steps, left, right);
                }
                return to_numpy(std::move(result));
            },
            py::arg("initial_temperatures"), py::arg("diffusivity"), py::arg("grid_spacing"), py::arg("time_step"),
            py::arg("steps"), py::arg("left_temperature"), py::arg("right_temperature"));
    m.def(
            "heat_equation_crank_nicolson",
            [](const DoubleArray &initial, const double diffusivity, const double dx, const double dt, const int steps,
               const double left, const double right) {
                const auto span = as_span(initial, "initial_temperatures");
                std::vector<double> result;
                {
                    py::gil_scoped_release release;
                    result = th::heat_equation_crank_nicolson(span, diffusivity, dx, dt, steps, left, right);
                }
                return to_numpy(std::move(result));
            },
            py::arg("initial_temperatures"), py::arg("diffusivity"), py::arg("grid_spacing"), py::arg("time_step"),
            py::arg("steps"), py::arg("left_temperature"), py::arg("right_temperature"));
}

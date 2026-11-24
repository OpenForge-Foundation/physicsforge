// python/bindings/electromagnetism.cpp - exposes physicsforge::electromagnetism
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include "binding_utils.h"
#include "physicsforge/electromagnetism/circuits.h"
#include "physicsforge/electromagnetism/electrostatics.h"
#include "physicsforge/electromagnetism/maxwell.h"

namespace py = pybind11;

void bind_electromagnetism(py::module_ &m) {
    namespace em = physicsforge::electromagnetism;
    using pfpy::as_span;
    using pfpy::DoubleArray;

    // Circuits
    m.def("ohms_law_voltage", py::vectorize(&em::ohms_law_voltage), py::arg("current"), py::arg("resistance"));
    m.def("ohms_law_current", py::vectorize(&em::ohms_law_current), py::arg("voltage"), py::arg("resistance"));
    m.def("ohms_law_resistance", py::vectorize(&em::ohms_law_resistance), py::arg("voltage"), py::arg("current"));
    m.def("electrical_power", py::vectorize(&em::electrical_power), py::arg("voltage"), py::arg("current"));
    m.def("power_from_current", py::vectorize(&em::power_from_current), py::arg("current"), py::arg("resistance"));
    m.def("power_from_voltage", py::vectorize(&em::power_from_voltage), py::arg("voltage"), py::arg("resistance"));
    m.def("voltage_divider", py::vectorize(&em::voltage_divider), py::arg("input_voltage"), py::arg("r1"),
          py::arg("r2"));
    m.def("current_divider", py::vectorize(&em::current_divider), py::arg("total_current"), py::arg("r_branch"),
          py::arg("r_other"));
    m.def(
            "series_resistance", [](const DoubleArray &r) { return em::series_resistance(as_span(r, "resistances")); },
            py::arg("resistances"));
    m.def(
            "parallel_resistance",
            [](const DoubleArray &r) { return em::parallel_resistance(as_span(r, "resistances")); },
            py::arg("resistances"));
    m.def(
            "series_capacitance",
            [](const DoubleArray &c) { return em::series_capacitance(as_span(c, "capacitances")); },
            py::arg("capacitances"));
    m.def(
            "parallel_capacitance",
            [](const DoubleArray &c) { return em::parallel_capacitance(as_span(c, "capacitances")); },
            py::arg("capacitances"));
    m.def("resistance_from_resistivity", py::vectorize(&em::resistance_from_resistivity), py::arg("resistivity"),
          py::arg("length"), py::arg("area"));
    m.def("rc_time_constant", py::vectorize(&em::rc_time_constant), py::arg("resistance"), py::arg("capacitance"));
    m.def("capacitor_charging_voltage", py::vectorize(&em::capacitor_charging_voltage), py::arg("source_voltage"),
          py::arg("time"), py::arg("resistance"), py::arg("capacitance"));
    m.def("capacitor_discharging_voltage", py::vectorize(&em::capacitor_discharging_voltage),
          py::arg("initial_voltage"), py::arg("time"), py::arg("resistance"), py::arg("capacitance"));

    // Electrostatics
    m.def("coulomb_force", py::vectorize(&em::coulomb_force), py::arg("charge1"), py::arg("charge2"),
          py::arg("distance"));
    m.def("coulomb_force_vector", &em::coulomb_force_vector, py::arg("charge1"), py::arg("charge2"),
          py::arg("position1"), py::arg("position2"));
    m.def("electric_field_point_charge", py::vectorize(&em::electric_field_point_charge), py::arg("charge"),
          py::arg("distance"));
    m.def("electric_field_vector", &em::electric_field_vector, py::arg("charge"), py::arg("source"), py::arg("point"));
    m.def("electric_potential_point_charge", py::vectorize(&em::electric_potential_point_charge), py::arg("charge"),
          py::arg("distance"));
    m.def("electric_potential_energy", py::vectorize(&em::electric_potential_energy), py::arg("charge1"),
          py::arg("charge2"), py::arg("distance"));
    m.def("parallel_plate_capacitance", py::vectorize(&em::parallel_plate_capacitance),
          py::arg("relative_permittivity"), py::arg("area"), py::arg("separation"));
    m.def("capacitor_energy", py::vectorize(&em::capacitor_energy), py::arg("capacitance"), py::arg("voltage"));

    // Maxwell
    m.def("gauss_electric_flux", py::vectorize(&em::gauss_electric_flux), py::arg("enclosed_charge"));
    m.def("gauss_enclosed_charge", py::vectorize(&em::gauss_enclosed_charge), py::arg("electric_flux"));
    m.def("magnetic_flux", py::vectorize(&em::magnetic_flux), py::arg("field_strength"), py::arg("area"),
          py::arg("angle"));
    m.def("faraday_emf", py::vectorize(&em::faraday_emf), py::arg("flux_rate_of_change"), py::arg("turns"));
    m.def("displacement_current", py::vectorize(&em::displacement_current), py::arg("electric_flux_rate_of_change"));
    m.def("ampere_maxwell_circulation", py::vectorize(&em::ampere_maxwell_circulation), py::arg("enclosed_current"),
          py::arg("electric_flux_rate_of_change"));
    m.def("magnetic_field_straight_wire", py::vectorize(&em::magnetic_field_straight_wire), py::arg("current"),
          py::arg("distance"));
    m.def("magnetic_field_solenoid", py::vectorize(&em::magnetic_field_solenoid), py::arg("turns_per_length"),
          py::arg("current"));
    m.def("magnetic_field_loop_center", py::vectorize(&em::magnetic_field_loop_center), py::arg("current"),
          py::arg("radius"));
    m.def("lorentz_force", &em::lorentz_force, py::arg("charge"), py::arg("electric_field"), py::arg("velocity"),
          py::arg("magnetic_field"));
    m.def("em_wave_speed", py::vectorize(&em::em_wave_speed), py::arg("relative_permittivity"),
          py::arg("relative_permeability"));
    m.def("speed_of_light_from_vacuum_constants", &em::speed_of_light_from_vacuum_constants);
    m.def("em_wave_intensity", py::vectorize(&em::em_wave_intensity), py::arg("electric_field_amplitude"));
}

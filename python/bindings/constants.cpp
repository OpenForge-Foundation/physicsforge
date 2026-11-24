// python/bindings/constants.cpp - exposes include/physicsforge/constants.h
#include <pybind11/pybind11.h>

#include "physicsforge/constants.h"

namespace py = pybind11;

void bind_constants(py::module_ &m) {
    namespace c = physicsforge::constants;
    py::module_ sub = m.def_submodule("constants", "Physical and mathematical constants (SI, CODATA 2022)");

    sub.attr("PI") = c::PI;
    sub.attr("TAU") = c::TAU;
    sub.attr("EULER") = c::EULER;

    sub.attr("SPEED_OF_LIGHT") = c::SPEED_OF_LIGHT;
    sub.attr("PLANCK") = c::PLANCK;
    sub.attr("REDUCED_PLANCK") = c::REDUCED_PLANCK;
    sub.attr("ELEMENTARY_CHARGE") = c::ELEMENTARY_CHARGE;
    sub.attr("BOLTZMANN") = c::BOLTZMANN;
    sub.attr("AVOGADRO") = c::AVOGADRO;

    sub.attr("GAS_CONSTANT") = c::GAS_CONSTANT;
    sub.attr("ELECTRON_VOLT") = c::ELECTRON_VOLT;
    sub.attr("STEFAN_BOLTZMANN") = c::STEFAN_BOLTZMANN;
    sub.attr("WIEN_DISPLACEMENT") = c::WIEN_DISPLACEMENT;
    sub.attr("STANDARD_GRAVITY") = c::STANDARD_GRAVITY;
    sub.attr("STANDARD_ATMOSPHERE") = c::STANDARD_ATMOSPHERE;

    sub.attr("GRAVITATIONAL_CONSTANT") = c::GRAVITATIONAL_CONSTANT;
    sub.attr("VACUUM_PERMITTIVITY") = c::VACUUM_PERMITTIVITY;
    sub.attr("VACUUM_PERMEABILITY") = c::VACUUM_PERMEABILITY;
    sub.attr("COULOMB_CONSTANT") = c::COULOMB_CONSTANT;
    sub.attr("ELECTRON_MASS") = c::ELECTRON_MASS;
    sub.attr("PROTON_MASS") = c::PROTON_MASS;
    sub.attr("NEUTRON_MASS") = c::NEUTRON_MASS;
}

// src/bindings/constants.cpp - exposes include/physicsforge/constants.h as module attributes
#include <pybind11/pybind11.h>

#include "physicsforge/constants.h"

namespace py = pybind11;

void bind_constants(py::module_ &m) {
    namespace c = physicsforge::constants;

    m.attr("PI") = c::PI;
    m.attr("TAU") = c::TAU;
    m.attr("EULER") = c::EULER;

    m.attr("SPEED_OF_LIGHT") = c::SPEED_OF_LIGHT;
    m.attr("PLANCK") = c::PLANCK;
    m.attr("REDUCED_PLANCK") = c::REDUCED_PLANCK;
    m.attr("ELEMENTARY_CHARGE") = c::ELEMENTARY_CHARGE;
    m.attr("BOLTZMANN") = c::BOLTZMANN;
    m.attr("AVOGADRO") = c::AVOGADRO;

    m.attr("GAS_CONSTANT") = c::GAS_CONSTANT;
    m.attr("ELECTRON_VOLT") = c::ELECTRON_VOLT;
    m.attr("STEFAN_BOLTZMANN") = c::STEFAN_BOLTZMANN;
    m.attr("WIEN_DISPLACEMENT") = c::WIEN_DISPLACEMENT;
    m.attr("STANDARD_GRAVITY") = c::STANDARD_GRAVITY;
    m.attr("STANDARD_ATMOSPHERE") = c::STANDARD_ATMOSPHERE;

    m.attr("GRAVITATIONAL_CONSTANT") = c::GRAVITATIONAL_CONSTANT;
    m.attr("VACUUM_PERMITTIVITY") = c::VACUUM_PERMITTIVITY;
    m.attr("VACUUM_PERMEABILITY") = c::VACUUM_PERMEABILITY;
    m.attr("COULOMB_CONSTANT") = c::COULOMB_CONSTANT;
    m.attr("ELECTRON_MASS") = c::ELECTRON_MASS;
    m.attr("PROTON_MASS") = c::PROTON_MASS;
    m.attr("NEUTRON_MASS") = c::NEUTRON_MASS;
}

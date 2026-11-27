// python/bindings/relativity.cpp - exposes physicsforge::relativity
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include "physicsforge/relativity.h"

namespace py = pybind11;

void bind_relativity(py::module_ &m) {
    namespace rel = physicsforge::relativity;

    m.def("lorentz_factor", py::vectorize(&rel::lorentz_factor), py::arg("velocity"));
    m.def("time_dilation", py::vectorize(&rel::time_dilation), py::arg("proper_time"), py::arg("velocity"));
    m.def("length_contraction", py::vectorize(&rel::length_contraction), py::arg("proper_length"), py::arg("velocity"));
    m.def("rest_energy", py::vectorize(&rel::rest_energy), py::arg("mass"));
    m.def("relativistic_momentum", py::vectorize(&rel::relativistic_momentum), py::arg("mass"), py::arg("velocity"));
    m.def("relativistic_kinetic_energy", py::vectorize(&rel::relativistic_kinetic_energy), py::arg("mass"),
          py::arg("velocity"));
}

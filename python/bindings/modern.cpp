// python/bindings/modern.cpp - exposes physicsforge::modern
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include "physicsforge/modern/quantum.h"
#include "physicsforge/modern/relativity.h"

namespace py = pybind11;

void bind_modern(py::module_ &m) {
    namespace mp = physicsforge::modern;

    // Quantum
    m.def("photon_energy_from_frequency", py::vectorize(&mp::photon_energy_from_frequency), py::arg("frequency"));
    m.def("photon_energy_from_wavelength", py::vectorize(&mp::photon_energy_from_wavelength), py::arg("wavelength"));
    m.def("photon_momentum", py::vectorize(&mp::photon_momentum), py::arg("wavelength"));
    m.def("wavelength_from_frequency", py::vectorize(&mp::wavelength_from_frequency), py::arg("frequency"));
    m.def("de_broglie_wavelength", py::vectorize(&mp::de_broglie_wavelength), py::arg("mass"), py::arg("velocity"));
    m.def("photoelectric_kinetic_energy", py::vectorize(&mp::photoelectric_kinetic_energy), py::arg("frequency"),
          py::arg("work_function"));
    m.def("planck_spectral_radiance", py::vectorize(&mp::planck_spectral_radiance), py::arg("wavelength"),
          py::arg("temperature"));
    m.def("wien_peak_wavelength", py::vectorize(&mp::wien_peak_wavelength), py::arg("temperature"));
    m.def("bohr_energy_level", py::vectorize(&mp::bohr_energy_level), py::arg("principal_quantum_number"));
    m.def("heisenberg_min_momentum_uncertainty", py::vectorize(&mp::heisenberg_min_momentum_uncertainty),
          py::arg("position_uncertainty"));

    // Relativity
    m.def("lorentz_factor", py::vectorize(&mp::lorentz_factor), py::arg("velocity"));
    m.def("time_dilation", py::vectorize(&mp::time_dilation), py::arg("proper_time"), py::arg("velocity"));
    m.def("length_contraction", py::vectorize(&mp::length_contraction), py::arg("proper_length"), py::arg("velocity"));
    m.def("rest_energy", py::vectorize(&mp::rest_energy), py::arg("mass"));
    m.def("relativistic_momentum", py::vectorize(&mp::relativistic_momentum), py::arg("mass"), py::arg("velocity"));
    m.def("relativistic_kinetic_energy", py::vectorize(&mp::relativistic_kinetic_energy), py::arg("mass"),
          py::arg("velocity"));
}

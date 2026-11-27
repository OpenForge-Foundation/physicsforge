// python/bindings/quantum.cpp - exposes physicsforge::quantum
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include "physicsforge/quantum.h"

namespace py = pybind11;

void bind_quantum(py::module_ &m) {
    namespace q = physicsforge::quantum;

    m.def("photon_energy_from_frequency", py::vectorize(&q::photon_energy_from_frequency), py::arg("frequency"));
    m.def("photon_energy_from_wavelength", py::vectorize(&q::photon_energy_from_wavelength), py::arg("wavelength"));
    m.def("photon_momentum", py::vectorize(&q::photon_momentum), py::arg("wavelength"));
    m.def("wavelength_from_frequency", py::vectorize(&q::wavelength_from_frequency), py::arg("frequency"));
    m.def("de_broglie_wavelength", py::vectorize(&q::de_broglie_wavelength), py::arg("mass"), py::arg("velocity"));
    m.def("photoelectric_kinetic_energy", py::vectorize(&q::photoelectric_kinetic_energy), py::arg("frequency"),
          py::arg("work_function"));
    m.def("planck_spectral_radiance", py::vectorize(&q::planck_spectral_radiance), py::arg("wavelength"),
          py::arg("temperature"));
    m.def("wien_peak_wavelength", py::vectorize(&q::wien_peak_wavelength), py::arg("temperature"));
    m.def("bohr_energy_level", py::vectorize(&q::bohr_energy_level), py::arg("principal_quantum_number"));
    m.def("heisenberg_min_momentum_uncertainty", py::vectorize(&q::heisenberg_min_momentum_uncertainty),
          py::arg("position_uncertainty"));
}

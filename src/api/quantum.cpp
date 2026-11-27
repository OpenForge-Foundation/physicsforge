// src/api/quantum.cpp - C API shims for physicsforge::quantum
#include "physicsforge/quantum.h"
#include "physicsforge/api.h"
#include "physicsforge/validation.h"

namespace q = physicsforge::quantum;
using physicsforge::detail::guarded;

extern "C" double pf_quantum_photon_energy_from_frequency(double f) {
    return guarded([&] { return q::photon_energy_from_frequency(f); });
}
extern "C" double pf_quantum_planck_spectral_radiance(double wavelength, double temperature) {
    return guarded([&] { return q::planck_spectral_radiance(wavelength, temperature); });
}

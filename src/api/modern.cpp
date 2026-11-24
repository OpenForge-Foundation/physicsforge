// src/api/modern.cpp - C API shims for physicsforge::modern
#include "api_guard.h"
#include "physicsforge/api.h"
#include "physicsforge/modern/quantum.h"
#include "physicsforge/modern/relativity.h"

namespace mp = physicsforge::modern;
using physicsforge::api::guarded;

extern "C" double pf_modern_photon_energy_from_frequency(double f) {
    return guarded([&] { return mp::photon_energy_from_frequency(f); });
}
extern "C" double pf_modern_planck_spectral_radiance(double wavelength, double temperature) {
    return guarded([&] { return mp::planck_spectral_radiance(wavelength, temperature); });
}
extern "C" double pf_modern_lorentz_factor(double v) {
    return guarded([&] { return mp::lorentz_factor(v); });
}
extern "C" double pf_modern_rest_energy(double mass) {
    return guarded([&] { return mp::rest_energy(mass); });
}

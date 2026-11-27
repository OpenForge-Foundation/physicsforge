// src/api/thermodynamics.cpp - C API shims for physicsforge::thermodynamics
#include "physicsforge/api.h"
#include "physicsforge/thermodynamics/heat.h"
#include "physicsforge/validation.h"

namespace th = physicsforge::thermodynamics;
using physicsforge::detail::guarded;

extern "C" double pf_thermo_conduction_rate(double k, double area, double dt, double thickness) {
    return guarded([&] { return th::conduction_rate(k, area, dt, thickness); });
}
extern "C" double pf_thermo_newton_cooling_temperature(double t0, double t_env, double k, double time) {
    return guarded([&] { return th::newton_cooling_temperature(t0, t_env, k, time); });
}
extern "C" double pf_thermo_radiated_power(double emissivity, double area, double temperature) {
    return guarded([&] { return th::radiated_power(emissivity, area, temperature); });
}
extern "C" double pf_thermo_carnot_efficiency(double t_cold, double t_hot) {
    return guarded([&] { return th::carnot_efficiency(t_cold, t_hot); });
}

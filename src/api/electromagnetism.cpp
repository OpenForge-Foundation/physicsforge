// src/api/electromagnetism.cpp - C API shims for physicsforge::electromagnetism
#include "physicsforge/api.h"
#include "physicsforge/electromagnetism/circuits.h"
#include "physicsforge/electromagnetism/electrostatics.h"
#include "physicsforge/electromagnetism/maxwell.h"
#include "physicsforge/validation.h"

namespace em = physicsforge::electromagnetism;
using physicsforge::detail::guarded;

extern "C" double pf_em_ohms_law_voltage(double i, double r) {
    return guarded([&] { return em::ohms_law_voltage(i, r); });
}
extern "C" double pf_em_ohms_law_current(double v, double r) {
    return guarded([&] { return em::ohms_law_current(v, r); });
}
extern "C" double pf_em_voltage_divider(double vin, double r1, double r2) {
    return guarded([&] { return em::voltage_divider(vin, r1, r2); });
}
extern "C" double pf_em_coulomb_force(double q1, double q2, double r) {
    return guarded([&] { return em::coulomb_force(q1, q2, r); });
}
extern "C" double pf_em_gauss_electric_flux(double q) {
    return guarded([&] { return em::gauss_electric_flux(q); });
}
extern "C" double pf_em_faraday_emf(double dphi_dt, double turns) {
    return guarded([&] { return em::faraday_emf(dphi_dt, turns); });
}
extern "C" double pf_em_magnetic_field_straight_wire(double i, double r) {
    return guarded([&] { return em::magnetic_field_straight_wire(i, r); });
}

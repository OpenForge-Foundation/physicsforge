// src/api/relativity.cpp - C API shims for physicsforge::relativity
#include "physicsforge/relativity.h"
#include "physicsforge/api.h"
#include "physicsforge/validation.h"

namespace rel = physicsforge::relativity;
using physicsforge::detail::guarded;

extern "C" double pf_relativity_lorentz_factor(double v) {
    return guarded([&] { return rel::lorentz_factor(v); });
}
extern "C" double pf_relativity_rest_energy(double mass) {
    return guarded([&] { return rel::rest_energy(mass); });
}

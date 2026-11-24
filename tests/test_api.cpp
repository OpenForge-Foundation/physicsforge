// tests/test_api.cpp - the C ABI returns values on success and NaN on invalid input.
#include <cmath>
#include <cstring>

#include "physicsforge/api.h"
#include "physicsforge/version.h"
#include "test_framework.h"

void test_api() {
    PF_CHECK(pf_version_major() == PHYSICSFORGE_VERSION_MAJOR);
    PF_CHECK(pf_version_minor() == PHYSICSFORGE_VERSION_MINOR);
    PF_CHECK(pf_version_patch() == PHYSICSFORGE_VERSION_PATCH);
    PF_CHECK(std::strcmp(pf_version_string(), PHYSICSFORGE_VERSION_STRING) == 0);

    PF_CHECK_CLOSE(pf_em_ohms_law_voltage(2.0, 5.0), 10.0);
    PF_CHECK_CLOSE(pf_em_voltage_divider(12.0, 1000.0, 2000.0), 8.0);
    PF_CHECK_CLOSE(pf_mechanics_net_force(2.0, 3.0), 6.0);
    PF_CHECK_CLOSE(pf_modern_lorentz_factor(0.0), 1.0);

    // Invalid inputs surface as NaN, never as an exception crossing the C boundary.
    PF_CHECK(std::isnan(pf_em_ohms_law_current(1.0, 0.0)));
    PF_CHECK(std::isnan(pf_mechanics_gravitational_force(1.0, 1.0, -1.0)));
    PF_CHECK(std::isnan(pf_modern_lorentz_factor(3e8)));
    PF_CHECK(std::isnan(pf_thermo_carnot_efficiency(600.0, 300.0)));
}

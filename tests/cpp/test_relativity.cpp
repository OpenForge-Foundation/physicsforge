// tests/cpp/test_relativity.cpp
#include <stdexcept>

#include "physicsforge/constants.h"
#include "physicsforge/relativity.h"
#include "test_framework.h"

void test_relativity() {
    using namespace physicsforge::relativity;
    namespace c = physicsforge::constants;

    PF_CHECK_CLOSE(lorentz_factor(0.0), 1.0);
    PF_CHECK_CLOSE(lorentz_factor(0.6 * c::SPEED_OF_LIGHT), 1.25);
    PF_CHECK_CLOSE(time_dilation(1.0, 0.6 * c::SPEED_OF_LIGHT), 1.25);
    PF_CHECK_CLOSE(length_contraction(1.0, 0.6 * c::SPEED_OF_LIGHT), 0.8);
    PF_CHECK_THROWS(lorentz_factor(c::SPEED_OF_LIGHT), std::invalid_argument);
    PF_CHECK_CLOSE(rest_energy(1.0), c::SPEED_OF_LIGHT * c::SPEED_OF_LIGHT);
    PF_CHECK_CLOSE(relativistic_kinetic_energy(1.0, 0.6 * c::SPEED_OF_LIGHT), 0.25 * rest_energy(1.0));
    // Low-speed limit matches m v^2 / 2 without cancellation error.
    PF_CHECK_CLOSE_TOL(relativistic_kinetic_energy(1.0, 1.0), 0.5, 1e-12, 0.0);
    PF_CHECK_CLOSE(relativistic_momentum(1.0, 0.6 * c::SPEED_OF_LIGHT), 1.25 * 0.6 * c::SPEED_OF_LIGHT);
}

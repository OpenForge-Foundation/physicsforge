// tests/test_constants.cpp
#include <cmath>

#include "physicsforge/constants.h"
#include "test_framework.h"

void test_constants() {
    namespace c = physicsforge::constants;

    // Internal consistency of the CODATA set.
    // c^2 eps_0 mu_0 = 1 to within the measurement uncertainty of eps_0 / mu_0 (~1.6e-10).
    PF_CHECK_CLOSE_TOL(c::SPEED_OF_LIGHT * c::SPEED_OF_LIGHT * c::VACUUM_PERMITTIVITY * c::VACUUM_PERMEABILITY, 1.0,
                       1e-9, 0.0);
    PF_CHECK_CLOSE(c::REDUCED_PLANCK * c::TAU, c::PLANCK);
    PF_CHECK_CLOSE(c::GAS_CONSTANT, 8.31446261815324); // exact since the 2019 SI redefinition
    // sigma = 2 pi^5 k^4 / (15 h^3 c^2)
    const double sigma = 2.0 * std::pow(c::PI, 5) * std::pow(c::BOLTZMANN, 4) /
                         (15.0 * std::pow(c::PLANCK, 3) * c::SPEED_OF_LIGHT * c::SPEED_OF_LIGHT);
    PF_CHECK_CLOSE_TOL(c::STEFAN_BOLTZMANN, sigma, 1e-9, 0.0);
    PF_CHECK_CLOSE(c::EULER, std::exp(1.0));
    PF_CHECK_CLOSE(c::PI, std::acos(-1.0));
    PF_CHECK_CLOSE_TOL(c::COULOMB_CONSTANT, 8.9875517862e9, 1e-10, 0.0);
}

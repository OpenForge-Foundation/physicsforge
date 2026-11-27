// tests/cpp/test_mechanics.cpp
#include <cmath>
#include <stdexcept>

#include "physicsforge/constants.h"
#include "physicsforge/mechanics/dynamics.h"
#include "physicsforge/mechanics/gravitation.h"
#include "physicsforge/mechanics/kinematics.h"
#include "physicsforge/mechanics/rocketry.h"
#include "physicsforge/mechanics/rotation.h"
#include "test_framework.h"

void test_mechanics() {
    using namespace physicsforge::mechanics;
    namespace c = physicsforge::constants;
    using physicsforge::math::Vec3;

    // --- Kinematics ---
    PF_CHECK_CLOSE(displacement(2.0, 3.0, 4.0), 2.0 * 4.0 + 0.5 * 3.0 * 16.0);
    PF_CHECK_CLOSE(final_velocity(2.0, 3.0, 4.0), 14.0);
    PF_CHECK_CLOSE(final_velocity_from_displacement(0.0, 9.80665, 10.0), std::sqrt(2.0 * 9.80665 * 10.0));
    PF_CHECK_THROWS(final_velocity_from_displacement(1.0, -1.0, 10.0), std::invalid_argument);
    PF_CHECK_THROWS(displacement(1.0, 1.0, -1.0), std::invalid_argument);
    PF_CHECK_CLOSE(average_velocity(100.0, 20.0), 5.0);
    // 45 degrees maximises range: R = v0^2 / g.
    PF_CHECK_CLOSE(projectile_range(10.0, c::PI / 4.0, 10.0), 10.0);
    PF_CHECK_CLOSE(projectile_time_of_flight(10.0, c::PI / 2.0, 10.0), 2.0);
    PF_CHECK_CLOSE(projectile_max_height(10.0, c::PI / 2.0, 10.0), 5.0);

    // --- Dynamics ---
    PF_CHECK_CLOSE(net_force(2.0, 3.0), 6.0);
    PF_CHECK_CLOSE(acceleration_from_force(6.0, 2.0), 3.0);
    PF_CHECK_THROWS(net_force(-1.0, 3.0), std::invalid_argument);
    PF_CHECK_CLOSE(kinetic_energy(2.0, 3.0), 9.0);
    PF_CHECK_CLOSE(momentum(2.0, -3.0), -6.0);
    PF_CHECK_CLOSE_TOL(work(10.0, 2.0, c::PI / 2.0), 0.0, 0.0, 1e-14);
    PF_CHECK_CLOSE(spring_force(100.0, 0.1), -10.0);
    PF_CHECK_CLOSE(spring_potential_energy(100.0, 0.1), 0.5);
    PF_CHECK_CLOSE(power(100.0, 4.0), 25.0);

    // --- Gravitation: Earth ---
    constexpr double earth_mass = 5.9722e24;
    constexpr double earth_radius = 6.371e6;
    PF_CHECK_CLOSE_TOL(gravitational_field_strength(earth_mass, earth_radius), 9.82, 1e-3, 0.0);
    PF_CHECK_CLOSE_TOL(escape_velocity(earth_mass, earth_radius), 11186.0, 1e-3, 0.0);
    PF_CHECK_CLOSE(escape_velocity(earth_mass, earth_radius),
                   std::sqrt(2.0) * orbital_velocity(earth_mass, earth_radius));
    // Moon's sidereal period ~27.3 days at a = 384,400 km (Earth+Moon mass).
    PF_CHECK_CLOSE_TOL(orbital_period(3.844e8, earth_mass + 7.342e22) / 86400.0, 27.3, 5e-3, 0.0);
    const Vec3 f12 = gravitational_force_vector(1.0, 1.0, Vec3{1, 0, 0}, Vec3{0, 0, 0});
    PF_CHECK_CLOSE(f12.x, -c::GRAVITATIONAL_CONSTANT); // attractive: points toward body 2
    PF_CHECK_THROWS(gravitational_force(1.0, 1.0, 0.0), std::invalid_argument);

    // --- Rotation ---
    PF_CHECK_CLOSE(torque(2.0, 10.0, c::PI / 2.0), 20.0);
    PF_CHECK((torque_vector(Vec3{1, 0, 0}, Vec3{0, 1, 0}) == Vec3{0, 0, 1}));
    PF_CHECK_CLOSE(moment_of_inertia_solid_sphere(5.0, 2.0), 8.0);
    PF_CHECK_CLOSE(moment_of_inertia_rod_end(3.0, 2.0), 4.0);
    // Parallel-axis theorem links the two rod formulas: I_end = I_cm + m (L/2)^2.
    PF_CHECK_CLOSE(parallel_axis(moment_of_inertia_rod_center(3.0, 2.0), 3.0, 1.0),
                   moment_of_inertia_rod_end(3.0, 2.0));
    PF_CHECK_CLOSE(rotational_kinetic_energy(2.0, 3.0), 9.0);
    PF_CHECK_CLOSE(centripetal_force(2.0, 3.0, 1.5), 12.0);

    // --- Rocketry ---
    PF_CHECK_CLOSE(rocket_delta_v(3000.0, c::EULER * 1000.0, 1000.0), 3000.0);
    PF_CHECK_CLOSE(rocket_delta_v(3000.0, 1000.0, 1000.0), 0.0);
    PF_CHECK_CLOSE(rocket_mass_ratio(3000.0, 3000.0), c::EULER);
    PF_CHECK_CLOSE(rocket_propellant_mass(3000.0, 3000.0, 1000.0), (c::EULER - 1.0) * 1000.0);
    // Round trip: the propellant needed for dv yields exactly dv.
    const double mp = rocket_propellant_mass(9400.0, 4400.0, 25000.0);
    PF_CHECK_CLOSE(rocket_delta_v(4400.0, 25000.0 + mp, 25000.0), 9400.0);
    PF_CHECK_CLOSE(exhaust_velocity_from_isp(300.0, c::STANDARD_GRAVITY), 2941.995);
    PF_CHECK_CLOSE(specific_impulse(2941.995, c::STANDARD_GRAVITY), 300.0);
    PF_CHECK_CLOSE(rocket_thrust(100.0, 3000.0, 0.0, 0.0, 0.0), 300000.0);
    PF_CHECK_THROWS(rocket_delta_v(3000.0, 500.0, 1000.0), std::invalid_argument);
    PF_CHECK_CLOSE(burn_time(1000.0, 50.0), 20.0);
}

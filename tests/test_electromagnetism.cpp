// tests/test_electromagnetism.cpp
#include <stdexcept>
#include <vector>

#include "physicsforge/constants.h"
#include "physicsforge/electromagnetism/circuits.h"
#include "physicsforge/electromagnetism/electrostatics.h"
#include "physicsforge/electromagnetism/maxwell.h"
#include "test_framework.h"

void test_electromagnetism() {
    using namespace physicsforge::electromagnetism;
    namespace c = physicsforge::constants;
    using physicsforge::math::Vec3;

    // --- Circuits ---
    PF_CHECK_CLOSE(ohms_law_voltage(2.0, 5.0), 10.0);
    PF_CHECK_CLOSE(ohms_law_current(10.0, 5.0), 2.0);
    PF_CHECK_CLOSE(ohms_law_resistance(10.0, 2.0), 5.0);
    PF_CHECK_THROWS(ohms_law_current(10.0, 0.0), std::invalid_argument);
    PF_CHECK_THROWS(ohms_law_resistance(10.0, 0.0), std::invalid_argument);
    PF_CHECK_CLOSE(power_from_current(2.0, 5.0), electrical_power(10.0, 2.0));
    PF_CHECK_CLOSE(power_from_voltage(10.0, 5.0), 20.0);
    PF_CHECK_CLOSE(voltage_divider(12.0, 1000.0, 2000.0), 8.0);
    PF_CHECK_CLOSE(voltage_divider(12.0, 0.0, 2000.0), 12.0);
    PF_CHECK_THROWS(voltage_divider(12.0, 0.0, 0.0), std::invalid_argument);
    PF_CHECK_CLOSE(current_divider(3.0, 10.0, 20.0), 2.0);
    const std::vector<double> rs{100.0, 200.0, 300.0};
    PF_CHECK_CLOSE(series_resistance(rs), 600.0);
    PF_CHECK_CLOSE(parallel_resistance(std::vector<double>{100.0, 100.0}), 50.0);
    PF_CHECK_CLOSE(series_capacitance(std::vector<double>{2e-6, 2e-6}), 1e-6);
    PF_CHECK_CLOSE(parallel_capacitance(std::vector<double>{1e-6, 2e-6}), 3e-6);
    PF_CHECK_THROWS(parallel_resistance(std::vector<double>{}), std::invalid_argument);
    PF_CHECK_CLOSE(rc_time_constant(1000.0, 1e-3), 1.0);
    PF_CHECK_CLOSE(capacitor_charging_voltage(5.0, 1.0, 1000.0, 1e-3), 5.0 * (1.0 - std::exp(-1.0)));
    PF_CHECK_CLOSE(capacitor_discharging_voltage(5.0, 1.0, 1000.0, 1e-3), 5.0 * std::exp(-1.0));
    PF_CHECK_CLOSE(resistance_from_resistivity(1.68e-8, 100.0, 1e-6), 1.68);

    // --- Electrostatics ---
    PF_CHECK_CLOSE(coulomb_force(1.0, 1.0, 1.0), c::COULOMB_CONSTANT);
    PF_CHECK(coulomb_force(1.0, -1.0, 1.0) < 0.0); // attractive
    const Vec3 f = coulomb_force_vector(1.0, 1.0, Vec3{2, 0, 0}, Vec3{0, 0, 0});
    PF_CHECK_CLOSE(f.x, c::COULOMB_CONSTANT / 4.0); // repulsive: points away from charge 2
    const Vec3 e = electric_field_vector(1e-9, Vec3{0, 0, 0}, Vec3{0, 3, 4});
    PF_CHECK_CLOSE(physicsforge::math::norm(e), electric_field_point_charge(1e-9, 5.0));
    PF_CHECK_CLOSE(electric_potential_point_charge(1e-9, 1.0) * 1e-9, electric_potential_energy(1e-9, 1e-9, 1.0));
    PF_CHECK_CLOSE(parallel_plate_capacitance(1.0, 1.0, 1.0), c::VACUUM_PERMITTIVITY);
    PF_CHECK_CLOSE(capacitor_energy(2e-6, 10.0), 1e-4);

    // --- Maxwell ---
    // Gauss's law agrees with Coulomb: flux of a point charge's field through a sphere.
    constexpr double q = 1e-9;
    constexpr double r = 0.5;
    PF_CHECK_CLOSE(electric_field_point_charge(q, r) * 4.0 * c::PI * r * r, gauss_electric_flux(q));
    PF_CHECK_CLOSE(gauss_enclosed_charge(gauss_electric_flux(q)), q);
    PF_CHECK_CLOSE(magnetic_flux(2.0, 3.0, 0.0), 6.0);
    PF_CHECK_CLOSE(faraday_emf(0.5, 100.0), -50.0);
    PF_CHECK_THROWS(faraday_emf(0.5, 0.0), std::invalid_argument);
    PF_CHECK_CLOSE(ampere_maxwell_circulation(1.0, 0.0), c::VACUUM_PERMEABILITY);
    PF_CHECK_CLOSE(ampere_maxwell_circulation(0.0, 1.0 / c::VACUUM_PERMITTIVITY), c::VACUUM_PERMEABILITY);
    // Ampere's law: circulation around a wire equals B * 2 pi r.
    PF_CHECK_CLOSE(magnetic_field_straight_wire(3.0, 0.2) * c::TAU * 0.2, ampere_maxwell_circulation(3.0, 0.0));
    PF_CHECK_CLOSE(magnetic_field_solenoid(1000.0, 2.0), c::VACUUM_PERMEABILITY * 2000.0);
    PF_CHECK_CLOSE(magnetic_field_loop_center(1.0, 0.5), c::VACUUM_PERMEABILITY);
    PF_CHECK_CLOSE_TOL(speed_of_light_from_vacuum_constants(), c::SPEED_OF_LIGHT, 1e-9, 0.0);
    PF_CHECK_CLOSE(em_wave_speed(4.0, 1.0), c::SPEED_OF_LIGHT / 2.0);
    // Magnetic force on +q moving along x in B along y points along z.
    const Vec3 lf = lorentz_force(2.0, Vec3{}, Vec3{3, 0, 0}, Vec3{0, 5, 0});
    PF_CHECK((lf == Vec3{0, 0, 30}));
    PF_CHECK_CLOSE(em_wave_intensity(1.0), 0.5 * c::SPEED_OF_LIGHT * c::VACUUM_PERMITTIVITY);
}

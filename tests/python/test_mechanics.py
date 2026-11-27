import math

import numpy as np
import pytest

from physicsforge import mechanics as m
from physicsforge.constants import GRAVITATIONAL_CONSTANT, STANDARD_GRAVITY

EARTH_MASS = 5.9722e24
EARTH_RADIUS = 6.371e6


class TestKinematics:
    def test_suvat_consistency(self):
        v0, a, t = 3.0, 2.0, 5.0
        s = m.displacement(v0, a, t)
        v = m.final_velocity(v0, a, t)
        assert m.final_velocity_from_displacement(v0, a, s) == pytest.approx(v)
        assert m.average_velocity(s, t) == pytest.approx((v0 + v) / 2)

    def test_projectile_45_degrees_maximises_range(self):
        angles = np.linspace(0.01, math.pi / 2 - 0.01, 1001)
        ranges = m.projectile_range(20.0, angles)
        assert angles[np.argmax(ranges)] == pytest.approx(math.pi / 4, abs=2e-3)
        assert m.projectile_range(20.0, math.pi / 4) == pytest.approx(400.0 / STANDARD_GRAVITY)

    def test_projectile_flight_relations(self):
        v0, theta = 30.0, 0.6
        t = m.projectile_time_of_flight(v0, theta)
        assert m.projectile_range(v0, theta) == pytest.approx(v0 * math.cos(theta) * t)
        assert m.projectile_max_height(v0, theta) == pytest.approx(
            v0 * math.sin(theta) * t / 2 - 0.5 * STANDARD_GRAVITY * (t / 2) ** 2
        )

    @pytest.mark.parametrize("bad_time", [-1.0, math.inf, math.nan])
    def test_rejects_invalid_time(self, bad_time):
        with pytest.raises(ValueError):
            m.displacement(1.0, 1.0, bad_time)

    def test_unreachable_displacement(self):
        with pytest.raises(ValueError, match="never reached"):
            m.final_velocity_from_displacement(1.0, -1.0, 10.0)


class TestDynamics:
    def test_newtons_second_law_round_trip(self):
        assert m.acceleration_from_force(m.net_force(4.0, 2.5), 4.0) == pytest.approx(2.5)

    def test_work_energy_theorem(self):
        mass, v0, force, d = 2.0, 1.0, 10.0, 3.0
        a = m.acceleration_from_force(force, mass)
        v = m.final_velocity_from_displacement(v0, a, d)
        assert m.kinetic_energy(mass, v) - m.kinetic_energy(mass, v0) == pytest.approx(m.work(force, d))

    def test_impulse_momentum(self):
        assert m.impulse(20.0, 0.5) == pytest.approx(m.momentum(2.0, 5.0))

    def test_spring_energy_is_work_against_spring(self):
        assert m.spring_potential_energy(50.0, 0.2) == pytest.approx(-0.5 * m.spring_force(50.0, 0.2) * 0.2)

    def test_weight_defaults_to_standard_gravity(self):
        assert m.weight(1.0) == STANDARD_GRAVITY

    def test_negative_mass_rejected(self):
        with pytest.raises(ValueError, match="mass"):
            m.kinetic_energy(-1.0, 1.0)


class TestGravitation:
    def test_surface_gravity_and_escape_velocity(self):
        assert m.gravitational_field_strength(EARTH_MASS, EARTH_RADIUS) == pytest.approx(9.82, rel=1e-3)
        assert m.escape_velocity(EARTH_MASS, EARTH_RADIUS) == pytest.approx(
            math.sqrt(2) * m.orbital_velocity(EARTH_MASS, EARTH_RADIUS)
        )

    def test_kepler_third_law_matches_circular_speed(self):
        r = 7.0e6
        period = m.orbital_period(r, EARTH_MASS)
        assert 2 * math.pi * r / period == pytest.approx(m.orbital_velocity(EARTH_MASS, r))

    def test_force_vector_is_newtons_third_law_pair(self):
        f12 = m.gravitational_force_vector(3.0, 5.0, (1, 2, 3), (4, 6, 3))
        f21 = m.gravitational_force_vector(5.0, 3.0, (4, 6, 3), (1, 2, 3))
        assert f12 == -f21
        assert f12.norm() == pytest.approx(m.gravitational_force(3.0, 5.0, 5.0))

    def test_potential_energy_relation(self):
        r = 2.0
        assert m.gravitational_potential_energy(1.0, 1.0, r) == pytest.approx(-GRAVITATIONAL_CONSTANT / r)

    def test_zero_distance_rejected(self):
        with pytest.raises(ValueError):
            m.gravitational_force(1.0, 1.0, 0.0)
        with pytest.raises(ValueError):
            m.gravitational_force_vector(1.0, 1.0, (1, 1, 1), (1, 1, 1))


class TestRotation:
    def test_torque_scalar_matches_vector(self):
        r, f = np.array([0.5, 0.0, 0.0]), np.array([3.0, 4.0, 0.0])
        theta = math.atan2(4.0, 3.0)
        assert m.torque(0.5, 5.0, theta) == pytest.approx(m.torque_vector(r, f).norm())

    @pytest.mark.parametrize(
        ("func", "factor"),
        [
            (m.moment_of_inertia_point_mass, 1.0),
            (m.moment_of_inertia_solid_sphere, 0.4),
            (m.moment_of_inertia_hollow_sphere, 2 / 3),
            (m.moment_of_inertia_solid_cylinder, 0.5),
            (m.moment_of_inertia_hoop, 1.0),
            (m.moment_of_inertia_rod_center, 1 / 12),
            (m.moment_of_inertia_rod_end, 1 / 3),
        ],
    )
    def test_moments_of_inertia(self, func, factor):
        assert func(3.0, 2.0) == pytest.approx(factor * 3.0 * 4.0)

    def test_parallel_axis_theorem(self):
        assert m.parallel_axis(m.moment_of_inertia_rod_center(2.0, 3.0), 2.0, 1.5) == pytest.approx(
            m.moment_of_inertia_rod_end(2.0, 3.0)
        )

    def test_rotational_analogues(self):
        assert m.rotational_kinetic_energy(4.0, 3.0) == pytest.approx(0.5 * m.angular_momentum(4.0, 3.0) * 3.0)
        v = m.tangential_velocity(2.0, 1.5)
        assert m.centripetal_acceleration(v, 1.5) == pytest.approx(2.0**2 * 1.5)

    def test_angular_momentum_vector(self):
        assert m.angular_momentum_vector((1, 0, 0), (0, 2, 0)).z == 2.0


class TestRocketry:
    def test_tsiolkovsky(self):
        assert m.rocket_delta_v(3000.0, math.e * 500.0, 500.0) == pytest.approx(3000.0)

    def test_round_trip_through_mass_ratio_and_propellant(self):
        dv, ve, dry = 9_400.0, 4_400.0, 25_000.0
        propellant = m.rocket_propellant_mass(dv, ve, dry)
        assert (dry + propellant) / dry == pytest.approx(m.rocket_mass_ratio(dv, ve))
        assert m.rocket_delta_v(ve, dry + propellant, dry) == pytest.approx(dv)

    def test_isp_conversions(self):
        assert m.specific_impulse(m.exhaust_velocity_from_isp(350.0)) == pytest.approx(350.0)

    def test_thrust_pressure_term(self):
        assert m.rocket_thrust(10.0, 2000.0, 50_000.0, 101_325.0, 0.5) == pytest.approx(
            20_000.0 + (50_000.0 - 101_325.0) * 0.5
        )

    def test_burn_time(self):
        assert m.burn_time(m.rocket_propellant_mass(1000.0, 3000.0, 100.0), 2.0) > 0

    def test_vectorised_over_stage_masses(self):
        dv = m.rocket_delta_v(3000.0, [2000.0, 4000.0, 8000.0], 1000.0)
        np.testing.assert_allclose(dv, 3000.0 * np.log([2.0, 4.0, 8.0]))

    def test_final_mass_cannot_exceed_initial(self):
        with pytest.raises(ValueError, match="initial_mass"):
            m.rocket_delta_v(3000.0, 500.0, 1000.0)

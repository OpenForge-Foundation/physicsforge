import math

import numpy as np
import pytest

from physicsforge import electromagnetism as em
from physicsforge import numerics
from physicsforge.constants import (
    COULOMB_CONSTANT,
    SPEED_OF_LIGHT,
    VACUUM_PERMEABILITY,
    VACUUM_PERMITTIVITY,
)


class TestCircuits:
    def test_ohms_law_is_self_consistent(self):
        v = em.ohms_law_voltage(2.5, 40.0)
        assert em.ohms_law_current(v, 40.0) == pytest.approx(2.5)
        assert em.ohms_law_resistance(v, 2.5) == pytest.approx(40.0)

    def test_three_power_forms_agree(self):
        i, r = 1.5, 12.0
        v = em.ohms_law_voltage(i, r)
        p = em.electrical_power(v, i)
        assert em.power_from_current(i, r) == pytest.approx(p)
        assert em.power_from_voltage(v, r) == pytest.approx(p)

    def test_voltage_divider(self):
        assert em.voltage_divider(12.0, 1000.0, 2000.0) == pytest.approx(8.0)
        np.testing.assert_allclose(em.voltage_divider(5.0, 10e3, [10e3, 30e3]), [2.5, 3.75])

    def test_divider_matches_series_current(self):
        vin, r1, r2 = 9.0, 470.0, 330.0
        i = em.ohms_law_current(vin, em.series_resistance([r1, r2]))
        assert em.voltage_divider(vin, r1, r2) == pytest.approx(em.ohms_law_voltage(i, r2))

    def test_current_divider_conserves_current(self):
        total, ra, rb = 3.0, 10.0, 20.0
        assert em.current_divider(total, ra, rb) + em.current_divider(total, rb, ra) == pytest.approx(total)

    def test_networks(self):
        assert em.series_resistance([1.0, 2.0, 3.0]) == 6.0
        assert em.parallel_resistance([6.0, 3.0]) == pytest.approx(2.0)
        assert em.series_capacitance([3e-6, 6e-6]) == pytest.approx(2e-6)
        assert em.parallel_capacitance([3e-6, 6e-6]) == pytest.approx(9e-6)

    @pytest.mark.parametrize("bad", [[], [10.0, 0.0], [10.0, -1.0]])
    def test_parallel_rejects_bad_lists(self, bad):
        with pytest.raises(ValueError):
            em.parallel_resistance(bad)

    def test_rc_transients(self):
        tau = em.rc_time_constant(1e3, 1e-3)
        assert em.capacitor_charging_voltage(5.0, tau, 1e3, 1e-3) == pytest.approx(5.0 * (1 - math.exp(-1)))
        assert em.capacitor_discharging_voltage(5.0, tau, 1e3, 1e-3) == pytest.approx(5.0 * math.exp(-1))
        # Charging + discharging curves always sum to the source voltage.
        t = np.linspace(0, 5, 11)
        np.testing.assert_allclose(
            em.capacitor_charging_voltage(5.0, t, 1e3, 1e-3) + em.capacitor_discharging_voltage(5.0, t, 1e3, 1e-3),
            5.0,
        )


class TestElectrostatics:
    def test_coulomb_sign_convention(self):
        assert em.coulomb_force(1e-6, 1e-6, 0.1) > 0
        assert em.coulomb_force(1e-6, -1e-6, 0.1) < 0

    def test_field_potential_and_energy_relations(self):
        q, r = 2e-9, 0.3
        assert em.electric_field_point_charge(q, r) * 1e-9 == pytest.approx(em.coulomb_force(q, 1e-9, r))
        assert em.electric_potential_point_charge(q, r) * 1e-9 == pytest.approx(
            em.electric_potential_energy(q, 1e-9, r)
        )
        # E = -dV/dr
        dv_dr = numerics.derivative(lambda x: em.electric_potential_point_charge(q, x), r)
        assert -dv_dr == pytest.approx(em.electric_field_point_charge(q, r), rel=1e-9)

    def test_vector_forms(self):
        f = em.coulomb_force_vector(1.0, 1.0, (3, 4, 0), (0, 0, 0))
        assert f.norm() == pytest.approx(COULOMB_CONSTANT / 25)
        assert f.x > 0 and f.y > 0
        e = em.electric_field_vector(1.0, (0, 0, 0), (0, 0, 2))
        assert e.z == pytest.approx(COULOMB_CONSTANT / 4)

    def test_capacitor(self):
        c = em.parallel_plate_capacitance(1.0, 1.0, 1.0)
        assert c == pytest.approx(VACUUM_PERMITTIVITY)
        assert em.capacitor_energy(c, 2.0) == pytest.approx(2 * VACUUM_PERMITTIVITY)
        with pytest.raises(ValueError):
            em.parallel_plate_capacitance(0.5, 1.0, 1.0)


class TestMaxwell:
    def test_gauss_law_equals_coulomb_flux(self):
        q = 3e-9
        for r in (0.1, 1.0, 10.0):
            assert em.electric_field_point_charge(q, r) * 4 * math.pi * r**2 == pytest.approx(em.gauss_electric_flux(q))
        assert em.gauss_enclosed_charge(em.gauss_electric_flux(q)) == pytest.approx(q)

    def test_gauss_differential_form_numerically(self):
        # Outside the source, a Coulomb field is divergence-free; a uniform charge density
        # rho gives div E = rho / eps0 inside a uniformly charged ball.
        rho = 1e-6
        inside = lambda r: rho * np.asarray(r) / (3 * VACUUM_PERMITTIVITY)  # noqa: E731
        assert numerics.divergence(inside, (0.1, 0.2, -0.1)) == pytest.approx(rho / VACUUM_PERMITTIVITY, rel=1e-9)

    def test_faraday_lenz_sign(self):
        assert em.faraday_emf(0.2) == pytest.approx(-0.2)
        assert em.faraday_emf(0.2, turns=50) == pytest.approx(-10.0)

    def test_ampere_law_for_wire(self):
        i, r = 4.0, 0.05
        assert em.magnetic_field_straight_wire(i, r) * 2 * math.pi * r == pytest.approx(
            em.ampere_maxwell_circulation(i)
        )

    def test_displacement_current_closes_ampere(self):
        # In a charging capacitor gap, displacement current equals the wire current.
        i = 0.7
        dphi_dt = i / VACUUM_PERMITTIVITY
        assert em.displacement_current(dphi_dt) == pytest.approx(i)
        assert em.ampere_maxwell_circulation(0.0, dphi_dt) == pytest.approx(em.ampere_maxwell_circulation(i, 0.0))

    def test_static_magnetic_field_is_divergence_free(self):
        # B of an infinite wire along z: B = mu0 I / (2 pi) * (-y, x, 0) / (x^2 + y^2).
        def wire_field(r):
            x, y, _ = r
            s2 = x * x + y * y
            return VACUUM_PERMEABILITY * 2.0 / (2 * math.pi) * np.array([-y, x, 0.0]) / s2

        assert abs(numerics.divergence(wire_field, (0.3, 0.4, 0.0))) < 1e-9
        assert np.allclose(np.asarray(numerics.curl(wire_field, (0.3, 0.4, 0.0))), 0.0, atol=1e-9)

    def test_lorentz_force(self):
        f = em.lorentz_force(2.0, electric_field=(1, 0, 0), velocity=(0, 3, 0), magnetic_field=(0, 0, 4))
        assert np.allclose(np.asarray(f), [2.0 + 24.0, 0.0, 0.0])

    def test_speed_of_light_emerges_from_maxwell(self):
        assert em.speed_of_light_from_vacuum_constants() == pytest.approx(SPEED_OF_LIGHT, rel=1e-9)
        assert em.em_wave_speed(2.25) == pytest.approx(SPEED_OF_LIGHT / 1.5)

    def test_field_geometries(self):
        assert em.magnetic_field_solenoid(1000.0, 2.0) == pytest.approx(VACUUM_PERMEABILITY * 2000.0)
        assert em.magnetic_field_loop_center(1.0, 0.5) == pytest.approx(VACUUM_PERMEABILITY)
        assert em.magnetic_flux(2.0, 3.0, math.pi / 2) == pytest.approx(0.0, abs=1e-15)

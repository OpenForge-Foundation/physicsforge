import math

import numpy as np
import pytest

from physicsforge import thermodynamics as th
from physicsforge.constants import GAS_CONSTANT, STANDARD_ATMOSPHERE, STEFAN_BOLTZMANN
from physicsforge.thermodynamics import HeatRod


class TestHeatLaws:
    def test_calorimetry(self):
        assert th.sensible_heat(2.0, 900.0, 5.0) == 9000.0
        assert th.latent_heat(0.5, 2.26e6) == pytest.approx(1.13e6)

    def test_fourier_law_forms_agree(self):
        k, area, dT, L = 0.8, 10.0, 20.0, 0.2
        flux = th.fourier_heat_flux(k, -dT / L)
        assert flux * area == pytest.approx(th.conduction_rate(k, area, dT, L))
        assert th.heat_rate_from_resistance(dT, th.conduction_resistance(L, k, area)) == pytest.approx(
            th.conduction_rate(k, area, dT, L)
        )

    def test_newton_cooling_solves_its_ode(self):
        # dT/dt = -k (T - T_env)
        k, t0, env = 0.05, 90.0, 20.0
        t = 7.0
        h = 1e-4
        dTdt = (th.newton_cooling_temperature(t0, env, k, t + h) - th.newton_cooling_temperature(t0, env, k, t - h)) / (
            2 * h
        )
        assert dTdt == pytest.approx(-k * (th.newton_cooling_temperature(t0, env, k, t) - env), rel=1e-7)

    def test_stefan_boltzmann(self):
        assert th.radiated_power(1.0, 2.0, 300.0) == pytest.approx(2 * STEFAN_BOLTZMANN * 300.0**4)
        assert th.net_radiated_power(0.9, 1.0, 300.0, 300.0) == pytest.approx(0.0, abs=1e-12)
        with pytest.raises(ValueError, match="emissivity"):
            th.radiated_power(1.2, 1.0, 300.0)

    def test_laws_of_thermodynamics(self):
        assert th.internal_energy_change(500.0, 200.0) == 300.0
        assert th.carnot_efficiency(300.0, 500.0) == pytest.approx(0.4)
        assert th.carnot_cop_refrigerator(270.0, 300.0) == pytest.approx(9.0)
        with pytest.raises(ValueError):
            th.carnot_efficiency(500.0, 300.0)

    def test_ideal_gas(self):
        v = th.ideal_gas_volume(1.0, 273.15, STANDARD_ATMOSPHERE)
        assert v == pytest.approx(0.022414, rel=1e-4)
        assert th.ideal_gas_pressure(2.0, 300.0, 0.05) == pytest.approx(2 * GAS_CONSTANT * 300 / 0.05)

    def test_temperature_conversion(self):
        np.testing.assert_allclose(th.kelvin_to_celsius(th.celsius_to_kelvin([-40.0, 0.0, 100.0])), [-40.0, 0.0, 100.0])
        with pytest.raises(ValueError):
            th.celsius_to_kelvin(-274.0)


class TestHeatEquation:
    LENGTH = 1.0
    ALPHA = 1e-2

    def exact_mode(self, n, x, t):
        k = n * math.pi / self.LENGTH
        return np.sin(k * x) * math.exp(-self.ALPHA * k * k * t)

    def test_fourier_series_is_exact_for_mode_superposition(self):
        x = np.linspace(0, self.LENGTH, 201)
        initial = self.exact_mode(1, x, 0) + 0.3 * self.exact_mode(4, x, 0)
        b = th.heat_fourier_coefficients(initial, 6)
        np.testing.assert_allclose(b, [1, 0, 0, 0.3, 0, 0], atol=1e-12)
        sol = th.heat_fourier_solution(b, self.LENGTH, self.ALPHA, 3.0, x)
        np.testing.assert_allclose(sol, self.exact_mode(1, x, 3.0) + 0.3 * self.exact_mode(4, x, 3.0), atol=1e-12)

    def test_heat_rod_methods_agree(self):
        rod = HeatRod(self.LENGTH, self.ALPHA, n_points=101)
        initial = lambda x: np.sin(np.pi * x) + 0.5 * np.sin(2 * np.pi * x)  # noqa: E731
        exact = self.exact_mode(1, rod.grid, 4.0) + 0.5 * self.exact_mode(2, rod.grid, 4.0)
        np.testing.assert_allclose(rod.fourier_solution(initial, 4.0, n_terms=10), exact, atol=1e-12)
        np.testing.assert_allclose(rod.solve(initial, 4.0), exact, atol=5e-4)
        np.testing.assert_allclose(rod.solve(initial, 4.0, method="explicit"), exact, atol=5e-4)

    def test_crank_nicolson_converges_second_order_in_time(self):
        rod = HeatRod(self.LENGTH, self.ALPHA, n_points=401)
        initial = np.sin(np.pi * rod.grid)
        reference = rod.fourier_solution(initial, 5.0, n_terms=1)
        errors = [np.max(np.abs(rod.solve(initial, 5.0, time_step=dt) - reference)) for dt in (1.0, 0.5)]
        assert errors[0] / errors[1] == pytest.approx(4.0, rel=0.1)

    def test_boundary_conditions_and_steady_state(self):
        rod = HeatRod(0.5, 1e-3, left_temperature=373.15, right_temperature=293.15, n_points=21)
        late = rod.solve(293.15, time=1e5, time_step=50.0)
        np.testing.assert_allclose(late, rod.steady_state(), atol=1e-6)
        fourier_late = rod.fourier_solution(293.15, time=1e5, n_terms=19)
        np.testing.assert_allclose(fourier_late, rod.steady_state(), atol=1e-9)
        assert late[0] == 373.15 and late[-1] == 293.15

    def test_explicit_rejects_unstable_step(self):
        rod = HeatRod(1.0, 1.0, n_points=11)
        with pytest.raises(ValueError, match="stability"):
            rod.solve(0.0, time=1.0, method="explicit", time_step=0.1)

    def test_crank_nicolson_stays_bounded_with_huge_steps(self):
        rod = HeatRod(1.0, 1.0, n_points=51)
        out = rod.solve(lambda x: np.sin(np.pi * x), time=10.0, time_step=5.0)
        assert np.all(np.abs(out) <= 1.0)

    def test_validation(self):
        with pytest.raises(ValueError):
            HeatRod(-1.0, 1.0)
        with pytest.raises(ValueError):
            HeatRod(1.0, 1.0, n_points=2)
        rod = HeatRod(1.0, 1.0, n_points=5)
        with pytest.raises(ValueError, match="shape"):
            rod.solve([1.0, 2.0], 1.0)
        with pytest.raises(ValueError, match="method"):
            rod.solve(0.0, 1.0, method="rk4")

    def test_zero_time_returns_initial_with_boundaries(self):
        rod = HeatRod(1.0, 1.0, left_temperature=1.0, n_points=5)
        np.testing.assert_array_equal(rod.solve(0.5, 0.0), [1.0, 0.5, 0.5, 0.5, 0.0])

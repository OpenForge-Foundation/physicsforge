import math

import numpy as np
import pytest

from physicsforge import modern as mp
from physicsforge import numerics
from physicsforge.constants import ELECTRON_VOLT, PLANCK, SPEED_OF_LIGHT, STEFAN_BOLTZMANN

C = SPEED_OF_LIGHT


class TestQuantum:
    def test_planck_einstein_relations(self):
        wavelength = 500e-9
        f = SPEED_OF_LIGHT / wavelength
        assert mp.photon_energy_from_wavelength(wavelength) == pytest.approx(mp.photon_energy_from_frequency(f))
        assert mp.photon_momentum(wavelength) * C == pytest.approx(mp.photon_energy_from_wavelength(wavelength))
        assert mp.wavelength_from_frequency(f) == pytest.approx(wavelength)

    def test_photoelectric_threshold(self):
        phi = 2.0 * ELECTRON_VOLT
        threshold = phi / PLANCK
        assert mp.photoelectric_kinetic_energy(0.99 * threshold, phi) == 0.0
        assert mp.photoelectric_kinetic_energy(2 * threshold, phi) == pytest.approx(phi)

    def test_planck_law_integrates_to_stefan_boltzmann(self):
        temperature = 1500.0
        # Integrate in log-wavelength for accuracy across the spectrum.
        total = numerics.integrate(
            lambda u: mp.planck_spectral_radiance(math.exp(u), temperature) * math.exp(u),
            math.log(1e-7),
            math.log(1e-2),
            tolerance=1e-6,
        )
        assert math.pi * total == pytest.approx(STEFAN_BOLTZMANN * temperature**4, rel=1e-6)

    def test_planck_peak_matches_wien(self):
        temperature = 3000.0
        peak = mp.wien_peak_wavelength(temperature)
        # dB/dlambda = 0 at the Wien peak. Sub-unit variables need an explicit step.
        slope = numerics.derivative(lambda wl: mp.planck_spectral_radiance(wl, temperature), peak, step=peak * 1e-3)
        scale = mp.planck_spectral_radiance(peak, temperature) / peak
        assert abs(slope) / scale < 1e-6

    def test_planck_limits(self):
        assert mp.planck_spectral_radiance(1e-9, 300.0) == 0.0
        assert mp.planck_spectral_radiance(1e-6, 0.0) == 0.0
        # Rayleigh-Jeans limit at long wavelength: B ~ 2 c k T / lambda^4.
        from physicsforge.constants import BOLTZMANN

        wl, T = 1.0, 5000.0
        assert mp.planck_spectral_radiance(wl, T) == pytest.approx(2 * C * BOLTZMANN * T / wl**4, rel=1e-5)

    def test_bohr_atom(self):
        e1 = mp.bohr_energy_level(1)
        assert e1 / ELECTRON_VOLT == pytest.approx(-13.6057, rel=1e-5)
        np.testing.assert_allclose(mp.bohr_energy_level([1, 2, 3]), e1 / np.array([1, 4, 9]))
        # Lyman-alpha (2 -> 1) is ~121.6 nm.
        photon = mp.bohr_energy_level(2) - e1
        assert PLANCK * C / photon == pytest.approx(121.5e-9, rel=2e-3)
        with pytest.raises(ValueError):
            mp.bohr_energy_level(0)

    def test_de_broglie_and_uncertainty(self):
        assert mp.de_broglie_wavelength(1.0, -2.0) == pytest.approx(PLANCK / 2.0)
        assert mp.heisenberg_min_momentum_uncertainty(1e-10) > 0


class TestRelativity:
    def test_lorentz_factor_values(self):
        assert mp.lorentz_factor(0.0) == 1.0
        assert mp.lorentz_factor(0.6 * C) == pytest.approx(1.25)
        assert mp.lorentz_factor(-0.8 * C) == pytest.approx(5 / 3)

    def test_dilation_and_contraction_are_reciprocal(self):
        v = 0.9 * C
        assert mp.time_dilation(1.0, v) * mp.length_contraction(1.0, v) == pytest.approx(1.0)

    def test_energy_momentum_relation(self):
        # E^2 = (pc)^2 + (mc^2)^2
        m, v = 2.0, 0.7 * C
        total = mp.rest_energy(m) + mp.relativistic_kinetic_energy(m, v)
        assert total**2 == pytest.approx((mp.relativistic_momentum(m, v) * C) ** 2 + mp.rest_energy(m) ** 2)

    def test_newtonian_limit_without_cancellation(self):
        for v in (1e-3, 1.0, 300.0):
            assert mp.relativistic_kinetic_energy(1.0, v) == pytest.approx(0.5 * v * v, rel=1e-9)

    @pytest.mark.parametrize("v", [C, 1.1 * C, -C])
    def test_superluminal_rejected(self, v):
        with pytest.raises(ValueError, match="speed of light"):
            mp.lorentz_factor(v)

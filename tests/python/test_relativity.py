import pytest

from physicsforge import relativity as mp
from physicsforge.constants import SPEED_OF_LIGHT

C = SPEED_OF_LIGHT


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

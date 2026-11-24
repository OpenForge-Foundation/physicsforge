"""Package-level contracts: version, public API surface, broadcasting and errors."""

import math

import numpy as np
import pytest

import physicsforge as pf
from physicsforge import constants


def test_version_matches_core():
    assert pf.__version__ == pf._physicsforge.__version__
    major, minor, patch = (int(p) for p in pf.__version__.split("."))
    assert (major, minor, patch) >= (0, 1, 0)


@pytest.mark.parametrize(
    "package", [pf.mechanics, pf.electromagnetism, pf.thermodynamics, pf.modern, pf.numerics, pf.vector]
)
def test_all_exports_exist_and_are_documented(package):
    for name in package.__all__:
        obj = getattr(package, name)
        if callable(obj) and not isinstance(obj, type(pf)):
            assert obj.__doc__, f"{package.__name__}.{name} has no docstring"


def test_scalar_inputs_return_python_float():
    result = pf.mechanics.kinetic_energy(2.0, 3.0)
    assert type(result) is float


def test_broadcasting_matches_numpy_semantics():
    masses = np.array([1.0, 2.0, 3.0])
    velocities = np.array([[1.0], [2.0]])
    result = pf.mechanics.kinetic_energy(masses, velocities)
    assert result.shape == (2, 3)
    np.testing.assert_allclose(result, 0.5 * masses * velocities**2)


def test_invalid_input_raises_value_error_with_message():
    with pytest.raises(ValueError, match="resistance must be finite and positive"):
        pf.electromagnetism.ohms_law_current(5.0, 0.0)


def test_invalid_element_in_array_raises():
    with pytest.raises(ValueError):
        pf.mechanics.gravitational_force(1.0, 1.0, [1.0, -1.0])


def test_nan_is_rejected():
    with pytest.raises(ValueError):
        pf.mechanics.net_force(math.nan, 1.0)


def test_constants_are_consistent():
    assert constants.c == constants.SPEED_OF_LIGHT == 299_792_458.0
    assert math.isclose(constants.c**2 * constants.epsilon_0 * constants.mu_0, 1.0, rel_tol=1e-9)
    assert math.isclose(constants.hbar * 2 * math.pi, constants.h, rel_tol=1e-15)
    assert constants.pi == math.pi
    assert constants.EULER == math.e

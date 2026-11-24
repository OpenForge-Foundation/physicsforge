import copy
import math
import pickle

import numpy as np
import pytest

from physicsforge.vector import (
    Vector3,
    angle_between,
    as_vector,
    cross,
    distance,
    dot,
    norm,
    normalize,
    project,
    triple_product,
)


def test_construction_and_components():
    v = Vector3(1, 2, 3)
    assert (v.x, v.y, v.z) == (1.0, 2.0, 3.0)
    assert Vector3() == Vector3(0, 0, 0)
    v.x = 5
    assert v[0] == 5.0
    assert v[-1] == 3.0
    with pytest.raises(IndexError):
        _ = v[3]


def test_operators():
    a, b = Vector3(1, 2, 3), Vector3(4, 5, 6)
    assert a + b == Vector3(5, 7, 9)
    assert b - a == Vector3(3, 3, 3)
    assert -a == Vector3(-1, -2, -3)
    assert 2 * a == a * 2 == Vector3(2, 4, 6)
    assert a / 2 == Vector3(0.5, 1, 1.5)
    assert a != b


def test_python_protocols():
    v = Vector3(1, 2, 3)
    assert len(v) == 3
    assert list(v) == [1.0, 2.0, 3.0]
    assert v.to_tuple() == (1.0, 2.0, 3.0)
    np.testing.assert_array_equal(np.asarray(v), [1.0, 2.0, 3.0])
    assert np.asarray(v, dtype=np.float32).dtype == np.float32
    assert repr(v) == "Vector3(1.0, 2.0, 3.0)"
    assert pickle.loads(pickle.dumps(v)) == v
    assert copy.deepcopy(v) == v
    with pytest.raises(TypeError):
        hash(v)


def test_products_match_numpy(rng):
    for _ in range(50):
        a, b, c = rng.normal(size=(3, 3))
        assert math.isclose(dot(a, b), float(np.dot(a, b)), rel_tol=1e-14, abs_tol=1e-14)
        np.testing.assert_allclose(np.asarray(cross(a, b)), np.cross(a, b), rtol=1e-14, atol=1e-14)
        assert math.isclose(triple_product(a, b, c), float(np.dot(a, np.cross(b, c))), rel_tol=1e-12, abs_tol=1e-14)


def test_cross_product_identities():
    i, j, k = (1, 0, 0), (0, 1, 0), (0, 0, 1)
    assert cross(i, j) == Vector3(*k)
    assert cross(j, i) == -Vector3(*k)
    a, b = Vector3(1, 2, 3), Vector3(-2, 0.5, 4)
    assert dot(cross(a, b), a) == pytest.approx(0.0, abs=1e-12)
    # |a x b|^2 + (a.b)^2 = |a|^2 |b|^2  (Lagrange's identity)
    assert norm(cross(a, b)) ** 2 + dot(a, b) ** 2 == pytest.approx(norm(a) ** 2 * norm(b) ** 2)


def test_geometry_helpers():
    assert norm((3, 4, 12)) == 13.0
    assert norm(normalize((1, 2, 3))) == pytest.approx(1.0)
    assert angle_between((1, 0, 0), (-1, 0, 0)) == pytest.approx(math.pi)
    assert angle_between((1, 0, 0), (1, 1e-9, 0)) == pytest.approx(1e-9, rel=1e-6)  # stable near parallel
    assert project((2, 3, 4), (0, 2, 0)) == Vector3(0, 3, 0)
    assert distance((0, 0, 0), (1, 2, 2)) == 3.0


def test_errors():
    with pytest.raises(ValueError):
        normalize((0, 0, 0))
    with pytest.raises(ValueError):
        angle_between((0, 0, 0), (1, 0, 0))
    with pytest.raises(ValueError, match="exactly 3 components"):
        as_vector((1, 2))

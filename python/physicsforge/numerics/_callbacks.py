"""Adapters between user callables (NumPy in, NumPy out) and the C++ core's signatures.

The C++ core passes ``std::vector<double>`` (a list in Python) and ``Vector3``
values to callbacks. Users should not have to care: their functions receive
``numpy.ndarray`` arguments and may return any array-like.
"""

from __future__ import annotations

from collections.abc import Callable
from typing import Any

import numpy as np

from physicsforge.vector import Vector3, as_vector


def multivariate(f: Callable[[np.ndarray], Any]) -> Callable[[list[float]], float]:
    """Wrap ``f: ndarray -> scalar`` for the core's ``vector<double> -> double`` signature."""

    def call(x: list[float]) -> float:
        return float(f(np.asarray(x, dtype=np.float64)))

    return call


def vector_valued(f: Callable[[np.ndarray], Any]) -> Callable[[list[float]], list[float]]:
    """Wrap ``F: ndarray -> array_like`` for the core's ``vector<double> -> vector<double>`` signature."""

    def call(x: list[float]) -> list[float]:
        return np.asarray(f(np.asarray(x, dtype=np.float64)), dtype=np.float64).ravel().tolist()

    return call


def vector_field(field: Callable[[np.ndarray], Any]) -> Callable[[Vector3], Vector3]:
    """Wrap ``F: ndarray(3,) -> array_like(3)`` for the core's ``Vec3 -> Vec3`` signature."""

    def call(point: Vector3) -> Vector3:
        return as_vector(field(np.array([point.x, point.y, point.z])), "field(point)")

    return call


def as_point(x: Any) -> list[float]:
    """Validate and convert a point in R^n to a list for the core."""
    array = np.asarray(x, dtype=np.float64)
    if array.ndim != 1 or array.size == 0:
        raise ValueError(f"x must be a non-empty 1-D array (got shape {array.shape})")
    return array.tolist()

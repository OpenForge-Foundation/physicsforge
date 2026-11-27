"""Three-dimensional vectors: :class:`Vector3`, dot and cross products.

:class:`Vector3` is implemented in C++ and supports ``+``, ``-``, scalar ``*`` and
``/``, ``==``, indexing, iteration, pickling and ``numpy.asarray``. The functions
here accept any :data:`~physicsforge._typing.VectorLike` (a ``Vector3``, a
length-3 sequence or a shape ``(3,)`` array).

For large batches of vectors, use NumPy directly (``numpy.cross`` /
``numpy.einsum``) on ``(N, 3)`` arrays; ``Vector3`` is designed for clear,
allocation-free single-vector physics.

Examples
--------
>>> from physicsforge.vector import Vector3, cross, dot
>>> cross([1, 0, 0], [0, 1, 0])
Vector3(0.0, 0.0, 1.0)
>>> dot(Vector3(1, 2, 3), (4, 5, 6))
32.0
>>> Vector3(3, 4, 0).norm()
5.0
"""

from __future__ import annotations

import numpy as np

from physicsforge import _physicsforge as _core
from physicsforge._physicsforge import Vector3
from physicsforge._typing import VectorLike


def as_vector(value: VectorLike, name: str = "vector") -> Vector3:
    """Convert a vector-like value to :class:`Vector3`.

    Parameters
    ----------
    value : Vector3 or array_like
        A ``Vector3``, a length-3 sequence or a shape ``(3,)`` array of numbers.
    name : str, optional
        Argument name used in error messages.

    Returns
    -------
    Vector3

    Raises
    ------
    ValueError
        If ``value`` does not have exactly three components.

    Examples
    --------
    >>> as_vector((1, 2, 3))
    Vector3(1.0, 2.0, 3.0)
    """
    if isinstance(value, Vector3):
        return value
    array = np.asarray(value, dtype=np.float64)
    if array.shape != (3,):
        raise ValueError(f"{name} must have exactly 3 components (got shape {array.shape})")
    return Vector3(float(array[0]), float(array[1]), float(array[2]))


def dot(a: VectorLike, b: VectorLike) -> float:
    r"""Dot (scalar) product :math:`\mathbf{a}\cdot\mathbf{b} = |a||b|\cos\theta`.

    Examples
    --------
    >>> dot((1, 0, 0), (0, 1, 0))
    0.0
    """
    return _core.vec_dot(as_vector(a, "a"), as_vector(b, "b"))


def cross(a: VectorLike, b: VectorLike) -> Vector3:
    r"""Cross (vector) product :math:`\mathbf{a}\times\mathbf{b}`.

    The result is perpendicular to both inputs, with magnitude
    :math:`|a||b|\sin\theta` and direction given by the right-hand rule.

    Examples
    --------
    >>> cross((0, 1, 0), (0, 0, 1))
    Vector3(1.0, 0.0, 0.0)
    """
    return _core.vec_cross(as_vector(a, "a"), as_vector(b, "b"))


def norm(a: VectorLike) -> float:
    r"""Euclidean length :math:`|\mathbf{a}|`, computed without intermediate overflow.

    Examples
    --------
    >>> norm((3, 4, 12))
    13.0
    """
    return _core.vec_norm(as_vector(a, "a"))


def normalize(a: VectorLike) -> Vector3:
    """Unit vector in the direction of ``a``.

    Raises
    ------
    ValueError
        If ``a`` is the zero vector.

    Examples
    --------
    >>> normalize((0, 0, 5))
    Vector3(0.0, 0.0, 1.0)
    """
    return _core.vec_normalized(as_vector(a, "a"))


def angle_between(a: VectorLike, b: VectorLike) -> float:
    r"""Angle between two vectors in radians, in :math:`[0, \pi]`.

    Computed as :math:`\operatorname{atan2}(|a\times b|, a\cdot b)`, which stays
    accurate for nearly (anti)parallel vectors where ``acos`` loses precision.

    Raises
    ------
    ValueError
        If either vector is zero.

    Examples
    --------
    >>> import math
    >>> math.degrees(angle_between((1, 0, 0), (0, 1, 0)))
    90.0
    """
    return _core.vec_angle_between(as_vector(a, "a"), as_vector(b, "b"))


def project(a: VectorLike, onto: VectorLike) -> Vector3:
    r"""Vector projection of ``a`` onto ``onto``: :math:`\frac{a\cdot b}{|b|^2}\,b`.

    Examples
    --------
    >>> project((2, 3, 4), (1, 0, 0))
    Vector3(2.0, 0.0, 0.0)
    """
    return _core.vec_project(as_vector(a, "a"), as_vector(onto, "onto"))


def triple_product(a: VectorLike, b: VectorLike, c: VectorLike) -> float:
    r"""Scalar triple product :math:`a\cdot(b\times c)` (signed parallelepiped volume).

    Examples
    --------
    >>> triple_product((1, 0, 0), (0, 1, 0), (0, 0, 1))
    1.0
    """
    return _core.vec_triple_product(as_vector(a, "a"), as_vector(b, "b"), as_vector(c, "c"))


def distance(a: VectorLike, b: VectorLike) -> float:
    """Euclidean distance between points ``a`` and ``b``.

    Examples
    --------
    >>> distance((0, 0, 0), (1, 2, 2))
    3.0
    """
    return _core.vec_distance(as_vector(a, "a"), as_vector(b, "b"))


__all__ = [
    "Vector3",
    "angle_between",
    "as_vector",
    "cross",
    "distance",
    "dot",
    "norm",
    "normalize",
    "project",
    "triple_product",
]

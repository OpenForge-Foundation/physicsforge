r"""Numerical differentiation and vector calculus.

Derivatives, gradients, Hessians, Jacobians, the Laplacian, and the divergence
and curl of vector fields (for example to check the differential form of
Maxwell's equations, :math:`\nabla\cdot\mathbf{E} = \rho/\varepsilon_0`).

Step sizes are chosen automatically to balance truncation and round-off error.
First derivatives use a fourth-order five-point stencil (relative error around
1e-12 for smooth functions); second derivatives use second-order stencils
(around 1e-7).

User functions receive ``numpy.ndarray`` arguments.
"""

from __future__ import annotations

from collections.abc import Callable
from typing import Any

import numpy as np

from physicsforge import _physicsforge as _core
from physicsforge._typing import ArrayLike, FloatArray, VectorLike
from physicsforge.vector import Vector3, as_vector

# --- Adapters between NumPy-style user callables and the core's list/Vector3 signatures ---


def _multivariate(f: Callable[[np.ndarray], Any]) -> Callable[[list[float]], float]:
    """Wrap ``f: ndarray -> scalar`` for the core's ``vector<double> -> double`` signature."""

    def call(x: list[float]) -> float:
        return float(f(np.asarray(x, dtype=np.float64)))

    return call


def _vector_valued(f: Callable[[np.ndarray], Any]) -> Callable[[list[float]], list[float]]:
    """Wrap ``F: ndarray -> array_like`` for the core's ``vector<double> -> vector<double>`` signature."""

    def call(x: list[float]) -> list[float]:
        return np.asarray(f(np.asarray(x, dtype=np.float64)), dtype=np.float64).ravel().tolist()

    return call


def _vector_field(field: Callable[[np.ndarray], Any]) -> Callable[[Vector3], Vector3]:
    """Wrap ``F: ndarray(3,) -> array_like(3)`` for the core's ``Vec3 -> Vec3`` signature."""

    def call(point: Vector3) -> Vector3:
        return as_vector(field(np.array([point.x, point.y, point.z])), "field(point)")

    return call


def _as_point(x: Any) -> list[float]:
    """Validate and convert a point in R^n to a list for the core."""
    array = np.asarray(x, dtype=np.float64)
    if array.ndim != 1 or array.size == 0:
        raise ValueError(f"x must be a non-empty 1-D array (got shape {array.shape})")
    return array.tolist()


def derivative(f: Callable[[float], float], x: float, step: float | None = None) -> float:
    r"""First derivative :math:`f'(x)` by the five-point central difference.

    .. math:: f'(x) \approx \frac{f(x-2h) - 8f(x-h) + 8f(x+h) - f(x+2h)}{12h}

    Parameters
    ----------
    f : callable
        Scalar function of one variable.
    x : float
        Evaluation point.
    step : float, optional
        Step :math:`h`. When omitted it is ``eps**(1/5) * max(1, |x|)``, an absolute
        floor of ~7e-4; pass an explicit step for quantities much smaller than 1
        (e.g. wavelengths in metres).

    Examples
    --------
    >>> import math
    >>> round(derivative(math.sin, 0.0), 12)
    1.0
    """
    return _core.derivative(f, x, 0.0 if step is None else step)


def second_derivative(f: Callable[[float], float], x: float, step: float | None = None) -> float:
    r"""Second derivative :math:`f''(x)` by the three-point central difference.

    Examples
    --------
    >>> round(second_derivative(lambda x: x**3, 2.0), 6)
    12.0
    """
    return _core.second_derivative(f, x, 0.0 if step is None else step)


def gradient(f: Callable[[np.ndarray], float], x: ArrayLike) -> FloatArray:
    r"""Gradient :math:`\nabla f(\mathbf{x})` of a scalar function of several variables.

    Parameters
    ----------
    f : callable
        ``f(x) -> float`` where ``x`` is a 1-D ``numpy.ndarray``.
    x : array_like
        Evaluation point, shape ``(n,)``.

    Returns
    -------
    numpy.ndarray
        Shape ``(n,)``.

    Examples
    --------
    >>> import numpy as np
    >>> np.round(gradient(lambda v: v[0]**2 + 3 * v[1], [2.0, 5.0]), 10)
    array([4., 3.])
    """
    return _core.gradient(_multivariate(f), _as_point(x))


def hessian(f: Callable[[np.ndarray], float], x: ArrayLike) -> FloatArray:
    r"""Hessian matrix :math:`H_{ij} = \partial^2 f / \partial x_i \partial x_j`.

    The result is exactly symmetric by construction. Its definiteness classifies
    critical points (positive definite: local minimum), and it drives Newton's
    method for optimisation (:func:`physicsforge.numerics.newton_minimize`).

    Returns
    -------
    numpy.ndarray
        Shape ``(n, n)``.

    Examples
    --------
    >>> import numpy as np
    >>> np.round(hessian(lambda v: v[0]**2 * v[1] + 3 * v[1]**2, [1.5, -2.0]), 5)
    array([[-4.,  3.],
           [ 3.,  6.]])
    """
    return _core.hessian(_multivariate(f), _as_point(x))


def jacobian(f: Callable[[np.ndarray], ArrayLike], x: ArrayLike) -> FloatArray:
    r"""Jacobian matrix :math:`J_{ij} = \partial F_i / \partial x_j` of :math:`F:\mathbb{R}^n\to\mathbb{R}^m`.

    Returns
    -------
    numpy.ndarray
        Shape ``(m, n)``.

    Examples
    --------
    >>> import numpy as np
    >>> np.round(jacobian(lambda v: [v[0] * v[1], v[0] + v[1]], [2.0, 3.0]), 10)
    array([[3., 2.],
           [1., 1.]])
    """
    return _core.jacobian(_vector_valued(f), _as_point(x))


def laplacian(f: Callable[[np.ndarray], float], x: ArrayLike) -> float:
    r"""Laplacian :math:`\nabla^2 f = \sum_i \partial^2 f/\partial x_i^2` (the operator in the heat equation).

    Examples
    --------
    >>> round(laplacian(lambda v: v[0]**2 + v[1]**2 + v[2]**2, [1.0, 2.0, 3.0]), 6)
    6.0
    """
    return _core.laplacian(_multivariate(f), _as_point(x))


def divergence(field: Callable[[np.ndarray], Any], point: VectorLike) -> float:
    r"""Divergence :math:`\nabla\cdot\mathbf{F}` of a 3-D vector field at a point.

    Parameters
    ----------
    field : callable
        ``field(r) -> array_like`` of 3 components, where ``r`` is a shape ``(3,)`` array.
    point : Vector3 or array_like
        Evaluation point.

    Examples
    --------
    Gauss's law in differential form: outside a point charge the Coulomb field is
    divergence-free.

    >>> import numpy as np
    >>> coulomb = lambda r: r / np.linalg.norm(r) ** 3
    >>> abs(divergence(coulomb, (1.0, 2.0, 2.0))) < 1e-9
    True
    """
    return _core.divergence(_vector_field(field), as_vector(point, "point"))


def curl(field: Callable[[np.ndarray], Any], point: VectorLike) -> Vector3:
    r"""Curl :math:`\nabla\times\mathbf{F}` of a 3-D vector field at a point.

    Examples
    --------
    A rigid rotation :math:`\mathbf{v} = \boldsymbol\omega\times\mathbf{r}` has curl :math:`2\boldsymbol\omega`.

    >>> import numpy as np
    >>> c = curl(lambda r: np.cross([0.0, 0.0, 1.0], r), (0.3, -0.7, 2.0))
    >>> [round(v, 9) for v in c]
    [0.0, 0.0, 2.0]
    """
    return _core.curl(_vector_field(field), as_vector(point, "point"))


__all__ = ["curl", "derivative", "divergence", "gradient", "hessian", "jacobian", "laplacian", "second_derivative"]

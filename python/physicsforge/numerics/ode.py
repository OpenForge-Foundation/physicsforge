r"""Fixed-step integrators for ordinary differential equations.

===========================  =====  ==============================================================
Method                       Order  Use it for
===========================  =====  ==============================================================
:func:`euler`                1      Teaching and quick sketches; energy drifts in oscillators.
:func:`semi_implicit_euler`  1      Mechanics :math:`\ddot x = a(t, x)`; symplectic, energy stays bounded.
:func:`runge_kutta4`         4      General smooth problems where accuracy matters.
===========================  =====  ==============================================================

Every solver takes ``steps`` equal steps over ``t_span = (t0, t_end)`` and returns the whole trajectory.
"""

from __future__ import annotations

from collections.abc import Callable
from typing import Any

import numpy as np

from physicsforge import _physicsforge as _core
from physicsforge._physicsforge import OdeSolution, SecondOrderSolution
from physicsforge._typing import ArrayLike


def _state(values: ArrayLike, name: str) -> list[float]:
    array = np.atleast_1d(np.asarray(values, dtype=np.float64))
    if array.ndim != 1:
        raise ValueError(f"{name} must be a scalar or 1-D array (got shape {array.shape})")
    return array.tolist()


def _span(t_span: tuple[float, float]) -> tuple[float, float]:
    t0, t_end = t_span
    return float(t0), float(t_end)


def _rhs(f: Callable[[float, np.ndarray], Any]) -> Callable[[float, list[float]], list[float]]:
    """Adapt ``f(t, y: ndarray) -> array_like`` to the core's ``(double, vector) -> vector`` signature."""

    def call(t: float, y: list[float]) -> list[float]:
        return np.atleast_1d(np.asarray(f(t, np.asarray(y, dtype=np.float64)), dtype=np.float64)).ravel().tolist()

    return call


def euler(
    f: Callable[[float, np.ndarray], ArrayLike], y0: ArrayLike, t_span: tuple[float, float], steps: int
) -> OdeSolution:
    r"""Forward Euler method, :math:`\mathbf y_{n+1} = \mathbf y_n + h\,\mathbf f(t_n, \mathbf y_n)`.

    First order: halving the step halves the error. Simple and instructive, but on an undamped
    oscillator it adds energy every step, so prefer :func:`semi_implicit_euler` for mechanics and
    :func:`runge_kutta4` when accuracy matters.

    Parameters
    ----------
    f : callable
        Right-hand side ``f(t, y)`` returning an array the same length as ``y``.
    y0 : float or array_like
        Initial state.
    t_span : tuple of float
        ``(t0, t_end)`` with ``t_end > t0``.
    steps : int
        Number of equal steps (at least 1).

    Returns
    -------
    OdeSolution
        ``t`` of shape ``(steps + 1,)`` and ``y`` of shape ``(steps + 1, len(y0))``.

    Raises
    ------
    ValueError
        If the arguments are invalid or ``f`` returns the wrong number of values.

    Examples
    --------
    Exponential decay :math:`y' = -y` over one time unit:

    >>> import math
    >>> sol = euler(lambda t, y: -y, 1.0, (0.0, 1.0), steps=1000)
    >>> round(float(sol.y[-1, 0]), 4), round(math.exp(-1), 4)
    (0.3677, 0.3679)
    """
    t0, t_end = _span(t_span)
    return _core.euler(_rhs(f), t0, _state(y0, "y0"), t_end, int(steps))


def runge_kutta4(
    f: Callable[[float, np.ndarray], ArrayLike], y0: ArrayLike, t_span: tuple[float, float], steps: int
) -> OdeSolution:
    r"""Classical fourth-order Runge–Kutta method.

    Four evaluations of ``f`` per step; halving the step reduces the error about sixteenfold.

    Parameters
    ----------
    f : callable
        Right-hand side ``f(t, y)`` returning an array the same length as ``y``.
    y0 : float or array_like
        Initial state.
    t_span : tuple of float
        ``(t0, t_end)`` with ``t_end > t0``.
    steps : int
        Number of equal steps (at least 1).

    Returns
    -------
    OdeSolution
        ``t`` of shape ``(steps + 1,)`` and ``y`` of shape ``(steps + 1, len(y0))``.

    Raises
    ------
    ValueError
        If the arguments are invalid or ``f`` returns the wrong number of values.

    Examples
    --------
    A harmonic oscillator written as the first-order system :math:`(x, v)' = (v, -x)`:

    >>> import math
    >>> sol = runge_kutta4(lambda t, y: [y[1], -y[0]], [1.0, 0.0], (0.0, math.pi), steps=100)
    >>> round(float(sol.y[-1, 0]), 8)
    -1.0
    """
    t0, t_end = _span(t_span)
    return _core.runge_kutta4(_rhs(f), t0, _state(y0, "y0"), t_end, int(steps))


def semi_implicit_euler(
    acceleration: Callable[[float, np.ndarray], ArrayLike],
    x0: ArrayLike,
    v0: ArrayLike,
    t_span: tuple[float, float],
    steps: int,
) -> SecondOrderSolution:
    r"""Semi-implicit (symplectic) Euler for :math:`\ddot{\mathbf x} = \mathbf a(t, \mathbf x)`.

    Updates velocity first, then position with the new velocity:

    .. math::
        \mathbf v_{n+1} = \mathbf v_n + h\,\mathbf a(t_n, \mathbf x_n), \qquad
        \mathbf x_{n+1} = \mathbf x_n + h\,\mathbf v_{n+1}.

    Same cost and order as :func:`euler`, but symplectic: with conservative forces the energy error
    stays bounded, so orbits and oscillators stay stable over long runs.

    Parameters
    ----------
    acceleration : callable
        ``acceleration(t, x)`` returning an array the same length as ``x``.
    x0 : float or array_like
        Initial position.
    v0 : float or array_like
        Initial velocity, same length as ``x0``.
    t_span : tuple of float
        ``(t0, t_end)`` with ``t_end > t0``.
    steps : int
        Number of equal steps (at least 1).

    Returns
    -------
    SecondOrderSolution
        ``t`` of shape ``(steps + 1,)``; ``x`` and ``v`` of shape ``(steps + 1, len(x0))``.

    Raises
    ------
    ValueError
        If the arguments are invalid or ``acceleration`` returns the wrong number of values.

    Examples
    --------
    A spring (:math:`\ddot x = -x`) run for 100 periods keeps its energy near 0.5:

    >>> import math
    >>> sol = semi_implicit_euler(lambda t, x: -x, 1.0, 0.0, (0.0, 200 * math.pi), steps=20_000)
    >>> energy = 0.5 * (sol.x[:, 0] ** 2 + sol.v[:, 0] ** 2)
    >>> bool(abs(energy - 0.5).max() < 0.05)
    True
    """
    t0, t_end = _span(t_span)
    return _core.semi_implicit_euler(_rhs(acceleration), t0, _state(x0, "x0"), _state(v0, "v0"), t_end, int(steps))


__all__ = ["OdeSolution", "SecondOrderSolution", "euler", "runge_kutta4", "semi_implicit_euler"]

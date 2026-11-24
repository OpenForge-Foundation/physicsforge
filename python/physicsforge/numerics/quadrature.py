"""Numerical quadrature for sampled data and for callables."""

from __future__ import annotations

from collections.abc import Callable

from physicsforge import _physicsforge as _core
from physicsforge._typing import ArrayLike


def trapezoid(samples: ArrayLike, spacing: float) -> float:
    r"""Composite trapezoidal rule over uniformly spaced samples.

    Parameters
    ----------
    samples : array_like
        Function values on a uniform grid (at least 2).
    spacing : float
        Grid spacing :math:`\Delta x` > 0.

    Examples
    --------
    >>> trapezoid([0.0, 1.0, 2.0], 1.0)
    2.0
    """
    return _core.trapezoid(samples, spacing)


def simpson(samples: ArrayLike, spacing: float) -> float:
    r"""Composite Simpson's rule over uniformly spaced samples (fourth-order accurate).

    Works for any sample count: an odd number of intervals is closed with
    Simpson's 3/8 rule, so the rule stays exact for cubics.

    Examples
    --------
    >>> simpson([x**3 for x in range(10)], 1.0)   # integral of x^3 on [0, 9]
    1640.25
    """
    return _core.simpson(samples, spacing)


def integrate(
    f: Callable[[float], float], a: float, b: float, *, tolerance: float = 1e-10, max_depth: int = 50
) -> float:
    r"""Adaptive Simpson quadrature of :math:`\int_a^b f(x)\,dx`.

    Parameters
    ----------
    f : callable
        Integrand of one variable.
    a, b : float
        Finite limits (``a > b`` gives the negated integral).
    tolerance : float, optional
        Absolute error target. Default 1e-10.
    max_depth : int, optional
        Maximum bisection depth. Default 50.

    Examples
    --------
    >>> import math
    >>> round(integrate(math.sin, 0.0, math.pi), 12)
    2.0
    """
    return _core.integrate(f, a, b, tolerance, max_depth)


__all__ = ["integrate", "simpson", "trapezoid"]

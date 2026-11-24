r"""Root finding: Newton-Raphson (scalar and multivariate) and bisection.

Newton's method iterates

.. math:: x_{k+1} = x_k - \frac{f(x_k)}{f'(x_k)}

and converges quadratically near a simple root. The multivariate form solves
:math:`J(\mathbf{x}_k)\,\Delta\mathbf{x} = -F(\mathbf{x}_k)` at each step.

By default a solver that misses its tolerance raises
:class:`~physicsforge.exceptions.ConvergenceError` (carrying the last result).
Pass ``raise_on_failure=False`` to inspect ``result.converged`` yourself.
"""

from __future__ import annotations

from collections.abc import Callable

import numpy as np

from physicsforge import _physicsforge as _core
from physicsforge._typing import ArrayLike
from physicsforge.exceptions import ConvergenceError
from physicsforge.numerics import _callbacks

RootResult = _core.RootResult
SystemRootResult = _core.SystemRootResult

DEFAULT_TOLERANCE = 1e-12
DEFAULT_MAX_ITERATIONS = 100


def _check(result, name: str, raise_on_failure: bool):
    if raise_on_failure and not result.converged:
        raise ConvergenceError(f"{name} did not converge after {result.iterations} iterations", result)
    return result


def newton(
    f: Callable[[float], float],
    x0: float,
    fprime: Callable[[float], float] | None = None,
    *,
    tolerance: float = DEFAULT_TOLERANCE,
    max_iterations: int = DEFAULT_MAX_ITERATIONS,
    raise_on_failure: bool = True,
) -> RootResult:
    r"""Find a root of :math:`f` with the Newton-Raphson method.

    Parameters
    ----------
    f : callable
        Scalar function of one variable.
    x0 : float
        Initial guess.
    fprime : callable, optional
        Derivative :math:`f'`. If omitted it is computed numerically (five-point stencil).
    tolerance : float, optional
        Stop when :math:`|\Delta x| \le \text{tol}\cdot\max(1, |x|)`. Default 1e-12.
    max_iterations : int, optional
        Iteration cap. Default 100.
    raise_on_failure : bool, optional
        Raise :class:`ConvergenceError` if not converged (default ``True``).

    Returns
    -------
    RootResult
        With attributes ``root``, ``residual``, ``iterations`` and ``converged``.

    Raises
    ------
    ConvergenceError
        If the iteration fails (e.g. hits a stationary point) and ``raise_on_failure``.

    Examples
    --------
    >>> result = newton(lambda x: x**2 - 2.0, x0=1.0, fprime=lambda x: 2.0 * x)
    >>> round(result.root, 15)
    1.414213562373095
    >>> result.converged
    True
    """
    if fprime is None:
        result = _core.newton_raphson_numeric(f, x0, tolerance, max_iterations)
    else:
        result = _core.newton_raphson(f, fprime, x0, tolerance, max_iterations)
    return _check(result, "newton", raise_on_failure)


def bisect(
    f: Callable[[float], float],
    a: float,
    b: float,
    *,
    tolerance: float = DEFAULT_TOLERANCE,
    max_iterations: int = 200,
    raise_on_failure: bool = True,
) -> RootResult:
    r"""Find a root of :math:`f` in a sign-changing bracket :math:`[a, b]` by bisection.

    Slower than Newton (linear convergence) but guaranteed to converge.

    Raises
    ------
    ValueError
        If :math:`f(a)` and :math:`f(b)` have the same sign.

    Examples
    --------
    >>> round(bisect(lambda x: x**3 - 2 * x - 5, 2.0, 3.0).root, 10)
    2.0945514815
    """
    return _check(_core.bisection(f, a, b, tolerance, max_iterations), "bisect", raise_on_failure)


def newton_system(
    f: Callable[[np.ndarray], ArrayLike],
    x0: ArrayLike,
    *,
    tolerance: float = DEFAULT_TOLERANCE,
    max_iterations: int = DEFAULT_MAX_ITERATIONS,
    raise_on_failure: bool = True,
) -> SystemRootResult:
    r"""Solve a square nonlinear system :math:`F(\mathbf{x}) = \mathbf{0}` with Newton's method.

    The Jacobian is evaluated numerically and each step solves
    :math:`J\,\Delta\mathbf{x} = -F` by Gaussian elimination with partial pivoting.

    Parameters
    ----------
    f : callable
        ``f(x) -> array_like`` with ``len(f(x)) == len(x)``.
    x0 : array_like
        Initial guess, shape ``(n,)``.

    Returns
    -------
    SystemRootResult
        With attributes ``root`` (ndarray), ``residual_norm``, ``iterations``, ``converged``.

    Examples
    --------
    Intersection of the circle :math:`x^2 + y^2 = 4` with the line :math:`y = x`:

    >>> import numpy as np
    >>> sol = newton_system(lambda v: [v[0]**2 + v[1]**2 - 4, v[1] - v[0]], [1.0, 0.5])
    >>> np.round(sol.root, 12)
    array([1.41421356, 1.41421356])
    """
    result = _core.newton_system(_callbacks.vector_valued(f), _callbacks.as_point(x0), tolerance, max_iterations)
    return _check(result, "newton_system", raise_on_failure)


__all__ = ["RootResult", "SystemRootResult", "bisect", "newton", "newton_system"]

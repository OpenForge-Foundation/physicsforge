r"""Unconstrained minimisation with Newton's method on the gradient and Hessian.

Each iteration solves

.. math:: H(\mathbf{x}_k)\,\mathbf{p} = -\nabla f(\mathbf{x}_k), \qquad
          \mathbf{x}_{k+1} = \mathbf{x}_k + t\,\mathbf{p}

with :math:`t` from an Armijo backtracking line search. If the Hessian is
singular or not positive definite, the step falls back to steepest descent, so
the method makes progress from poor starting points and converges quadratically
near a strict local minimum.
"""

from __future__ import annotations

from collections.abc import Callable

import numpy as np

from physicsforge import _physicsforge as _core
from physicsforge._typing import ArrayLike
from physicsforge.exceptions import ConvergenceError
from physicsforge.numerics import _callbacks

OptimizeResult = _core.OptimizeResult


def newton_minimize(
    f: Callable[[np.ndarray], float],
    x0: ArrayLike,
    *,
    tolerance: float = 1e-8,
    max_iterations: int = 100,
    raise_on_failure: bool = True,
) -> OptimizeResult:
    r"""Minimise :math:`f:\mathbb{R}^n\to\mathbb{R}` with a damped Newton method.

    Parameters
    ----------
    f : callable
        Objective ``f(x) -> float`` where ``x`` is a 1-D ``numpy.ndarray``.
    x0 : array_like
        Starting point, shape ``(n,)``.
    tolerance : float, optional
        Converged when :math:`\|\nabla f\| \le \text{tol}` or the step is below
        :math:`\text{tol}\cdot\max(1, \|x\|)`. Default 1e-8 (gradients are numerical).
    max_iterations : int, optional
        Iteration cap. Default 100.
    raise_on_failure : bool, optional
        Raise :class:`~physicsforge.exceptions.ConvergenceError` if not converged.

    Returns
    -------
    OptimizeResult
        With attributes ``x`` (ndarray), ``value``, ``gradient``, ``iterations``, ``converged``.

    Examples
    --------
    The Rosenbrock function has its minimum at (1, 1):

    >>> import numpy as np
    >>> rosen = lambda v: (1 - v[0])**2 + 100 * (v[1] - v[0]**2)**2
    >>> result = newton_minimize(rosen, x0=[-1.2, 1.0])
    >>> np.round(result.x, 6)
    array([1., 1.])
    """
    result = _core.newton_minimize(_callbacks.multivariate(f), _callbacks.as_point(x0), tolerance, max_iterations)
    if raise_on_failure and not result.converged:
        raise ConvergenceError(f"newton_minimize did not converge after {result.iterations} iterations", result)
    return result


__all__ = ["OptimizeResult", "newton_minimize"]

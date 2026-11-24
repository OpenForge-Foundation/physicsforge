"""Exception types raised by physicsforge.

Invalid physical inputs (negative mass, zero distance, speeds at or above *c*, ...)
raise the built-in :class:`ValueError`. Iterative solvers that fail to converge
raise :class:`ConvergenceError` unless called with ``raise_on_failure=False``.
"""

from __future__ import annotations

from typing import Any


class PhysicsForgeError(Exception):
    """Base class for physicsforge-specific errors."""


class ConvergenceError(PhysicsForgeError, RuntimeError):
    """An iterative solver stopped without meeting its tolerance.

    Attributes
    ----------
    result : object
        The solver's result object (e.g. :class:`~physicsforge.numerics.RootResult`)
        holding the last iterate, so callers can inspect or reuse it.
    """

    def __init__(self, message: str, result: Any) -> None:
        super().__init__(message)
        self.result = result


__all__ = ["ConvergenceError", "PhysicsForgeError"]

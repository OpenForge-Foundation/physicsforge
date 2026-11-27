"""Special relativity: Lorentz factor, time dilation, length contraction, mass-energy.

Speeds must satisfy :math:`|v| < c`; otherwise ``ValueError`` is raised.
"""

from __future__ import annotations

from physicsforge import _physicsforge as _core
from physicsforge._typing import ArrayLike, FloatOrArray


def lorentz_factor(velocity: ArrayLike) -> FloatOrArray:
    r"""Lorentz factor :math:`\gamma = 1/\sqrt{1 - v^2/c^2}`.

    Examples
    --------
    >>> from physicsforge.constants import SPEED_OF_LIGHT as c
    >>> lorentz_factor(0.6 * c)
    1.25
    """
    return _core.lorentz_factor(velocity)


def time_dilation(proper_time: ArrayLike, velocity: ArrayLike) -> FloatOrArray:
    r"""Dilated time interval :math:`\Delta t = \gamma\,\Delta t_0`.

    Examples
    --------
    >>> from physicsforge.constants import SPEED_OF_LIGHT as c
    >>> time_dilation(1.0, 0.6 * c)
    1.25
    """
    return _core.time_dilation(proper_time, velocity)


def length_contraction(proper_length: ArrayLike, velocity: ArrayLike) -> FloatOrArray:
    r"""Contracted length :math:`L = L_0/\gamma`.

    Examples
    --------
    >>> from physicsforge.constants import SPEED_OF_LIGHT as c
    >>> round(length_contraction(1.0, 0.6 * c), 12)
    0.8
    """
    return _core.length_contraction(proper_length, velocity)


def rest_energy(mass: ArrayLike) -> FloatOrArray:
    r"""Mass-energy equivalence, :math:`E = mc^2`.

    Examples
    --------
    >>> rest_energy(1.0)
    8.987551787368176e+16
    """
    return _core.rest_energy(mass)


def relativistic_momentum(mass: ArrayLike, velocity: ArrayLike) -> FloatOrArray:
    r"""Relativistic momentum :math:`p = \gamma m v`."""
    return _core.relativistic_momentum(mass, velocity)


def relativistic_kinetic_energy(mass: ArrayLike, velocity: ArrayLike) -> FloatOrArray:
    r"""Relativistic kinetic energy :math:`K = (\gamma - 1) mc^2`.

    Evaluated in a cancellation-free form, so it matches :math:`\tfrac12 mv^2` at
    everyday speeds to full precision.

    Examples
    --------
    >>> relativistic_kinetic_energy(1.0, 1.0)
    0.5
    """
    return _core.relativistic_kinetic_energy(mass, velocity)


__all__ = [
    "length_contraction",
    "lorentz_factor",
    "relativistic_kinetic_energy",
    "relativistic_momentum",
    "rest_energy",
    "time_dilation",
]

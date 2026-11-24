"""Newtonian dynamics: Newton's laws, momentum, work, energy and power.

All functions broadcast over NumPy arrays (computed in C++) and return a ``float``
for scalar inputs. Units are SI; angles are in radians.
"""

from __future__ import annotations

from physicsforge import _physicsforge as _core
from physicsforge._typing import ArrayLike, FloatOrArray
from physicsforge.constants import STANDARD_GRAVITY


def net_force(mass: ArrayLike, acceleration: ArrayLike) -> FloatOrArray:
    r"""Newton's second law, :math:`F = m a`.

    Parameters
    ----------
    mass : float or array_like
        Mass :math:`m` [kg]. Must be non-negative.
    acceleration : float or array_like
        Acceleration :math:`a` [m/s^2].

    Returns
    -------
    float or numpy.ndarray
        Net force [N].

    Raises
    ------
    ValueError
        If ``mass`` is negative or any input is not finite.

    Examples
    --------
    >>> net_force(mass=2.0, acceleration=3.0)
    6.0
    """
    return _core.net_force(mass, acceleration)


def acceleration_from_force(force: ArrayLike, mass: ArrayLike) -> FloatOrArray:
    r"""Acceleration from Newton's second law, :math:`a = F / m`.

    Parameters
    ----------
    force : float or array_like
        Net force [N].
    mass : float or array_like
        Mass [kg]. Must be positive.

    Examples
    --------
    >>> acceleration_from_force(10.0, 4.0)
    2.5
    """
    return _core.acceleration_from_force(force, mass)


def weight(mass: ArrayLike, gravity: ArrayLike = STANDARD_GRAVITY) -> FloatOrArray:
    r"""Weight near a planetary surface, :math:`W = m g`.

    Parameters
    ----------
    mass : float or array_like
        Mass [kg]. Must be non-negative.
    gravity : float or array_like, optional
        Local gravitational acceleration [m/s^2]. Defaults to standard gravity.

    Examples
    --------
    >>> round(weight(10.0), 10)
    98.0665
    """
    return _core.weight(mass, gravity)


def momentum(mass: ArrayLike, velocity: ArrayLike) -> FloatOrArray:
    r"""Linear momentum :math:`p = m v`.

    Examples
    --------
    >>> momentum(2.0, -3.0)
    -6.0
    """
    return _core.momentum(mass, velocity)


def impulse(force: ArrayLike, duration: ArrayLike) -> FloatOrArray:
    r"""Impulse of a constant force, :math:`J = F \Delta t` (equal to the change in momentum).

    Examples
    --------
    >>> impulse(50.0, 0.1)
    5.0
    """
    return _core.impulse(force, duration)


def kinetic_energy(mass: ArrayLike, velocity: ArrayLike) -> FloatOrArray:
    r"""Translational kinetic energy :math:`K = \tfrac{1}{2} m v^2`.

    Examples
    --------
    >>> kinetic_energy(2.0, 3.0)
    9.0
    """
    return _core.kinetic_energy(mass, velocity)


def potential_energy(mass: ArrayLike, height: ArrayLike, gravity: ArrayLike = STANDARD_GRAVITY) -> FloatOrArray:
    r"""Gravitational potential energy in a uniform field, :math:`U = m g h`.

    For the two-body (inverse-square) form use
    :func:`physicsforge.mechanics.gravitation.gravitational_potential_energy`.

    Examples
    --------
    >>> potential_energy(2.0, 10.0, gravity=10.0)
    200.0
    """
    return _core.potential_energy(mass, height, gravity)


def work(force: ArrayLike, distance: ArrayLike, angle: ArrayLike = 0.0) -> FloatOrArray:
    r"""Work done by a constant force, :math:`W = F d \cos\theta`.

    Parameters
    ----------
    force : float or array_like
        Force magnitude [N].
    distance : float or array_like
        Distance moved [m]. Must be non-negative.
    angle : float or array_like, optional
        Angle between force and displacement [rad]. Defaults to 0 (parallel).

    Examples
    --------
    >>> work(10.0, 3.0)
    30.0
    """
    return _core.work(force, distance, angle)


def power(work: ArrayLike, time: ArrayLike) -> FloatOrArray:
    r"""Average power :math:`P = W / t`.

    Examples
    --------
    >>> power(100.0, 4.0)
    25.0
    """
    return _core.power(work, time)


def spring_force(spring_constant: ArrayLike, extension: ArrayLike) -> FloatOrArray:
    r"""Hooke's law restoring force, :math:`F = -k x`.

    Examples
    --------
    >>> spring_force(100.0, 0.1)
    -10.0
    """
    return _core.spring_force(spring_constant, extension)


def spring_potential_energy(spring_constant: ArrayLike, extension: ArrayLike) -> FloatOrArray:
    r"""Elastic potential energy of an ideal spring, :math:`U = \tfrac{1}{2} k x^2`.

    Examples
    --------
    >>> spring_potential_energy(100.0, 0.1)
    0.5
    """
    return _core.spring_potential_energy(spring_constant, extension)


def friction_force(coefficient: ArrayLike, normal_force: ArrayLike) -> FloatOrArray:
    r"""Coulomb (kinetic) friction magnitude, :math:`f = \mu N`.

    Examples
    --------
    >>> friction_force(0.5, 20.0)
    10.0
    """
    return _core.friction_force(coefficient, normal_force)


__all__ = [
    "acceleration_from_force",
    "friction_force",
    "impulse",
    "kinetic_energy",
    "momentum",
    "net_force",
    "potential_energy",
    "power",
    "spring_force",
    "spring_potential_energy",
    "weight",
    "work",
]

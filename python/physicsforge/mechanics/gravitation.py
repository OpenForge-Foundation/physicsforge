"""Newton's law of universal gravitation and two-body orbital mechanics.

Orbital helpers treat the orbiting body's mass as negligible relative to the
central mass :math:`M` (restricted two-body problem). Units are SI.
"""

from __future__ import annotations

from physicsforge import _physicsforge as _core
from physicsforge._typing import ArrayLike, FloatOrArray, VectorLike
from physicsforge.vector import Vector3, as_vector


def gravitational_force(mass1: ArrayLike, mass2: ArrayLike, distance: ArrayLike) -> FloatOrArray:
    r"""Magnitude of the gravitational attraction between two point masses.

    .. math:: F = G \frac{m_1 m_2}{r^2}

    Parameters
    ----------
    mass1, mass2 : float or array_like
        Masses [kg]. Must be non-negative.
    distance : float or array_like
        Centre-to-centre separation :math:`r` [m]. Must be positive.

    Returns
    -------
    float or numpy.ndarray
        Force magnitude [N].

    Raises
    ------
    ValueError
        If a mass is negative or ``distance`` is not positive.

    Examples
    --------
    >>> f = gravitational_force(5.9722e24, 70.0, 6.371e6)   # a person on Earth
    >>> round(f, 1)
    687.4
    """
    return _core.gravitational_force(mass1, mass2, distance)


def gravitational_force_vector(mass1: float, mass2: float, position1: VectorLike, position2: VectorLike) -> Vector3:
    r"""Gravitational force on body 1 due to body 2.

    .. math:: \mathbf{F}_{12} = -G m_1 m_2 \frac{\mathbf{r}_1 - \mathbf{r}_2}{|\mathbf{r}_1 - \mathbf{r}_2|^3}

    Raises
    ------
    ValueError
        If the bodies coincide or a mass is negative.

    Examples
    --------
    >>> f = gravitational_force_vector(1.0, 1.0, (2.0, 0.0, 0.0), (0.0, 0.0, 0.0))
    >>> f.x < 0 and f.y == 0.0   # points toward body 2
    True
    """
    return _core.gravitational_force_vector(
        mass1, mass2, as_vector(position1, "position1"), as_vector(position2, "position2")
    )


def gravitational_field_strength(mass: ArrayLike, distance: ArrayLike) -> FloatOrArray:
    r"""Gravitational field strength of a point or spherical mass, :math:`g = G M / r^2`.

    Examples
    --------
    >>> round(gravitational_field_strength(5.9722e24, 6.371e6), 2)
    9.82
    """
    return _core.gravitational_field_strength(mass, distance)


def gravitational_potential_energy(mass1: ArrayLike, mass2: ArrayLike, distance: ArrayLike) -> FloatOrArray:
    r"""Two-body gravitational potential energy, :math:`U = -G m_1 m_2 / r` (zero at infinity).

    Examples
    --------
    >>> gravitational_potential_energy(1.0, 1.0, 1.0) < 0
    True
    """
    return _core.gravitational_potential_energy(mass1, mass2, distance)


def escape_velocity(central_mass: ArrayLike, distance: ArrayLike) -> FloatOrArray:
    r"""Escape velocity from radius :math:`r`, :math:`v_{esc} = \sqrt{2GM/r}`.

    Examples
    --------
    >>> round(escape_velocity(5.9722e24, 6.371e6))   # Earth's surface, m/s
    11186
    """
    return _core.escape_velocity(central_mass, distance)


def orbital_velocity(central_mass: ArrayLike, distance: ArrayLike) -> FloatOrArray:
    r"""Circular orbital speed at radius :math:`r`, :math:`v = \sqrt{GM/r}`.

    Examples
    --------
    >>> round(orbital_velocity(5.9722e24, 6.771e6))   # ~400 km altitude (ISS), m/s
    7673
    """
    return _core.orbital_velocity(central_mass, distance)


def orbital_period(semi_major_axis: ArrayLike, central_mass: ArrayLike) -> FloatOrArray:
    r"""Orbital period from Kepler's third law, :math:`T = 2\pi\sqrt{a^3 / GM}`.

    Parameters
    ----------
    semi_major_axis : float or array_like
        Semi-major axis :math:`a` [m]. Must be positive.
    central_mass : float or array_like
        Central mass :math:`M` [kg]. Must be positive.

    Examples
    --------
    >>> days = orbital_period(1.496e11, 1.989e30) / 86400   # Earth around the Sun
    >>> round(days)
    365
    """
    return _core.orbital_period(semi_major_axis, central_mass)


__all__ = [
    "escape_velocity",
    "gravitational_field_strength",
    "gravitational_force",
    "gravitational_force_vector",
    "gravitational_potential_energy",
    "orbital_period",
    "orbital_velocity",
]

"""Rigid-body rotation: torque, angular momentum, moments of inertia, circular motion.

Units are SI; angles are in radians and angular velocities in rad/s.
"""

from __future__ import annotations

from physicsforge import _physicsforge as _core
from physicsforge._typing import ArrayLike, FloatOrArray, VectorLike
from physicsforge.vector import Vector3, as_vector


def torque(lever_arm: ArrayLike, force: ArrayLike, angle: ArrayLike) -> FloatOrArray:
    r"""Torque magnitude :math:`\tau = r F \sin\theta`.

    Parameters
    ----------
    lever_arm : float or array_like
        Distance from the pivot to the point of application :math:`r` [m].
    force : float or array_like
        Force magnitude :math:`F` [N].
    angle : float or array_like
        Angle between the lever arm and the force :math:`\theta` [rad].

    Examples
    --------
    >>> import math
    >>> torque(2.0, 10.0, math.pi / 2)
    20.0
    """
    return _core.torque(lever_arm, force, angle)


def torque_vector(lever_arm: VectorLike, force: VectorLike) -> Vector3:
    r"""Torque vector :math:`\boldsymbol{\tau} = \mathbf{r}\times\mathbf{F}`.

    Examples
    --------
    >>> torque_vector((1, 0, 0), (0, 5, 0))
    Vector3(0.0, 0.0, 5.0)
    """
    return _core.torque_vector(as_vector(lever_arm, "lever_arm"), as_vector(force, "force"))


def angular_acceleration(torque: ArrayLike, moment_of_inertia: ArrayLike) -> FloatOrArray:
    r"""Newton's second law for rotation, :math:`\alpha = \tau / I`.

    Examples
    --------
    >>> angular_acceleration(10.0, 4.0)
    2.5
    """
    return _core.angular_acceleration(torque, moment_of_inertia)


def angular_momentum(moment_of_inertia: ArrayLike, angular_velocity: ArrayLike) -> FloatOrArray:
    r"""Angular momentum about a fixed axis, :math:`L = I\omega`.

    Examples
    --------
    >>> angular_momentum(2.0, 3.0)
    6.0
    """
    return _core.angular_momentum(moment_of_inertia, angular_velocity)


def angular_momentum_vector(position: VectorLike, linear_momentum: VectorLike) -> Vector3:
    r"""Angular momentum of a particle, :math:`\mathbf{L} = \mathbf{r}\times\mathbf{p}`.

    Examples
    --------
    >>> angular_momentum_vector((0, 2, 0), (3, 0, 0))
    Vector3(0.0, 0.0, -6.0)
    """
    return _core.angular_momentum_vector(as_vector(position, "position"), as_vector(linear_momentum, "linear_momentum"))


def rotational_kinetic_energy(moment_of_inertia: ArrayLike, angular_velocity: ArrayLike) -> FloatOrArray:
    r"""Rotational kinetic energy :math:`K = \tfrac{1}{2} I \omega^2`.

    Examples
    --------
    >>> rotational_kinetic_energy(2.0, 3.0)
    9.0
    """
    return _core.rotational_kinetic_energy(moment_of_inertia, angular_velocity)


def tangential_velocity(angular_velocity: ArrayLike, radius: ArrayLike) -> FloatOrArray:
    r"""Tangential speed of a point at radius :math:`r`, :math:`v = \omega r`.

    Examples
    --------
    >>> tangential_velocity(4.0, 0.5)
    2.0
    """
    return _core.tangential_velocity(angular_velocity, radius)


def centripetal_acceleration(speed: ArrayLike, radius: ArrayLike) -> FloatOrArray:
    r"""Centripetal acceleration :math:`a_c = v^2 / r`.

    Examples
    --------
    >>> centripetal_acceleration(3.0, 1.5)
    6.0
    """
    return _core.centripetal_acceleration(speed, radius)


def centripetal_force(mass: ArrayLike, speed: ArrayLike, radius: ArrayLike) -> FloatOrArray:
    r"""Centripetal force :math:`F_c = m v^2 / r`.

    Examples
    --------
    >>> centripetal_force(2.0, 3.0, 1.5)
    12.0
    """
    return _core.centripetal_force(mass, speed, radius)


def moment_of_inertia_point_mass(mass: ArrayLike, radius: ArrayLike) -> FloatOrArray:
    r"""Point mass at distance :math:`r` from the axis: :math:`I = m r^2`."""
    return _core.moment_of_inertia_point_mass(mass, radius)


def moment_of_inertia_solid_sphere(mass: ArrayLike, radius: ArrayLike) -> FloatOrArray:
    r"""Solid sphere about a diameter: :math:`I = \tfrac{2}{5} m r^2`.

    Examples
    --------
    >>> moment_of_inertia_solid_sphere(5.0, 2.0)
    8.0
    """
    return _core.moment_of_inertia_solid_sphere(mass, radius)


def moment_of_inertia_hollow_sphere(mass: ArrayLike, radius: ArrayLike) -> FloatOrArray:
    r"""Thin spherical shell about a diameter: :math:`I = \tfrac{2}{3} m r^2`."""
    return _core.moment_of_inertia_hollow_sphere(mass, radius)


def moment_of_inertia_solid_cylinder(mass: ArrayLike, radius: ArrayLike) -> FloatOrArray:
    r"""Solid cylinder or disc about its symmetry axis: :math:`I = \tfrac{1}{2} m r^2`."""
    return _core.moment_of_inertia_solid_cylinder(mass, radius)


def moment_of_inertia_hoop(mass: ArrayLike, radius: ArrayLike) -> FloatOrArray:
    r"""Thin hoop or ring about its symmetry axis: :math:`I = m r^2`."""
    return _core.moment_of_inertia_hoop(mass, radius)


def moment_of_inertia_rod_center(mass: ArrayLike, length: ArrayLike) -> FloatOrArray:
    r"""Thin rod about its centre, perpendicular to its length: :math:`I = \tfrac{1}{12} m L^2`."""
    return _core.moment_of_inertia_rod_center(mass, length)


def moment_of_inertia_rod_end(mass: ArrayLike, length: ArrayLike) -> FloatOrArray:
    r"""Thin rod about one end, perpendicular to its length: :math:`I = \tfrac{1}{3} m L^2`."""
    return _core.moment_of_inertia_rod_end(mass, length)


def parallel_axis(moment_of_inertia_cm: ArrayLike, mass: ArrayLike, distance: ArrayLike) -> FloatOrArray:
    r"""Parallel-axis theorem, :math:`I = I_{cm} + m d^2`.

    Examples
    --------
    >>> # A rod about its end, derived from the centre-of-mass value.
    >>> parallel_axis(moment_of_inertia_rod_center(3.0, 2.0), 3.0, 1.0)
    4.0
    """
    return _core.parallel_axis(moment_of_inertia_cm, mass, distance)


__all__ = [
    "angular_acceleration",
    "angular_momentum",
    "angular_momentum_vector",
    "centripetal_acceleration",
    "centripetal_force",
    "moment_of_inertia_hollow_sphere",
    "moment_of_inertia_hoop",
    "moment_of_inertia_point_mass",
    "moment_of_inertia_rod_center",
    "moment_of_inertia_rod_end",
    "moment_of_inertia_solid_cylinder",
    "moment_of_inertia_solid_sphere",
    "parallel_axis",
    "rotational_kinetic_energy",
    "tangential_velocity",
    "torque",
    "torque_vector",
]

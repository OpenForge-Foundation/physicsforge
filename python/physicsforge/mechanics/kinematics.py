"""Constant-acceleration kinematics and ideal projectile motion.

All functions broadcast over NumPy arrays (computed in C++) and return a ``float``
for scalar inputs. Units are SI; angles are in radians.
"""

from __future__ import annotations

from physicsforge import _physicsforge as _core
from physicsforge._typing import ArrayLike, FloatOrArray
from physicsforge.constants import STANDARD_GRAVITY


def displacement(initial_velocity: ArrayLike, acceleration: ArrayLike, time: ArrayLike) -> FloatOrArray:
    r"""Displacement under constant acceleration.

    .. math:: s = v_0 t + \tfrac{1}{2} a t^2

    Parameters
    ----------
    initial_velocity : float or array_like
        Initial velocity :math:`v_0` [m/s].
    acceleration : float or array_like
        Constant acceleration :math:`a` [m/s^2].
    time : float or array_like
        Elapsed time :math:`t` [s]. Must be non-negative.

    Returns
    -------
    float or numpy.ndarray
        Displacement :math:`s` [m].

    Raises
    ------
    ValueError
        If ``time`` is negative or any input is not finite.

    Examples
    --------
    >>> displacement(initial_velocity=2.0, acceleration=3.0, time=4.0)
    32.0
    >>> displacement(0.0, 9.80665, [0.0, 1.0, 2.0])
    array([ 0.      ,  4.903325, 19.6133  ])
    """
    return _core.displacement(initial_velocity, acceleration, time)


def final_velocity(initial_velocity: ArrayLike, acceleration: ArrayLike, time: ArrayLike) -> FloatOrArray:
    r"""Velocity after time :math:`t` under constant acceleration, :math:`v = v_0 + a t`.

    Parameters
    ----------
    initial_velocity : float or array_like
        Initial velocity :math:`v_0` [m/s].
    acceleration : float or array_like
        Constant acceleration :math:`a` [m/s^2].
    time : float or array_like
        Elapsed time :math:`t` [s]. Must be non-negative.

    Returns
    -------
    float or numpy.ndarray
        Final velocity [m/s].

    Examples
    --------
    >>> final_velocity(2.0, 3.0, 4.0)
    14.0
    """
    return _core.final_velocity(initial_velocity, acceleration, time)


def final_velocity_from_displacement(
    initial_velocity: ArrayLike, acceleration: ArrayLike, displacement: ArrayLike
) -> FloatOrArray:
    r"""Speed after a displacement under constant acceleration (time-independent form).

    .. math:: |v| = \sqrt{v_0^2 + 2 a s}

    Parameters
    ----------
    initial_velocity : float or array_like
        Initial velocity :math:`v_0` [m/s].
    acceleration : float or array_like
        Constant acceleration :math:`a` [m/s^2].
    displacement : float or array_like
        Displacement :math:`s` [m].

    Returns
    -------
    float or numpy.ndarray
        Final speed :math:`|v|` [m/s].

    Raises
    ------
    ValueError
        If :math:`v_0^2 + 2as < 0`, i.e. the displacement is never reached.

    Examples
    --------
    >>> final_velocity_from_displacement(0.0, 2.0, 4.0)
    4.0
    """
    return _core.final_velocity_from_displacement(initial_velocity, acceleration, displacement)


def average_velocity(displacement: ArrayLike, time: ArrayLike) -> FloatOrArray:
    r"""Average velocity :math:`\bar v = \Delta x / \Delta t`.

    Parameters
    ----------
    displacement : float or array_like
        Net displacement :math:`\Delta x` [m].
    time : float or array_like
        Elapsed time :math:`\Delta t` [s]. Must be positive.

    Examples
    --------
    >>> average_velocity(100.0, 20.0)
    5.0
    """
    return _core.average_velocity(displacement, time)


def average_acceleration(velocity_change: ArrayLike, time: ArrayLike) -> FloatOrArray:
    r"""Average acceleration :math:`\bar a = \Delta v / \Delta t`.

    Parameters
    ----------
    velocity_change : float or array_like
        Change in velocity :math:`\Delta v` [m/s].
    time : float or array_like
        Elapsed time :math:`\Delta t` [s]. Must be positive.

    Examples
    --------
    >>> average_acceleration(27.0, 9.0)
    3.0
    """
    return _core.average_acceleration(velocity_change, time)


def projectile_range(
    initial_speed: ArrayLike, launch_angle: ArrayLike, gravity: ArrayLike = STANDARD_GRAVITY
) -> FloatOrArray:
    r"""Horizontal range of a drag-free projectile on level ground.

    .. math:: R = \frac{v_0^2 \sin 2\theta}{g}

    Parameters
    ----------
    initial_speed : float or array_like
        Launch speed :math:`v_0` [m/s]. Must be non-negative.
    launch_angle : float or array_like
        Elevation angle :math:`\theta` [rad].
    gravity : float or array_like, optional
        Gravitational acceleration :math:`g` [m/s^2]. Defaults to standard gravity.

    Examples
    --------
    >>> import math
    >>> projectile_range(10.0, math.pi / 4, gravity=10.0)
    10.0
    """
    return _core.projectile_range(initial_speed, launch_angle, gravity)


def projectile_time_of_flight(
    initial_speed: ArrayLike, launch_angle: ArrayLike, gravity: ArrayLike = STANDARD_GRAVITY
) -> FloatOrArray:
    r"""Time of flight of a drag-free projectile on level ground, :math:`T = 2 v_0 \sin\theta / g`.

    Examples
    --------
    >>> import math
    >>> projectile_time_of_flight(10.0, math.pi / 2, gravity=10.0)
    2.0
    """
    return _core.projectile_time_of_flight(initial_speed, launch_angle, gravity)


def projectile_max_height(
    initial_speed: ArrayLike, launch_angle: ArrayLike, gravity: ArrayLike = STANDARD_GRAVITY
) -> FloatOrArray:
    r"""Apex height of a drag-free projectile, :math:`H = (v_0 \sin\theta)^2 / 2g`.

    Examples
    --------
    >>> import math
    >>> projectile_max_height(10.0, math.pi / 2, gravity=10.0)
    5.0
    """
    return _core.projectile_max_height(initial_speed, launch_angle, gravity)


__all__ = [
    "average_acceleration",
    "average_velocity",
    "displacement",
    "final_velocity",
    "final_velocity_from_displacement",
    "projectile_max_height",
    "projectile_range",
    "projectile_time_of_flight",
]

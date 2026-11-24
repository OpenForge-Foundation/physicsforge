r"""Rocket propulsion: the Tsiolkovsky rocket equation and related quantities.

The ideal rocket equation relates the velocity change of a rocket to its
effective exhaust velocity and mass ratio, neglecting gravity and drag losses:

.. math:: \Delta v = v_e \ln\frac{m_0}{m_f}

Examples
--------
How much propellant does a 25 t (dry) upper stage need for 9.4 km/s with a
kerosene/LOX engine (:math:`I_{sp}` ≈ 311 s)?

>>> from physicsforge.mechanics import rocketry
>>> ve = rocketry.exhaust_velocity_from_isp(311.0)
>>> round(rocketry.rocket_propellant_mass(9400.0, ve, 25_000.0) / 1000, 1)   # tonnes
520.1
"""

from __future__ import annotations

from physicsforge import _physicsforge as _core
from physicsforge._typing import ArrayLike, FloatOrArray
from physicsforge.constants import STANDARD_GRAVITY


def rocket_delta_v(exhaust_velocity: ArrayLike, initial_mass: ArrayLike, final_mass: ArrayLike) -> FloatOrArray:
    r"""Tsiolkovsky rocket equation, :math:`\Delta v = v_e \ln(m_0 / m_f)`.

    Parameters
    ----------
    exhaust_velocity : float or array_like
        Effective exhaust velocity :math:`v_e` [m/s]. Must be non-negative.
    initial_mass : float or array_like
        Initial (wet) mass :math:`m_0` [kg]. Must satisfy :math:`m_0 \ge m_f`.
    final_mass : float or array_like
        Final (dry) mass :math:`m_f` [kg]. Must be positive.

    Returns
    -------
    float or numpy.ndarray
        Ideal velocity change :math:`\Delta v` [m/s].

    Raises
    ------
    ValueError
        If ``final_mass`` is not positive or exceeds ``initial_mass``.

    Examples
    --------
    >>> import math
    >>> rocket_delta_v(3000.0, math.e * 1000.0, 1000.0)
    3000.0
    """
    return _core.rocket_delta_v(exhaust_velocity, initial_mass, final_mass)


def rocket_mass_ratio(delta_v: ArrayLike, exhaust_velocity: ArrayLike) -> FloatOrArray:
    r"""Mass ratio required for a velocity change, :math:`m_0/m_f = e^{\Delta v / v_e}`.

    Examples
    --------
    >>> round(rocket_mass_ratio(3000.0, 3000.0), 6)
    2.718282
    """
    return _core.rocket_mass_ratio(delta_v, exhaust_velocity)


def rocket_propellant_mass(delta_v: ArrayLike, exhaust_velocity: ArrayLike, final_mass: ArrayLike) -> FloatOrArray:
    r"""Propellant mass needed for :math:`\Delta v`, :math:`m_p = m_f (e^{\Delta v / v_e} - 1)`.

    Computed with ``expm1`` so it stays accurate for small manoeuvres.

    Examples
    --------
    >>> round(rocket_propellant_mass(1.0, 3000.0, 1000.0), 6)   # tiny burn, no cancellation
    0.333389
    """
    return _core.rocket_propellant_mass(delta_v, exhaust_velocity, final_mass)


def exhaust_velocity_from_isp(
    specific_impulse: ArrayLike, standard_gravity: ArrayLike = STANDARD_GRAVITY
) -> FloatOrArray:
    r"""Effective exhaust velocity from specific impulse, :math:`v_e = I_{sp} g_0`.

    Examples
    --------
    >>> exhaust_velocity_from_isp(300.0)
    2941.995
    """
    return _core.exhaust_velocity_from_isp(specific_impulse, standard_gravity)


def specific_impulse(exhaust_velocity: ArrayLike, standard_gravity: ArrayLike = STANDARD_GRAVITY) -> FloatOrArray:
    r"""Specific impulse from effective exhaust velocity, :math:`I_{sp} = v_e / g_0` [s].

    Examples
    --------
    >>> specific_impulse(2941.995)
    300.0
    """
    return _core.specific_impulse(exhaust_velocity, standard_gravity)


def rocket_thrust(
    mass_flow_rate: ArrayLike,
    exhaust_velocity: ArrayLike,
    exit_pressure: ArrayLike = 0.0,
    ambient_pressure: ArrayLike = 0.0,
    exit_area: ArrayLike = 0.0,
) -> FloatOrArray:
    r"""Rocket thrust including the pressure term, :math:`F = \dot m v_e + (p_e - p_a) A_e`.

    Parameters
    ----------
    mass_flow_rate : float or array_like
        Propellant mass flow rate :math:`\dot m` [kg/s].
    exhaust_velocity : float or array_like
        Exhaust velocity :math:`v_e` [m/s].
    exit_pressure, ambient_pressure : float or array_like, optional
        Nozzle exit and ambient pressures [Pa]. Default 0 (vacuum, matched nozzle).
    exit_area : float or array_like, optional
        Nozzle exit area :math:`A_e` [m^2]. Default 0.

    Examples
    --------
    >>> rocket_thrust(100.0, 3000.0)
    300000.0
    """
    return _core.rocket_thrust(mass_flow_rate, exhaust_velocity, exit_pressure, ambient_pressure, exit_area)


def burn_time(propellant_mass: ArrayLike, mass_flow_rate: ArrayLike) -> FloatOrArray:
    r"""Burn time at constant mass flow, :math:`t = m_p / \dot m`.

    Examples
    --------
    >>> burn_time(1000.0, 50.0)
    20.0
    """
    return _core.burn_time(propellant_mass, mass_flow_rate)


__all__ = [
    "burn_time",
    "exhaust_velocity_from_isp",
    "rocket_delta_v",
    "rocket_mass_ratio",
    "rocket_propellant_mass",
    "rocket_thrust",
    "specific_impulse",
]

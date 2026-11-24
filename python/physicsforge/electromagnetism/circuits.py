"""Lumped-element circuit analysis: Ohm's law, dividers, networks and RC transients.

Scalar laws broadcast over NumPy arrays (computed in C++) and return a ``float``
for scalar inputs. Units: V, A, Ω, F, W, s.
"""

from __future__ import annotations

from physicsforge import _physicsforge as _core
from physicsforge._typing import ArrayLike, FloatOrArray


def ohms_law_voltage(current: ArrayLike, resistance: ArrayLike) -> FloatOrArray:
    r"""Ohm's law solved for voltage, :math:`V = I R`.

    Parameters
    ----------
    current : float or array_like
        Current :math:`I` [A].
    resistance : float or array_like
        Resistance :math:`R` [Ω]. Must be non-negative.

    Returns
    -------
    float or numpy.ndarray
        Voltage [V].

    Examples
    --------
    >>> ohms_law_voltage(current=2.0, resistance=5.0)
    10.0
    """
    return _core.ohms_law_voltage(current, resistance)


def ohms_law_current(voltage: ArrayLike, resistance: ArrayLike) -> FloatOrArray:
    r"""Ohm's law solved for current, :math:`I = V / R`.

    Raises
    ------
    ValueError
        If ``resistance`` is not positive.

    Examples
    --------
    >>> ohms_law_current(12.0, 4.0)
    3.0
    """
    return _core.ohms_law_current(voltage, resistance)


def ohms_law_resistance(voltage: ArrayLike, current: ArrayLike) -> FloatOrArray:
    r"""Ohm's law solved for resistance, :math:`R = V / I`.

    Raises
    ------
    ValueError
        If ``current`` is zero.

    Examples
    --------
    >>> ohms_law_resistance(12.0, 3.0)
    4.0
    """
    return _core.ohms_law_resistance(voltage, current)


def electrical_power(voltage: ArrayLike, current: ArrayLike) -> FloatOrArray:
    r"""Electrical power :math:`P = V I`.

    Examples
    --------
    >>> electrical_power(230.0, 2.0)
    460.0
    """
    return _core.electrical_power(voltage, current)


def power_from_current(current: ArrayLike, resistance: ArrayLike) -> FloatOrArray:
    r"""Joule heating :math:`P = I^2 R`.

    Examples
    --------
    >>> power_from_current(2.0, 5.0)
    20.0
    """
    return _core.power_from_current(current, resistance)


def power_from_voltage(voltage: ArrayLike, resistance: ArrayLike) -> FloatOrArray:
    r"""Joule heating :math:`P = V^2 / R`.

    Examples
    --------
    >>> power_from_voltage(10.0, 5.0)
    20.0
    """
    return _core.power_from_voltage(voltage, resistance)


def voltage_divider(input_voltage: ArrayLike, r1: ArrayLike, r2: ArrayLike) -> FloatOrArray:
    r"""Output of an unloaded two-resistor voltage divider, taken across ``r2``.

    .. math:: V_{out} = V_{in} \frac{R_2}{R_1 + R_2}

    Parameters
    ----------
    input_voltage : float or array_like
        Source voltage :math:`V_{in}` [V].
    r1 : float or array_like
        Upper (series) resistor [Ω].
    r2 : float or array_like
        Lower resistor, across which the output is measured [Ω].

    Raises
    ------
    ValueError
        If a resistance is negative or both are zero.

    Examples
    --------
    >>> voltage_divider(12.0, r1=1_000.0, r2=2_000.0)
    8.0
    >>> voltage_divider(5.0, 10_000.0, [10_000.0, 30_000.0])
    array([2.5 , 3.75])
    """
    return _core.voltage_divider(input_voltage, r1, r2)


def current_divider(total_current: ArrayLike, r_branch: ArrayLike, r_other: ArrayLike) -> FloatOrArray:
    r"""Current through one branch of a two-branch parallel divider.

    .. math:: I_{branch} = I_{total} \frac{R_{other}}{R_{branch} + R_{other}}

    Examples
    --------
    >>> current_divider(3.0, r_branch=10.0, r_other=20.0)
    2.0
    """
    return _core.current_divider(total_current, r_branch, r_other)


def series_resistance(resistances: ArrayLike) -> float:
    r"""Equivalent series resistance :math:`R = \sum_i R_i`.

    Examples
    --------
    >>> series_resistance([100.0, 220.0, 330.0])
    650.0
    """
    return _core.series_resistance(resistances)


def parallel_resistance(resistances: ArrayLike) -> float:
    r"""Equivalent parallel resistance :math:`1/R = \sum_i 1/R_i`.

    Raises
    ------
    ValueError
        If the list is empty or any resistance is not positive.

    Examples
    --------
    >>> parallel_resistance([100.0, 100.0])
    50.0
    """
    return _core.parallel_resistance(resistances)


def series_capacitance(capacitances: ArrayLike) -> float:
    r"""Equivalent series capacitance :math:`1/C = \sum_i 1/C_i`.

    Examples
    --------
    >>> series_capacitance([2e-6, 2e-6])
    1e-06
    """
    return _core.series_capacitance(capacitances)


def parallel_capacitance(capacitances: ArrayLike) -> float:
    r"""Equivalent parallel capacitance :math:`C = \sum_i C_i`.

    Examples
    --------
    >>> parallel_capacitance([1e-6, 2e-6])
    3e-06
    """
    return _core.parallel_capacitance(capacitances)


def resistance_from_resistivity(resistivity: ArrayLike, length: ArrayLike, area: ArrayLike) -> FloatOrArray:
    r"""Resistance of a uniform conductor, :math:`R = \rho L / A`.

    Examples
    --------
    >>> round(resistance_from_resistivity(1.68e-8, 100.0, 1e-6), 6)   # 100 m of 1 mm² copper
    1.68
    """
    return _core.resistance_from_resistivity(resistivity, length, area)


def rc_time_constant(resistance: ArrayLike, capacitance: ArrayLike) -> FloatOrArray:
    r"""RC time constant :math:`\tau = R C` [s].

    Examples
    --------
    >>> rc_time_constant(1_000.0, 1e-3)
    1.0
    """
    return _core.rc_time_constant(resistance, capacitance)


def capacitor_charging_voltage(
    source_voltage: ArrayLike, time: ArrayLike, resistance: ArrayLike, capacitance: ArrayLike
) -> FloatOrArray:
    r"""Capacitor voltage while charging from 0 V, :math:`V(t) = V_s (1 - e^{-t/RC})`.

    Examples
    --------
    >>> round(capacitor_charging_voltage(5.0, time=1.0, resistance=1_000.0, capacitance=1e-3), 4)
    3.1606
    """
    return _core.capacitor_charging_voltage(source_voltage, time, resistance, capacitance)


def capacitor_discharging_voltage(
    initial_voltage: ArrayLike, time: ArrayLike, resistance: ArrayLike, capacitance: ArrayLike
) -> FloatOrArray:
    r"""Capacitor voltage while discharging, :math:`V(t) = V_0 e^{-t/RC}`.

    Examples
    --------
    >>> round(capacitor_discharging_voltage(5.0, 1.0, 1_000.0, 1e-3), 4)
    1.8394
    """
    return _core.capacitor_discharging_voltage(initial_voltage, time, resistance, capacitance)


__all__ = [
    "capacitor_charging_voltage",
    "capacitor_discharging_voltage",
    "current_divider",
    "electrical_power",
    "ohms_law_current",
    "ohms_law_resistance",
    "ohms_law_voltage",
    "parallel_capacitance",
    "parallel_resistance",
    "power_from_current",
    "power_from_voltage",
    "rc_time_constant",
    "resistance_from_resistivity",
    "series_capacitance",
    "series_resistance",
    "voltage_divider",
]

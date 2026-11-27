"""Heat transfer and the laws of thermodynamics.

Covers calorimetry, conduction (Fourier's law), convection (Newton's law of
cooling), radiation (Stefan-Boltzmann law), the first law, Carnot limits and the
ideal gas law. Temperatures are absolute (K) except where a function works with
temperature *differences*, which are the same in K and °C.

All functions broadcast over NumPy arrays and return ``float`` for scalar inputs.
"""

from __future__ import annotations

from physicsforge import _physicsforge as _core
from physicsforge._typing import ArrayLike, FloatOrArray

# --- Calorimetry ---------------------------------------------------------------------------------


def sensible_heat(mass: ArrayLike, specific_heat: ArrayLike, temperature_change: ArrayLike) -> FloatOrArray:
    r"""Heat to change a body's temperature without phase change, :math:`Q = m c \Delta T`.

    Parameters
    ----------
    mass : float or array_like
        Mass [kg].
    specific_heat : float or array_like
        Specific heat capacity :math:`c` [J kg^-1 K^-1].
    temperature_change : float or array_like
        :math:`\Delta T` [K].

    Examples
    --------
    >>> sensible_heat(1.0, 4186.0, 10.0)   # warm 1 kg of water by 10 K
    41860.0
    """
    return _core.sensible_heat(mass, specific_heat, temperature_change)


def latent_heat(mass: ArrayLike, specific_latent_heat: ArrayLike) -> FloatOrArray:
    r"""Heat absorbed or released in a phase change, :math:`Q = m L`.

    Examples
    --------
    >>> latent_heat(2.0, 334_000.0)   # melt 2 kg of ice
    668000.0
    """
    return _core.latent_heat(mass, specific_latent_heat)


# --- Conduction ----------------------------------------------------------------------------------


def fourier_heat_flux(thermal_conductivity: ArrayLike, temperature_gradient: ArrayLike) -> FloatOrArray:
    r"""Fourier's law of heat conduction in one dimension, :math:`q = -k\, dT/dx`.

    The minus sign expresses that heat flows from hot to cold, i.e. down the
    temperature gradient.

    Parameters
    ----------
    thermal_conductivity : float or array_like
        :math:`k` [W m^-1 K^-1].
    temperature_gradient : float or array_like
        :math:`dT/dx` [K/m].

    Returns
    -------
    float or numpy.ndarray
        Heat flux :math:`q` [W/m^2].

    Examples
    --------
    >>> fourier_heat_flux(200.0, -50.0)
    10000.0
    """
    return _core.fourier_heat_flux(thermal_conductivity, temperature_gradient)


def conduction_rate(
    thermal_conductivity: ArrayLike, area: ArrayLike, temperature_difference: ArrayLike, thickness: ArrayLike
) -> FloatOrArray:
    r"""Steady heat conduction through a plane wall, :math:`\dot Q = k A \Delta T / L`.

    Parameters
    ----------
    thermal_conductivity : float or array_like
        :math:`k` [W m^-1 K^-1].
    area : float or array_like
        Wall area :math:`A` [m^2].
    temperature_difference : float or array_like
        :math:`T_{hot} - T_{cold}` [K].
    thickness : float or array_like
        Wall thickness :math:`L` [m]. Must be positive.

    Returns
    -------
    float or numpy.ndarray
        Heat flow rate [W].

    Examples
    --------
    >>> conduction_rate(0.8, area=10.0, temperature_difference=20.0, thickness=0.2)   # brick wall
    800.0
    """
    return _core.conduction_rate(thermal_conductivity, area, temperature_difference, thickness)


def conduction_resistance(thickness: ArrayLike, thermal_conductivity: ArrayLike, area: ArrayLike) -> FloatOrArray:
    r"""Conductive thermal resistance of a plane wall, :math:`R = L / kA` [K/W].

    Resistances in series add, so multilayer walls can be analysed like circuits
    with :func:`heat_rate_from_resistance`.
    """
    return _core.conduction_resistance(thickness, thermal_conductivity, area)


def convection_resistance(heat_transfer_coefficient: ArrayLike, area: ArrayLike) -> FloatOrArray:
    r"""Convective thermal resistance, :math:`R = 1 / hA` [K/W].

    Examples
    --------
    >>> convection_resistance(25.0, 2.0)
    0.02
    """
    return _core.convection_resistance(heat_transfer_coefficient, area)


def heat_rate_from_resistance(temperature_difference: ArrayLike, thermal_resistance: ArrayLike) -> FloatOrArray:
    r"""Heat flow across a thermal resistance network, :math:`\dot Q = \Delta T / R_{total}`.

    Examples
    --------
    >>> wall = conduction_resistance(0.2, 0.8, 10.0) + convection_resistance(25.0, 10.0)
    >>> round(heat_rate_from_resistance(20.0, wall), 3)
    689.655
    """
    return _core.heat_rate_from_resistance(temperature_difference, thermal_resistance)


def thermal_diffusivity(thermal_conductivity: ArrayLike, density: ArrayLike, specific_heat: ArrayLike) -> FloatOrArray:
    r"""Thermal diffusivity :math:`\alpha = k / \rho c_p` [m^2/s], the coefficient in the heat equation.

    Examples
    --------
    >>> round(thermal_diffusivity(237.0, 2700.0, 897.0), 8)   # aluminium
    9.786e-05
    """
    return _core.thermal_diffusivity(thermal_conductivity, density, specific_heat)


# --- Convection ----------------------------------------------------------------------------------


def newton_cooling_rate(
    heat_transfer_coefficient: ArrayLike,
    area: ArrayLike,
    surface_temperature: ArrayLike,
    ambient_temperature: ArrayLike,
) -> FloatOrArray:
    r"""Newton's law of cooling (convective heat rate), :math:`\dot Q = h A (T_s - T_\infty)`.

    Examples
    --------
    >>> newton_cooling_rate(10.0, 2.0, surface_temperature=80.0, ambient_temperature=20.0)
    1200.0
    """
    return _core.newton_cooling_rate(heat_transfer_coefficient, area, surface_temperature, ambient_temperature)


def newton_cooling_temperature(
    initial_temperature: ArrayLike, ambient_temperature: ArrayLike, cooling_constant: ArrayLike, time: ArrayLike
) -> FloatOrArray:
    r"""Lumped-capacitance temperature history, :math:`T(t) = T_{env} + (T_0 - T_{env}) e^{-kt}`.

    Parameters
    ----------
    initial_temperature : float or array_like
        :math:`T_0`.
    ambient_temperature : float or array_like
        :math:`T_{env}`.
    cooling_constant : float or array_like
        :math:`k = hA / (m c)` [1/s].
    time : float or array_like
        Elapsed time [s].

    Examples
    --------
    >>> round(newton_cooling_temperature(90.0, 20.0, 0.1, time=10.0), 3)   # coffee cooling
    45.752
    """
    return _core.newton_cooling_temperature(initial_temperature, ambient_temperature, cooling_constant, time)


# --- Radiation -----------------------------------------------------------------------------------


def radiated_power(emissivity: ArrayLike, area: ArrayLike, temperature: ArrayLike) -> FloatOrArray:
    r"""Stefan-Boltzmann law for a grey body, :math:`P = \varepsilon \sigma A T^4`.

    Parameters
    ----------
    emissivity : float or array_like
        Emissivity in [0, 1] (1 for a black body).
    area : float or array_like
        Radiating area [m^2].
    temperature : float or array_like
        Absolute temperature [K].

    Examples
    --------
    >>> round(radiated_power(1.0, 1.0, 1000.0))
    56704
    """
    return _core.radiated_power(emissivity, area, temperature)


def net_radiated_power(
    emissivity: ArrayLike, area: ArrayLike, temperature: ArrayLike, surroundings_temperature: ArrayLike
) -> FloatOrArray:
    r"""Net radiative exchange with large surroundings, :math:`P = \varepsilon\sigma A (T^4 - T_{surr}^4)`."""
    return _core.net_radiated_power(emissivity, area, temperature, surroundings_temperature)


# --- Laws of thermodynamics ----------------------------------------------------------------------


def internal_energy_change(heat_added: ArrayLike, work_done_by_system: ArrayLike) -> FloatOrArray:
    r"""First law of thermodynamics, :math:`\Delta U = Q - W` (W done *by* the system).

    Examples
    --------
    >>> internal_energy_change(100.0, 40.0)
    60.0
    """
    return _core.internal_energy_change(heat_added, work_done_by_system)


def carnot_efficiency(cold_temperature: ArrayLike, hot_temperature: ArrayLike) -> FloatOrArray:
    r"""Maximum (Carnot) heat-engine efficiency, :math:`\eta = 1 - T_c / T_h`.

    Raises
    ------
    ValueError
        If a temperature is not positive or ``cold_temperature > hot_temperature``.

    Examples
    --------
    >>> carnot_efficiency(300.0, 600.0)
    0.5
    """
    return _core.carnot_efficiency(cold_temperature, hot_temperature)


def carnot_cop_refrigerator(cold_temperature: ArrayLike, hot_temperature: ArrayLike) -> FloatOrArray:
    r"""Maximum refrigerator coefficient of performance, :math:`COP = T_c / (T_h - T_c)`.

    Examples
    --------
    >>> carnot_cop_refrigerator(250.0, 300.0)
    5.0
    """
    return _core.carnot_cop_refrigerator(cold_temperature, hot_temperature)


def ideal_gas_pressure(amount: ArrayLike, temperature: ArrayLike, volume: ArrayLike) -> FloatOrArray:
    r"""Ideal gas law solved for pressure, :math:`P = nRT/V` [Pa].

    Examples
    --------
    >>> round(ideal_gas_pressure(1.0, 273.15, 0.022414))
    101325
    """
    return _core.ideal_gas_pressure(amount, temperature, volume)


def ideal_gas_volume(amount: ArrayLike, temperature: ArrayLike, pressure: ArrayLike) -> FloatOrArray:
    r"""Ideal gas law solved for volume, :math:`V = nRT/P` [m^3]."""
    return _core.ideal_gas_volume(amount, temperature, pressure)


def linear_thermal_expansion(
    expansion_coefficient: ArrayLike, initial_length: ArrayLike, temperature_change: ArrayLike
) -> FloatOrArray:
    r"""Linear thermal expansion, :math:`\Delta L = \alpha L_0 \Delta T` [m]."""
    return _core.linear_thermal_expansion(expansion_coefficient, initial_length, temperature_change)


def celsius_to_kelvin(celsius: ArrayLike) -> FloatOrArray:
    """Convert °C to K (exact offset 273.15). Raises ``ValueError`` below absolute zero.

    Examples
    --------
    >>> celsius_to_kelvin(25.0)
    298.15
    """
    return _core.celsius_to_kelvin(celsius)


def kelvin_to_celsius(kelvin: ArrayLike) -> FloatOrArray:
    """Convert K to °C (exact offset 273.15). Raises ``ValueError`` for negative kelvin."""
    return _core.kelvin_to_celsius(kelvin)


__all__ = [
    "carnot_cop_refrigerator",
    "carnot_efficiency",
    "celsius_to_kelvin",
    "conduction_rate",
    "conduction_resistance",
    "convection_resistance",
    "fourier_heat_flux",
    "heat_rate_from_resistance",
    "ideal_gas_pressure",
    "ideal_gas_volume",
    "internal_energy_change",
    "kelvin_to_celsius",
    "latent_heat",
    "linear_thermal_expansion",
    "net_radiated_power",
    "newton_cooling_rate",
    "newton_cooling_temperature",
    "radiated_power",
    "sensible_heat",
    "thermal_diffusivity",
]

r"""Maxwell's equations and their standard closed-form consequences.

Integral form (implemented here for the symmetric configurations where they give
closed-form answers):

.. math::

    \oint \mathbf{E}\cdot d\mathbf{A} &= Q_{enc} / \varepsilon_0          && \text{(Gauss)} \\
    \oint \mathbf{B}\cdot d\mathbf{A} &= 0                                && \text{(Gauss, magnetic)} \\
    \oint \mathbf{E}\cdot d\boldsymbol{\ell} &= -\frac{d\Phi_B}{dt}         && \text{(Faraday)} \\
    \oint \mathbf{B}\cdot d\boldsymbol{\ell} &= \mu_0 \left(I_{enc} + \varepsilon_0 \frac{d\Phi_E}{dt}\right)
        && \text{(Ampère-Maxwell)}

The differential forms (:math:`\nabla\cdot\mathbf{E} = \rho/\varepsilon_0`,
:math:`\nabla\times\mathbf{E} = -\partial_t\mathbf{B}`, ...) can be checked
numerically for any field with :func:`physicsforge.numerics.divergence` and
:func:`physicsforge.numerics.curl`.

Examples
--------
Gauss's law reproduces Coulomb's law: the flux of a point charge's field through
any enclosing sphere is :math:`q/\varepsilon_0`.

>>> import math
>>> from physicsforge.electromagnetism import electric_field_point_charge, gauss_electric_flux
>>> q, r = 1e-9, 0.25
>>> flux_through_sphere = electric_field_point_charge(q, r) * 4 * math.pi * r**2
>>> math.isclose(flux_through_sphere, gauss_electric_flux(q))
True
"""

from __future__ import annotations

from physicsforge import _physicsforge as _core
from physicsforge._typing import ArrayLike, FloatOrArray, VectorLike
from physicsforge.vector import Vector3, as_vector


def gauss_electric_flux(enclosed_charge: ArrayLike) -> FloatOrArray:
    r"""Gauss's law: electric flux through a closed surface, :math:`\Phi_E = Q_{enc}/\varepsilon_0`.

    Parameters
    ----------
    enclosed_charge : float or array_like
        Net charge enclosed by the surface [C].

    Returns
    -------
    float or numpy.ndarray
        Electric flux [V m].

    Examples
    --------
    >>> round(gauss_electric_flux(1e-9), 3)
    112.941
    """
    return _core.gauss_electric_flux(enclosed_charge)


def gauss_enclosed_charge(electric_flux: ArrayLike) -> FloatOrArray:
    r"""Gauss's law inverted: enclosed charge from flux, :math:`Q_{enc} = \varepsilon_0 \Phi_E`."""
    return _core.gauss_enclosed_charge(electric_flux)


def magnetic_flux(field_strength: ArrayLike, area: ArrayLike, angle: ArrayLike = 0.0) -> FloatOrArray:
    r"""Magnetic flux through a flat surface in a uniform field, :math:`\Phi_B = B A \cos\theta`.

    Parameters
    ----------
    field_strength : float or array_like
        Magnetic flux density :math:`B` [T].
    area : float or array_like
        Surface area :math:`A` [m^2].
    angle : float or array_like, optional
        Angle between :math:`\mathbf{B}` and the surface normal [rad]. Default 0.

    Examples
    --------
    >>> magnetic_flux(0.5, 0.2)
    0.1
    """
    return _core.magnetic_flux(field_strength, area, angle)


def faraday_emf(flux_rate_of_change: ArrayLike, turns: ArrayLike = 1.0) -> FloatOrArray:
    r"""Faraday's law of induction with Lenz's sign, :math:`\mathcal{E} = -N\, d\Phi_B/dt`.

    Parameters
    ----------
    flux_rate_of_change : float or array_like
        :math:`d\Phi_B/dt` through one turn [Wb/s].
    turns : float or array_like, optional
        Number of turns :math:`N \ge 1`. Default 1.

    Examples
    --------
    >>> faraday_emf(0.5, turns=100)
    -50.0
    """
    return _core.faraday_emf(flux_rate_of_change, turns)


def displacement_current(electric_flux_rate_of_change: ArrayLike) -> FloatOrArray:
    r"""Maxwell's displacement current, :math:`I_d = \varepsilon_0\, d\Phi_E/dt` [A]."""
    return _core.displacement_current(electric_flux_rate_of_change)


def ampere_maxwell_circulation(
    enclosed_current: ArrayLike, electric_flux_rate_of_change: ArrayLike = 0.0
) -> FloatOrArray:
    r"""Ampère-Maxwell law, :math:`\oint \mathbf{B}\cdot d\boldsymbol\ell = \mu_0 (I + \varepsilon_0\, d\Phi_E/dt)`.

    Parameters
    ----------
    enclosed_current : float or array_like
        Conduction current through the loop [A].
    electric_flux_rate_of_change : float or array_like, optional
        :math:`d\Phi_E/dt` through the loop [V m/s]. Default 0 (magnetostatics).

    Returns
    -------
    float or numpy.ndarray
        Circulation of :math:`\mathbf{B}` [T m].
    """
    return _core.ampere_maxwell_circulation(enclosed_current, electric_flux_rate_of_change)


def magnetic_field_straight_wire(current: ArrayLike, distance: ArrayLike) -> FloatOrArray:
    r"""Field of a long straight wire (Ampère's law), :math:`B = \mu_0 I / 2\pi r`.

    Examples
    --------
    >>> round(magnetic_field_straight_wire(10.0, 0.05) * 1e6, 3)   # microtesla
    40.0
    """
    return _core.magnetic_field_straight_wire(current, distance)


def magnetic_field_solenoid(turns_per_length: ArrayLike, current: ArrayLike) -> FloatOrArray:
    r"""Field inside an ideal long solenoid, :math:`B = \mu_0 n I` (``n`` in turns per metre)."""
    return _core.magnetic_field_solenoid(turns_per_length, current)


def magnetic_field_loop_center(current: ArrayLike, radius: ArrayLike) -> FloatOrArray:
    r"""Field at the centre of a circular current loop, :math:`B = \mu_0 I / 2R`."""
    return _core.magnetic_field_loop_center(current, radius)


def lorentz_force(
    charge: float,
    electric_field: VectorLike = (0.0, 0.0, 0.0),
    velocity: VectorLike = (0.0, 0.0, 0.0),
    magnetic_field: VectorLike = (0.0, 0.0, 0.0),
) -> Vector3:
    r"""Lorentz force on a point charge, :math:`\mathbf{F} = q(\mathbf{E} + \mathbf{v}\times\mathbf{B})`.

    Examples
    --------
    >>> lorentz_force(2.0, velocity=(3, 0, 0), magnetic_field=(0, 5, 0))
    Vector3(0.0, 0.0, 30.0)
    """
    return _core.lorentz_force(
        charge,
        as_vector(electric_field, "electric_field"),
        as_vector(velocity, "velocity"),
        as_vector(magnetic_field, "magnetic_field"),
    )


def em_wave_speed(relative_permittivity: ArrayLike, relative_permeability: ArrayLike = 1.0) -> FloatOrArray:
    r"""Phase velocity of light in a linear medium, :math:`v = c / \sqrt{\varepsilon_r \mu_r}`.

    Examples
    --------
    >>> em_wave_speed(4.0)
    149896229.0
    """
    return _core.em_wave_speed(relative_permittivity, relative_permeability)


def speed_of_light_from_vacuum_constants() -> float:
    r"""Speed of light predicted by Maxwell's equations, :math:`c = 1/\sqrt{\mu_0\varepsilon_0}`.

    Examples
    --------
    >>> round(speed_of_light_from_vacuum_constants())
    299792458
    """
    return _core.speed_of_light_from_vacuum_constants()


def em_wave_intensity(electric_field_amplitude: ArrayLike) -> FloatOrArray:
    r"""Time-averaged plane-wave intensity in vacuum, :math:`I = \tfrac{1}{2} c \varepsilon_0 E_0^2` [W/m^2]."""
    return _core.em_wave_intensity(electric_field_amplitude)


__all__ = [
    "ampere_maxwell_circulation",
    "displacement_current",
    "em_wave_intensity",
    "em_wave_speed",
    "faraday_emf",
    "gauss_electric_flux",
    "gauss_enclosed_charge",
    "lorentz_force",
    "magnetic_field_loop_center",
    "magnetic_field_solenoid",
    "magnetic_field_straight_wire",
    "magnetic_flux",
    "speed_of_light_from_vacuum_constants",
]

"""Electrostatics: Coulomb's law, point-charge fields and potentials, capacitors.

Units: C, m, N, V, F. Scalar laws broadcast over NumPy arrays.
"""

from __future__ import annotations

from physicsforge import _physicsforge as _core
from physicsforge._typing import ArrayLike, FloatOrArray, VectorLike
from physicsforge.vector import Vector3, as_vector


def coulomb_force(charge1: ArrayLike, charge2: ArrayLike, distance: ArrayLike) -> FloatOrArray:
    r"""Coulomb's law (signed), :math:`F = k_e q_1 q_2 / r^2`.

    Positive values are repulsive, negative values attractive.

    Parameters
    ----------
    charge1, charge2 : float or array_like
        Charges [C].
    distance : float or array_like
        Separation [m]. Must be positive.

    Examples
    --------
    >>> round(coulomb_force(1e-6, 1e-6, 0.1), 4)
    0.8988
    >>> coulomb_force(1e-6, -1e-6, 0.1) < 0   # opposite charges attract
    True
    """
    return _core.coulomb_force(charge1, charge2, distance)


def coulomb_force_vector(charge1: float, charge2: float, position1: VectorLike, position2: VectorLike) -> Vector3:
    r"""Electrostatic force on charge 1 due to charge 2.

    .. math:: \mathbf{F}_{12} = k_e q_1 q_2 \frac{\mathbf{r}_1 - \mathbf{r}_2}{|\mathbf{r}_1 - \mathbf{r}_2|^3}

    Examples
    --------
    >>> coulomb_force_vector(1e-6, 1e-6, (1, 0, 0), (0, 0, 0)).x > 0   # pushed away
    True
    """
    return _core.coulomb_force_vector(
        charge1, charge2, as_vector(position1, "position1"), as_vector(position2, "position2")
    )


def electric_field_point_charge(charge: ArrayLike, distance: ArrayLike) -> FloatOrArray:
    r"""Radial electric field of a point charge (signed), :math:`E = k_e q / r^2` [V/m].

    Examples
    --------
    >>> round(electric_field_point_charge(1e-9, 1.0), 4)
    8.9876
    """
    return _core.electric_field_point_charge(charge, distance)


def electric_field_vector(charge: float, source: VectorLike, point: VectorLike) -> Vector3:
    r"""Electric field vector at ``point`` due to a point charge at ``source``.

    Examples
    --------
    >>> e = electric_field_vector(1e-9, (0, 0, 0), (0, 0, 2))
    >>> e.x == 0.0 and e.z > 0
    True
    """
    return _core.electric_field_vector(charge, as_vector(source, "source"), as_vector(point, "point"))


def electric_potential_point_charge(charge: ArrayLike, distance: ArrayLike) -> FloatOrArray:
    r"""Electric potential of a point charge (zero at infinity), :math:`V = k_e q / r`."""
    return _core.electric_potential_point_charge(charge, distance)


def electric_potential_energy(charge1: ArrayLike, charge2: ArrayLike, distance: ArrayLike) -> FloatOrArray:
    r"""Electric potential energy of two point charges, :math:`U = k_e q_1 q_2 / r`."""
    return _core.electric_potential_energy(charge1, charge2, distance)


def parallel_plate_capacitance(
    relative_permittivity: ArrayLike, area: ArrayLike, separation: ArrayLike
) -> FloatOrArray:
    r"""Ideal parallel-plate capacitance, :math:`C = \varepsilon_r \varepsilon_0 A / d`.

    Parameters
    ----------
    relative_permittivity : float or array_like
        Dielectric constant :math:`\varepsilon_r \ge 1` (1 for vacuum).
    area : float or array_like
        Plate area [m^2].
    separation : float or array_like
        Plate separation [m]. Must be positive.

    Examples
    --------
    >>> round(parallel_plate_capacitance(1.0, 1e-2, 1e-3) * 1e12, 2)   # pF
    88.54
    """
    return _core.parallel_plate_capacitance(relative_permittivity, area, separation)


def capacitor_energy(capacitance: ArrayLike, voltage: ArrayLike) -> FloatOrArray:
    r"""Energy stored in a capacitor, :math:`U = \tfrac{1}{2} C V^2`.

    Examples
    --------
    >>> round(capacitor_energy(2e-6, 10.0), 12)
    0.0001
    """
    return _core.capacitor_energy(capacitance, voltage)


__all__ = [
    "capacitor_energy",
    "coulomb_force",
    "coulomb_force_vector",
    "electric_field_point_charge",
    "electric_field_vector",
    "electric_potential_energy",
    "electric_potential_point_charge",
    "parallel_plate_capacitance",
]

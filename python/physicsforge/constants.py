"""Fundamental physical and mathematical constants (SI units, CODATA 2022).

Values come from the compiled C++ core (``include/physicsforge/constants.h``), so
Python and C++ always agree exactly. Since the 2019 SI redefinition, ``c``, ``h``,
``e``, ``k_B`` and ``N_A`` are exact by definition.

Examples
--------
>>> from physicsforge import constants
>>> constants.SPEED_OF_LIGHT
299792458.0
>>> constants.PLANCK
6.62607015e-34
>>> round(constants.ELEMENTARY_CHARGE / constants.ELECTRON_VOLT)
1
"""

from __future__ import annotations

from typing import Final

from physicsforge import _physicsforge as _c

# Mathematical
PI: Final[float] = _c.PI
"""Archimedes' constant, :math:`\\pi`."""
TAU: Final[float] = _c.TAU
"""Full turn, :math:`2\\pi`."""
EULER: Final[float] = _c.EULER
"""Euler's number, :math:`e`."""

# Exact SI defining constants
SPEED_OF_LIGHT: Final[float] = _c.SPEED_OF_LIGHT
"""Speed of light in vacuum, :math:`c` [m s^-1] (exact)."""
PLANCK: Final[float] = _c.PLANCK
"""Planck constant, :math:`h` [J s] (exact)."""
REDUCED_PLANCK: Final[float] = _c.REDUCED_PLANCK
"""Reduced Planck constant, :math:`\\hbar = h / 2\\pi` [J s]."""
ELEMENTARY_CHARGE: Final[float] = _c.ELEMENTARY_CHARGE
"""Elementary charge, :math:`e` [C] (exact)."""
BOLTZMANN: Final[float] = _c.BOLTZMANN
"""Boltzmann constant, :math:`k_B` [J K^-1] (exact)."""
AVOGADRO: Final[float] = _c.AVOGADRO
"""Avogadro constant, :math:`N_A` [mol^-1] (exact)."""

# Derived exact constants
GAS_CONSTANT: Final[float] = _c.GAS_CONSTANT
"""Molar gas constant, :math:`R = N_A k_B` [J mol^-1 K^-1]."""
ELECTRON_VOLT: Final[float] = _c.ELECTRON_VOLT
"""Electron volt [J]."""
STEFAN_BOLTZMANN: Final[float] = _c.STEFAN_BOLTZMANN
"""Stefan-Boltzmann constant, :math:`\\sigma` [W m^-2 K^-4]."""
WIEN_DISPLACEMENT: Final[float] = _c.WIEN_DISPLACEMENT
"""Wien wavelength displacement constant, :math:`b` [m K]."""
STANDARD_GRAVITY: Final[float] = _c.STANDARD_GRAVITY
"""Standard acceleration of gravity, :math:`g_0` [m s^-2] (exact by convention)."""
STANDARD_ATMOSPHERE: Final[float] = _c.STANDARD_ATMOSPHERE
"""Standard atmosphere [Pa] (exact by convention)."""

# Measured constants (CODATA 2022)
GRAVITATIONAL_CONSTANT: Final[float] = _c.GRAVITATIONAL_CONSTANT
"""Newtonian constant of gravitation, :math:`G` [m^3 kg^-1 s^-2]."""
VACUUM_PERMITTIVITY: Final[float] = _c.VACUUM_PERMITTIVITY
"""Vacuum electric permittivity, :math:`\\varepsilon_0` [F m^-1]."""
VACUUM_PERMEABILITY: Final[float] = _c.VACUUM_PERMEABILITY
"""Vacuum magnetic permeability, :math:`\\mu_0` [N A^-2]."""
COULOMB_CONSTANT: Final[float] = _c.COULOMB_CONSTANT
"""Coulomb constant, :math:`k_e = 1 / 4\\pi\\varepsilon_0` [N m^2 C^-2]."""
ELECTRON_MASS: Final[float] = _c.ELECTRON_MASS
"""Electron rest mass, :math:`m_e` [kg]."""
PROTON_MASS: Final[float] = _c.PROTON_MASS
"""Proton rest mass, :math:`m_p` [kg]."""
NEUTRON_MASS: Final[float] = _c.NEUTRON_MASS
"""Neutron rest mass, :math:`m_n` [kg]."""

# Conventional short aliases
c: Final[float] = SPEED_OF_LIGHT
h: Final[float] = PLANCK
hbar: Final[float] = REDUCED_PLANCK
e: Final[float] = ELEMENTARY_CHARGE
k_B: Final[float] = BOLTZMANN
N_A: Final[float] = AVOGADRO
G: Final[float] = GRAVITATIONAL_CONSTANT
g0: Final[float] = STANDARD_GRAVITY
epsilon_0: Final[float] = VACUUM_PERMITTIVITY
mu_0: Final[float] = VACUUM_PERMEABILITY
pi: Final[float] = PI

__all__ = [
    "AVOGADRO",
    "BOLTZMANN",
    "COULOMB_CONSTANT",
    "ELECTRON_MASS",
    "ELECTRON_VOLT",
    "ELEMENTARY_CHARGE",
    "EULER",
    "GAS_CONSTANT",
    "GRAVITATIONAL_CONSTANT",
    "NEUTRON_MASS",
    "PI",
    "PLANCK",
    "PROTON_MASS",
    "REDUCED_PLANCK",
    "SPEED_OF_LIGHT",
    "STANDARD_ATMOSPHERE",
    "STANDARD_GRAVITY",
    "STEFAN_BOLTZMANN",
    "TAU",
    "VACUUM_PERMEABILITY",
    "VACUUM_PERMITTIVITY",
    "WIEN_DISPLACEMENT",
    "G",
    "N_A",
    "c",
    "e",
    "epsilon_0",
    "g0",
    "h",
    "hbar",
    "k_B",
    "mu_0",
    "pi",
]

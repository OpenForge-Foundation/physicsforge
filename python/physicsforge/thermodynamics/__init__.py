"""Thermodynamics and heat transfer.

Submodules
----------
heat
    Calorimetry, Fourier conduction, Newton cooling, Stefan-Boltzmann radiation,
    the first law, Carnot limits and the ideal gas law.
heat_equation
    The 1-D heat equation: Fourier-series, FTCS and Crank-Nicolson solvers and
    the :class:`HeatRod` convenience class.

Every public name is also re-exported here.
"""

from physicsforge.thermodynamics import heat, heat_equation
from physicsforge.thermodynamics.heat import *  # noqa: F403
from physicsforge.thermodynamics.heat_equation import *  # noqa: F403

__all__ = ["heat", "heat_equation", *heat.__all__, *heat_equation.__all__]

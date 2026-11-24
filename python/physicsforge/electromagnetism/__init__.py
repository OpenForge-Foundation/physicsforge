"""Electromagnetism.

Submodules
----------
circuits
    Ohm's law, power, voltage/current dividers, R/C networks, RC transients.
electrostatics
    Coulomb's law, fields and potentials of point charges, capacitors.
maxwell
    Maxwell's equations (Gauss, Faraday, Ampère-Maxwell), magnetic fields of
    standard geometries, the Lorentz force and EM waves.

Every public function is also re-exported here.
"""

from physicsforge.electromagnetism import circuits, electrostatics, maxwell
from physicsforge.electromagnetism.circuits import *  # noqa: F403
from physicsforge.electromagnetism.electrostatics import *  # noqa: F403
from physicsforge.electromagnetism.maxwell import *  # noqa: F403

__all__ = [
    "circuits",
    "electrostatics",
    "maxwell",
    *circuits.__all__,
    *electrostatics.__all__,
    *maxwell.__all__,
]

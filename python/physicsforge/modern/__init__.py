"""Modern physics.

Submodules
----------
quantum
    Planck-Einstein relation, photons, de Broglie waves, photoelectric effect,
    Planck's law, Wien's law, the Bohr atom and the uncertainty principle.
relativity
    Lorentz factor, time dilation, length contraction and mass-energy.

Every public function is also re-exported here.
"""

from physicsforge.modern import quantum, relativity
from physicsforge.modern.quantum import *  # noqa: F403
from physicsforge.modern.relativity import *  # noqa: F403

__all__ = ["quantum", "relativity", *quantum.__all__, *relativity.__all__]

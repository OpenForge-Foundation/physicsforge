"""Classical (Newtonian) mechanics.

Submodules
----------
kinematics
    Constant-acceleration motion and projectiles.
dynamics
    Newton's laws, momentum, work, energy, power, springs and friction.
gravitation
    Universal gravitation and orbital mechanics.
rotation
    Torque, angular momentum, moments of inertia and circular motion.
rocketry
    The Tsiolkovsky rocket equation, specific impulse and thrust.

Every public function is also re-exported here, e.g.
``physicsforge.mechanics.net_force``.
"""

from physicsforge.mechanics import dynamics, gravitation, kinematics, rocketry, rotation
from physicsforge.mechanics.dynamics import *  # noqa: F403
from physicsforge.mechanics.gravitation import *  # noqa: F403
from physicsforge.mechanics.kinematics import *  # noqa: F403
from physicsforge.mechanics.rocketry import *  # noqa: F403
from physicsforge.mechanics.rotation import *  # noqa: F403

__all__ = [
    "dynamics",
    "gravitation",
    "kinematics",
    "rocketry",
    "rotation",
    *dynamics.__all__,
    *gravitation.__all__,
    *kinematics.__all__,
    *rocketry.__all__,
    *rotation.__all__,
]

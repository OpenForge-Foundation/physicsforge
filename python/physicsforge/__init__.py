"""physicsforge: high-performance physics with a C++20 core.

physicsforge is the OpenForge Foundation's physics library, a sibling of
probforge. Every law is implemented once in C++ and exposed through a typed,
documented Python API that broadcasts over NumPy arrays.

Packages
--------
constants
    CODATA 2022 physical constants and mathematical constants (π, e, ...).
vector
    :class:`Vector3`, dot and cross products.
mechanics
    Kinematics, Newton's laws, gravitation, rotation/torque, rocketry.
electromagnetism
    Ohm's law and circuits, electrostatics, Maxwell's equations.
thermodynamics
    Heat laws (Fourier, Newton, Stefan-Boltzmann) and the heat equation.
modern
    Quantum (Planck) relations and special relativity.
numerics
    Newton's method, Hessians, vector calculus, quadrature, FFT, Fourier
    series and Laplace transforms.

Conventions
-----------
* SI units everywhere; angles in radians; absolute temperatures in kelvin.
* Scalar laws accept floats or array-likes, broadcast like NumPy ufuncs, and
  return ``float`` for scalar input.
* Physically invalid input raises :class:`ValueError` with a descriptive message.

Examples
--------
>>> import physicsforge as pf
>>> pf.electromagnetism.voltage_divider(12.0, r1=1_000.0, r2=2_000.0)
8.0
>>> pf.cross([1, 0, 0], [0, 1, 0])
Vector3(0.0, 0.0, 1.0)
"""

from __future__ import annotations

try:
    from physicsforge import _physicsforge
except ImportError as exc:  # pragma: no cover - only hit on a broken install
    raise ImportError(
        "physicsforge's compiled core (_physicsforge) could not be imported. "
        "Reinstall the wheel for your Python version or rebuild from source (see README)."
    ) from exc

from physicsforge import constants, electromagnetism, mechanics, modern, numerics, thermodynamics, vector
from physicsforge.exceptions import ConvergenceError, PhysicsForgeError
from physicsforge.vector import Vector3, cross, dot

__version__: str = _physicsforge.__version__

__all__ = [
    "ConvergenceError",
    "PhysicsForgeError",
    "Vector3",
    "__version__",
    "constants",
    "cross",
    "dot",
    "electromagnetism",
    "mechanics",
    "modern",
    "numerics",
    "thermodynamics",
    "vector",
]

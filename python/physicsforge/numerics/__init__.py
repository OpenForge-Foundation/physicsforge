"""Numerical methods used across physics.

Submodules
----------
calculus
    Derivatives, gradient, Hessian, Jacobian, Laplacian, divergence and curl.
roots
    Newton-Raphson (scalar and systems) and bisection.
optimize
    Newton's method for minimisation using the Hessian.
quadrature
    Trapezoid, Simpson and adaptive integration.
ode
    Euler, semi-implicit (symplectic) Euler and Runge-Kutta 4 integrators.
transforms
    FFT, Fourier series, and forward/inverse Laplace transforms.

Every public name is also re-exported here.
"""

from physicsforge.numerics import calculus, ode, optimize, quadrature, roots, transforms
from physicsforge.numerics.calculus import *  # noqa: F403
from physicsforge.numerics.ode import *  # noqa: F403
from physicsforge.numerics.optimize import *  # noqa: F403
from physicsforge.numerics.quadrature import *  # noqa: F403
from physicsforge.numerics.roots import *  # noqa: F403
from physicsforge.numerics.transforms import *  # noqa: F403

__all__ = [
    "calculus",
    "ode",
    "optimize",
    "quadrature",
    "roots",
    "transforms",
    *calculus.__all__,
    *ode.__all__,
    *optimize.__all__,
    *quadrature.__all__,
    *roots.__all__,
    *transforms.__all__,
]

r"""The one-dimensional heat equation and Fourier's series solution.

.. math:: \frac{\partial T}{\partial t} = \alpha \frac{\partial^2 T}{\partial x^2},
          \qquad T(0,t) = T_L,\quad T(L,t) = T_R

Fourier invented his series to solve exactly this problem. Subtracting the
linear steady state :math:`T_s(x)` leaves a transient with zero end values whose
sine-series modes decay independently:

.. math:: T(x,t) = T_s(x) + \sum_{n\ge1} b_n \sin\frac{n\pi x}{L}\, e^{-\alpha (n\pi/L)^2 t}

Three solvers are provided, all implemented in C++:

* :func:`heat_fourier_coefficients` / :func:`heat_fourier_solution`: the
  analytic series. Exact apart from truncation; evaluates any time directly.
* :func:`heat_equation_explicit`: FTCS finite differences. Cheap per step but
  stable only for Fourier number :math:`r = \alpha\,\Delta t/\Delta x^2 \le 1/2`.
* :func:`heat_equation_crank_nicolson`: implicit, second order, unconditionally
  stable. Each step is an :math:`O(N)` tridiagonal solve.

:class:`HeatRod` wraps all three behind one object.

Examples
--------
>>> import numpy as np
>>> from physicsforge.thermodynamics import HeatRod
>>> rod = HeatRod(length=1.0, diffusivity=1e-2, n_points=51)
>>> initial = np.sin(np.pi * rod.grid)                      # single Fourier mode
>>> exact = np.exp(-1e-2 * np.pi**2 * 5.0)                  # its decay after 5 s
>>> bool(abs(rod.fourier_solution(initial, time=5.0)[25] - exact) < 1e-12)
True
>>> bool(abs(rod.solve(initial, time=5.0)[25] - exact) < 1e-3)
True
"""

from __future__ import annotations

import math
from collections.abc import Callable
from typing import Literal

import numpy as np

from physicsforge import _physicsforge as _core
from physicsforge._typing import ArrayLike, FloatArray, FloatOrArray

#: An initial temperature profile: a constant, samples on the rod's grid, or a
#: callable mapping positions (array) to temperatures (array).
InitialProfile = float | ArrayLike | Callable[[FloatArray], ArrayLike]


def fourier_number(diffusivity: ArrayLike, time_step: ArrayLike, grid_spacing: ArrayLike) -> FloatOrArray:
    r"""Grid Fourier number :math:`r = \alpha\,\Delta t / \Delta x^2`.

    FTCS is stable only for :math:`r \le 1/2`.

    Examples
    --------
    >>> fourier_number(1e-4, 0.5, 0.01)
    0.5
    """
    return _core.fourier_number(diffusivity, time_step, grid_spacing)


def heat_fourier_coefficients(
    initial_temperatures: ArrayLike,
    n_terms: int,
    left_temperature: float = 0.0,
    right_temperature: float = 0.0,
) -> FloatArray:
    r"""Sine-series coefficients :math:`b_1 \ldots b_N` of the transient initial profile.

    .. math:: b_n = \frac{2}{L}\int_0^L \big(T(x,0) - T_s(x)\big)\sin\frac{n\pi x}{L}\,dx

    Parameters
    ----------
    initial_temperatures : array_like
        Initial temperatures on a uniform grid that includes both ends (at least 3
        points). The coefficients do not depend on the rod length.
    n_terms : int
        Number of modes :math:`N` to compute (at least 1).
    left_temperature, right_temperature : float, optional
        Fixed end temperatures :math:`T_L, T_R`. Default 0.

    Returns
    -------
    numpy.ndarray
        Coefficients, shape ``(n_terms,)``.

    Examples
    --------
    >>> import numpy as np
    >>> x = np.linspace(0.0, 1.0, 101)
    >>> np.round(heat_fourier_coefficients(np.sin(np.pi * x) + 0.5 * np.sin(3 * np.pi * x), 3), 12) + 0.0
    array([1. , 0. , 0.5])
    """
    return _core.heat_fourier_coefficients(initial_temperatures, int(n_terms), left_temperature, right_temperature)


def heat_fourier_solution(
    coefficients: ArrayLike,
    length: float,
    diffusivity: float,
    time: float,
    positions: ArrayLike,
    left_temperature: float = 0.0,
    right_temperature: float = 0.0,
) -> FloatArray:
    r"""Evaluate the Fourier-series solution :math:`T(x, t)` at arbitrary positions.

    Parameters
    ----------
    coefficients : array_like
        Mode amplitudes from :func:`heat_fourier_coefficients`.
    length : float
        Rod length :math:`L` [m].
    diffusivity : float
        Thermal diffusivity :math:`\alpha` [m^2/s].
    time : float
        Time :math:`t \ge 0` [s].
    positions : array_like
        Positions in :math:`[0, L]` [m].
    left_temperature, right_temperature : float, optional
        End temperatures; must match those used for the coefficients.

    Returns
    -------
    numpy.ndarray
        Temperatures at ``positions``.
    """
    return _core.heat_fourier_solution(
        coefficients, length, diffusivity, time, positions, left_temperature, right_temperature
    )


def heat_equation_explicit(
    initial_temperatures: ArrayLike,
    diffusivity: float,
    grid_spacing: float,
    time_step: float,
    steps: int,
    left_temperature: float = 0.0,
    right_temperature: float = 0.0,
) -> FloatArray:
    r"""Advance a temperature profile with the explicit FTCS scheme.

    .. math:: T_i^{k+1} = T_i^k + r\,(T_{i-1}^k - 2T_i^k + T_{i+1}^k)

    Raises
    ------
    ValueError
        If the Fourier number :math:`r` exceeds 1/2 (the scheme would be unstable).
    """
    return _core.heat_equation_explicit(
        initial_temperatures, diffusivity, grid_spacing, time_step, int(steps), left_temperature, right_temperature
    )


def heat_equation_crank_nicolson(
    initial_temperatures: ArrayLike,
    diffusivity: float,
    grid_spacing: float,
    time_step: float,
    steps: int,
    left_temperature: float = 0.0,
    right_temperature: float = 0.0,
) -> FloatArray:
    r"""Advance a temperature profile with the Crank-Nicolson scheme.

    Second order in space and time and unconditionally stable. The constant
    tridiagonal system is factorised once, so each step costs :math:`O(N)`.
    """
    return _core.heat_equation_crank_nicolson(
        initial_temperatures, diffusivity, grid_spacing, time_step, int(steps), left_temperature, right_temperature
    )


class HeatRod:
    r"""A one-dimensional rod with fixed end temperatures, discretised on a uniform grid.

    Parameters
    ----------
    length : float
        Rod length :math:`L` [m]. Must be positive.
    diffusivity : float
        Thermal diffusivity :math:`\alpha` [m^2/s]. Must be positive. See
        :func:`physicsforge.thermodynamics.thermal_diffusivity`.
    left_temperature, right_temperature : float, optional
        Fixed (Dirichlet) end temperatures. Default 0.
    n_points : int, optional
        Number of grid points including both ends (at least 3). Default 101.

    Examples
    --------
    >>> rod = HeatRod(length=0.5, diffusivity=1e-4, left_temperature=100.0, n_points=11)
    >>> rod
    HeatRod(length=0.5, diffusivity=0.0001, left_temperature=100.0, right_temperature=0.0, n_points=11)
    >>> float(rod.steady_state()[5])
    50.0
    """

    __slots__ = ("_diffusivity", "_left", "_length", "_n_points", "_right")

    def __init__(
        self,
        length: float,
        diffusivity: float,
        *,
        left_temperature: float = 0.0,
        right_temperature: float = 0.0,
        n_points: int = 101,
    ) -> None:
        if not (math.isfinite(length) and length > 0):
            raise ValueError(f"length must be finite and positive (got {length})")
        if not (math.isfinite(diffusivity) and diffusivity > 0):
            raise ValueError(f"diffusivity must be finite and positive (got {diffusivity})")
        if not (math.isfinite(left_temperature) and math.isfinite(right_temperature)):
            raise ValueError("end temperatures must be finite")
        if int(n_points) != n_points or n_points < 3:
            raise ValueError(f"n_points must be an integer >= 3 (got {n_points})")
        self._length = float(length)
        self._diffusivity = float(diffusivity)
        self._left = float(left_temperature)
        self._right = float(right_temperature)
        self._n_points = int(n_points)

    # -- Properties ---------------------------------------------------------------------------

    @property
    def length(self) -> float:
        """Rod length [m]."""
        return self._length

    @property
    def diffusivity(self) -> float:
        """Thermal diffusivity [m^2/s]."""
        return self._diffusivity

    @property
    def left_temperature(self) -> float:
        """Fixed temperature at x = 0."""
        return self._left

    @property
    def right_temperature(self) -> float:
        """Fixed temperature at x = L."""
        return self._right

    @property
    def n_points(self) -> int:
        """Number of grid points, ends included."""
        return self._n_points

    @property
    def grid_spacing(self) -> float:
        r"""Grid spacing :math:`\Delta x` [m]."""
        return self._length / (self._n_points - 1)

    @property
    def grid(self) -> FloatArray:
        """Grid positions, shape ``(n_points,)``, from 0 to ``length`` inclusive."""
        return np.linspace(0.0, self._length, self._n_points)

    def __repr__(self) -> str:
        return (
            f"HeatRod(length={self._length!r}, diffusivity={self._diffusivity!r}, "
            f"left_temperature={self._left!r}, right_temperature={self._right!r}, n_points={self._n_points})"
        )

    # -- Solvers ------------------------------------------------------------------------------

    def steady_state(self) -> FloatArray:
        """The linear steady-state profile between the end temperatures."""
        return np.linspace(self._left, self._right, self._n_points)

    def fourier_solution(self, initial: InitialProfile, time: float, n_terms: int = 100) -> FloatArray:
        """Temperature on the grid at ``time`` from the analytic Fourier series.

        Parameters
        ----------
        initial : float, array_like or callable
            Initial profile: a constant, ``n_points`` samples, or ``f(x) -> T``.
        time : float
            Time [s], non-negative.
        n_terms : int, optional
            Number of sine modes retained (default 100). Must be below ``n_points``
            for the modes to be resolved by the grid.
        """
        samples = self._sample(initial)
        coefficients = heat_fourier_coefficients(samples, n_terms, self._left, self._right)
        return heat_fourier_solution(
            coefficients, self._length, self._diffusivity, time, self.grid, self._left, self._right
        )

    def solve(
        self,
        initial: InitialProfile,
        time: float,
        *,
        method: Literal["crank-nicolson", "explicit"] = "crank-nicolson",
        time_step: float | None = None,
    ) -> FloatArray:
        r"""Temperature on the grid at ``time`` by finite differences.

        Parameters
        ----------
        initial : float, array_like or callable
            Initial profile: a constant, ``n_points`` samples, or ``f(x) -> T``.
        time : float
            Final time [s], non-negative.
        method : {"crank-nicolson", "explicit"}, optional
            Time integrator. Default ``"crank-nicolson"``.
        time_step : float, optional
            Maximum time step [s]. Defaults to the value giving Fourier number
            0.4 (stable for both methods). The step is shrunk slightly so that an
            integer number of steps lands exactly on ``time``.

        Raises
        ------
        ValueError
            If ``method`` is unknown, or ``"explicit"`` is used with an unstable step.
        """
        if not (math.isfinite(time) and time >= 0):
            raise ValueError(f"time must be finite and non-negative (got {time})")
        dx = self.grid_spacing
        max_step = 0.4 * dx * dx / self._diffusivity if time_step is None else float(time_step)
        if not (math.isfinite(max_step) and max_step > 0):
            raise ValueError(f"time_step must be finite and positive (got {time_step})")
        steps = max(1, math.ceil(time / max_step)) if time > 0 else 0
        dt = time / steps if steps else 0.0

        samples = self._sample(initial)
        if method == "crank-nicolson":
            solver = heat_equation_crank_nicolson
        elif method == "explicit":
            solver = heat_equation_explicit
        else:
            raise ValueError(f"method must be 'crank-nicolson' or 'explicit' (got {method!r})")
        return solver(samples, self._diffusivity, dx, dt, steps, self._left, self._right)

    def _sample(self, initial: InitialProfile) -> FloatArray:
        if callable(initial):
            values = np.asarray(initial(self.grid), dtype=np.float64)
        else:
            values = np.asarray(initial, dtype=np.float64)
        if values.ndim == 0:
            values = np.full(self._n_points, float(values))
        if values.shape != (self._n_points,):
            raise ValueError(f"initial profile must have shape ({self._n_points},) (got {values.shape})")
        return values


__all__ = [
    "HeatRod",
    "InitialProfile",
    "fourier_number",
    "heat_equation_crank_nicolson",
    "heat_equation_explicit",
    "heat_fourier_coefficients",
    "heat_fourier_solution",
]

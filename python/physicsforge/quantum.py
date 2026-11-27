"""Foundational quantum relations built on Planck's constant.

Photon energy and momentum, matter waves, the photoelectric effect, black-body
radiation (Planck's law and Wien's law), the Bohr hydrogen atom and the
Heisenberg bound. All functions broadcast over NumPy arrays.
"""

from __future__ import annotations

from physicsforge import _physicsforge as _core
from physicsforge._typing import ArrayLike, FloatOrArray


def photon_energy_from_frequency(frequency: ArrayLike) -> FloatOrArray:
    r"""Planck-Einstein relation, :math:`E = h f`.

    Parameters
    ----------
    frequency : float or array_like
        Frequency :math:`f` [Hz]. Must be non-negative.

    Returns
    -------
    float or numpy.ndarray
        Photon energy [J].

    Examples
    --------
    >>> photon_energy_from_frequency(1.0)   # Planck's constant
    6.62607015e-34
    """
    return _core.photon_energy_from_frequency(frequency)


def photon_energy_from_wavelength(wavelength: ArrayLike) -> FloatOrArray:
    r"""Photon energy from wavelength, :math:`E = hc/\lambda`.

    Examples
    --------
    >>> from physicsforge.constants import ELECTRON_VOLT
    >>> round(photon_energy_from_wavelength(550e-9) / ELECTRON_VOLT, 3)   # green light, eV
    2.254
    """
    return _core.photon_energy_from_wavelength(wavelength)


def photon_momentum(wavelength: ArrayLike) -> FloatOrArray:
    r"""Photon momentum :math:`p = h/\lambda` [kg m/s]."""
    return _core.photon_momentum(wavelength)


def wavelength_from_frequency(frequency: ArrayLike) -> FloatOrArray:
    r"""Vacuum wavelength of light, :math:`\lambda = c/f` [m].

    Examples
    --------
    >>> wavelength_from_frequency(299_792_458.0)
    1.0
    """
    return _core.wavelength_from_frequency(frequency)


def de_broglie_wavelength(mass: ArrayLike, velocity: ArrayLike) -> FloatOrArray:
    r"""de Broglie wavelength of a non-relativistic particle, :math:`\lambda = h / mv`.

    Examples
    --------
    >>> from physicsforge.constants import ELECTRON_MASS
    >>> round(de_broglie_wavelength(ELECTRON_MASS, 1e6) * 1e9, 4)   # nm
    0.7274
    """
    return _core.de_broglie_wavelength(mass, velocity)


def photoelectric_kinetic_energy(frequency: ArrayLike, work_function: ArrayLike) -> FloatOrArray:
    r"""Einstein's photoelectric equation, :math:`K_{max} = hf - \phi`.

    Returns 0 below the threshold frequency :math:`\phi/h` (no photoemission).

    Parameters
    ----------
    frequency : float or array_like
        Incident light frequency [Hz].
    work_function : float or array_like
        Work function of the surface :math:`\phi` [J].
    """
    return _core.photoelectric_kinetic_energy(frequency, work_function)


def planck_spectral_radiance(wavelength: ArrayLike, temperature: ArrayLike) -> FloatOrArray:
    r"""Planck's law: black-body spectral radiance per unit wavelength.

    .. math:: B_\lambda(T) = \frac{2hc^2}{\lambda^5}\,\frac{1}{e^{hc/\lambda k_B T} - 1}

    Parameters
    ----------
    wavelength : float or array_like
        Wavelength :math:`\lambda` [m]. Must be positive.
    temperature : float or array_like
        Absolute temperature :math:`T` [K].

    Returns
    -------
    float or numpy.ndarray
        Spectral radiance [W sr^-1 m^-3].

    Notes
    -----
    Uses ``expm1`` so the long-wavelength (Rayleigh-Jeans) limit is computed
    without cancellation, and the deep short-wavelength tail underflows cleanly to 0.

    Examples
    --------
    >>> import numpy as np
    >>> wl = np.linspace(100e-9, 3000e-9, 2901)
    >>> peak = wl[np.argmax(planck_spectral_radiance(wl, 5772.0))]
    >>> round(peak * 1e9)   # the Sun peaks near 502 nm (Wien's law)
    502
    """
    return _core.planck_spectral_radiance(wavelength, temperature)


def wien_peak_wavelength(temperature: ArrayLike) -> FloatOrArray:
    r"""Wien's displacement law, :math:`\lambda_{max} = b/T`.

    Examples
    --------
    >>> round(wien_peak_wavelength(5772.0) * 1e9, 1)
    502.0
    """
    return _core.wien_peak_wavelength(temperature)


def bohr_energy_level(principal_quantum_number: ArrayLike) -> FloatOrArray:
    r"""Energy of hydrogen level :math:`n` in the Bohr model (infinite nuclear mass).

    .. math:: E_n = -\frac{m_e e^4}{8\varepsilon_0^2 h^2}\,\frac{1}{n^2}

    Examples
    --------
    >>> from physicsforge.constants import ELECTRON_VOLT
    >>> round(bohr_energy_level(1) / ELECTRON_VOLT, 4)
    -13.6057
    """
    return _core.bohr_energy_level(principal_quantum_number)


def heisenberg_min_momentum_uncertainty(position_uncertainty: ArrayLike) -> FloatOrArray:
    r"""Heisenberg bound on momentum uncertainty, :math:`\Delta p \ge \hbar / 2\Delta x`."""
    return _core.heisenberg_min_momentum_uncertainty(position_uncertainty)


__all__ = [
    "bohr_energy_level",
    "de_broglie_wavelength",
    "heisenberg_min_momentum_uncertainty",
    "photoelectric_kinetic_energy",
    "photon_energy_from_frequency",
    "photon_energy_from_wavelength",
    "photon_momentum",
    "planck_spectral_radiance",
    "wavelength_from_frequency",
    "wien_peak_wavelength",
]

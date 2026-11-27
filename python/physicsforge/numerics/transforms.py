r"""Fourier and Laplace transforms.

* :func:`fft` / :func:`ifft`: discrete Fourier transform of **any** length in
  :math:`O(n\log n)` (radix-2 Cooley-Tukey, plus Bluestein's algorithm for
  non-power-of-two sizes). Same convention as :func:`numpy.fft.fft`.
* :func:`fourier_series` / :func:`fourier_series_evaluate`: real trigonometric
  series :math:`f(x) \approx a_0/2 + \sum_k a_k\cos(2\pi kx/P) + b_k\sin(2\pi kx/P)`.
* :func:`laplace_transform`: :math:`F(s) = \int_0^\infty f(t)\,e^{-st}\,dt` for real :math:`s > 0`.
* :func:`inverse_laplace_transform`: the Gaver-Stehfest algorithm, which needs
  :math:`F` only on the positive real axis.
"""

from __future__ import annotations

from collections.abc import Callable

import numpy as np
from numpy.typing import NDArray

from physicsforge import _physicsforge as _core
from physicsforge._physicsforge import FourierSeries
from physicsforge._typing import ArrayLike, FloatArray


def fft(x: ArrayLike) -> NDArray[np.complex128]:
    r"""Discrete Fourier transform, :math:`X_k = \sum_n x_n e^{-2\pi i kn/N}`.

    Parameters
    ----------
    x : array_like
        Real or complex 1-D input of any length.

    Returns
    -------
    numpy.ndarray
        Complex spectrum, same length as ``x``.

    Examples
    --------
    >>> import numpy as np
    >>> x = np.random.default_rng(0).normal(size=12)      # not a power of two
    >>> bool(np.allclose(fft(x), np.fft.fft(x)))
    True
    """
    return _core.fft(x)


def ifft(x: ArrayLike) -> NDArray[np.complex128]:
    r"""Inverse discrete Fourier transform (includes the :math:`1/N` factor).

    Examples
    --------
    >>> import numpy as np
    >>> x = np.arange(7.0)
    >>> bool(np.allclose(ifft(fft(x)).real, x))
    True
    """
    return _core.ifft(x)


def fft_frequencies(n: int, sample_spacing: float = 1.0) -> FloatArray:
    """Sample frequencies matching :func:`fft` output order (like :func:`numpy.fft.fftfreq`).

    Examples
    --------
    >>> fft_frequencies(4, 0.25)
    array([ 0.,  1., -2., -1.])
    """
    return _core.fft_frequencies(int(n), sample_spacing)


def fourier_series(samples: ArrayLike, n_terms: int, period: float = 2 * np.pi) -> FourierSeries:
    r"""Real Fourier-series coefficients from one period of uniformly spaced samples.

    Parameters
    ----------
    samples : array_like
        :math:`N` samples over one period, end point excluded
        (e.g. ``f(np.linspace(0, P, N, endpoint=False))``).
    n_terms : int
        Highest harmonic :math:`K` to return. Requires :math:`2K < N` (Nyquist).
    period : float, optional
        Period :math:`P`. Default :math:`2\pi`.

    Returns
    -------
    FourierSeries
        With arrays ``a`` and ``b`` of length :math:`K + 1` (``b[0] == 0``) and ``period``.

    Examples
    --------
    >>> import numpy as np
    >>> x = np.linspace(0, 2 * np.pi, 64, endpoint=False)
    >>> s = fourier_series(1 + 3 * np.cos(x) - 2 * np.sin(3 * x), n_terms=4)
    >>> np.round(s.a, 12) + 0.0, np.round(s.b, 12) + 0.0   # + 0.0 normalises -0.0
    (array([2., 3., 0., 0., 0.]), array([ 0.,  0.,  0., -2.,  0.]))
    """
    return _core.fourier_series_coefficients(samples, int(n_terms), period)


def fourier_series_evaluate(series: FourierSeries, x: ArrayLike) -> FloatArray:
    """Evaluate a truncated Fourier series at the points ``x``.

    Examples
    --------
    >>> import numpy as np
    >>> s = FourierSeries(a=[2.0, 3.0], b=[0.0, 0.0], period=2 * np.pi)
    >>> fourier_series_evaluate(s, [0.0, np.pi])
    array([ 4., -2.])
    """
    return _core.fourier_series_evaluate(series, np.atleast_1d(np.asarray(x, dtype=np.float64)))


def laplace_transform(f: Callable[[float], float], s: float, *, tolerance: float = 1e-10) -> float:
    r"""Numerical Laplace transform :math:`F(s) = \int_0^\infty f(t)e^{-st}\,dt` for real :math:`s > 0`.

    The half-line is mapped onto :math:`[0, 1)` and integrated adaptively, with
    time rescaled by :math:`s` so accuracy does not depend on the size of :math:`s`.

    Examples
    --------
    >>> import math
    >>> round(laplace_transform(lambda t: math.sin(2 * t), 1.0), 9)   # 2 / (s^2 + 4)
    0.4
    """
    return _core.laplace_transform(f, s, tolerance)


def inverse_laplace_transform(transform: Callable[[float], float], t: float, *, n_terms: int = 14) -> float:
    r"""Numerical inverse Laplace transform by the Gaver-Stehfest algorithm.

    .. math:: f(t) \approx \frac{\ln 2}{t}\sum_{k=1}^{N} V_k\, F\!\left(\frac{k\ln 2}{t}\right)

    Parameters
    ----------
    transform : callable
        :math:`F(s)` for real :math:`s > 0`.
    t : float
        Time :math:`t > 0`.
    n_terms : int, optional
        Even number of terms in [2, 30]; default 14. Larger values amplify
        round-off error in double precision.

    Notes
    -----
    Accurate to about 4-6 significant digits for smooth, non-oscillatory
    :math:`f(t)`. Not suitable for strongly oscillatory or discontinuous signals.

    Examples
    --------
    >>> import math
    >>> abs(inverse_laplace_transform(lambda s: 1 / (s + 1), 1.0) - math.exp(-1)) < 1e-5
    True
    """
    return _core.inverse_laplace_transform(transform, t, int(n_terms))


__all__ = [
    "FourierSeries",
    "fft",
    "fft_frequencies",
    "fourier_series",
    "fourier_series_evaluate",
    "ifft",
    "inverse_laplace_transform",
    "laplace_transform",
]

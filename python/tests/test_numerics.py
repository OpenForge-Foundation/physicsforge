import math

import numpy as np
import pytest

from physicsforge import ConvergenceError
from physicsforge import numerics as nm


class TestCalculus:
    def test_derivatives(self):
        assert nm.derivative(math.exp, 1.0) == pytest.approx(math.e, rel=1e-11)
        assert nm.derivative(math.exp, 1.0, step=1e-3) == pytest.approx(math.e, rel=1e-10)
        assert nm.second_derivative(math.cos, 0.3) == pytest.approx(-math.cos(0.3), rel=1e-7)

    def test_gradient_and_hessian_of_quadratic_form(self, rng):
        a = rng.normal(size=(4, 4))
        q = a @ a.T + 4 * np.eye(4)
        b = rng.normal(size=4)
        f = lambda x: 0.5 * x @ q @ x - b @ x  # noqa: E731
        x = rng.normal(size=4)
        np.testing.assert_allclose(nm.gradient(f, x), q @ x - b, rtol=1e-8, atol=1e-8)
        h = nm.hessian(f, x)
        np.testing.assert_allclose(h, q, rtol=1e-5, atol=1e-5)
        np.testing.assert_array_equal(h, h.T)

    def test_jacobian_non_square(self):
        j = nm.jacobian(lambda v: [v[0] * v[1] * v[2], v[0] ** 2], [1.0, 2.0, 3.0])
        assert j.shape == (2, 3)
        np.testing.assert_allclose(j, [[6.0, 3.0, 2.0], [2.0, 0.0, 0.0]], atol=1e-9)

    def test_laplacian_of_harmonic_function_is_zero(self):
        assert nm.laplacian(lambda v: v[0] ** 2 - v[1] ** 2, [0.4, -1.3]) == pytest.approx(0.0, abs=1e-6)

    def test_vector_calculus_identities(self):
        field = lambda r: np.array([r[0] * r[1], r[1] * r[2], r[2] * r[0]])  # noqa: E731
        assert nm.divergence(field, (1, 2, 3)) == pytest.approx(6.0, rel=1e-10)
        np.testing.assert_allclose(np.asarray(nm.curl(field, (1, 2, 3))), [-2.0, -3.0, -1.0], rtol=1e-10)
        # div(curl F) = 0 and curl(grad f) = 0
        curl_field = lambda r: np.asarray(nm.curl(field, r))  # noqa: E731
        assert abs(nm.divergence(curl_field, (0.5, -0.2, 0.8))) < 1e-6
        grad_field = lambda r: nm.gradient(lambda v: v[0] ** 2 * v[1] + np.sin(v[2]), r)  # noqa: E731
        assert np.allclose(np.asarray(nm.curl(grad_field, (0.5, -0.2, 0.8))), 0.0, atol=1e-6)

    def test_rejects_empty_point(self):
        with pytest.raises(ValueError):
            nm.gradient(lambda v: 0.0, [])


class TestRoots:
    def test_newton_analytic_and_numeric_derivative(self):
        f = lambda x: x**3 - 2 * x - 5  # noqa: E731
        exact = 2.0945514815423265
        r1 = nm.newton(f, 2.0, fprime=lambda x: 3 * x**2 - 2)
        r2 = nm.newton(f, 2.0)
        assert r1.root == pytest.approx(exact, rel=1e-15)
        assert r2.root == pytest.approx(exact, rel=1e-13)
        assert r1.converged and r1.iterations <= 6

    def test_newton_quadratic_convergence(self):
        # From a nearby start, Newton on sqrt(2) needs very few iterations.
        assert nm.newton(lambda x: x * x - 2, 1.5, fprime=lambda x: 2 * x).iterations <= 5

    def test_convergence_error_carries_result(self):
        with pytest.raises(ConvergenceError) as info:
            nm.newton(lambda x: x * x + 1, 0.0, fprime=lambda x: 2 * x)
        assert info.value.result.converged is False
        result = nm.newton(lambda x: x * x + 1, 0.0, fprime=lambda x: 2 * x, raise_on_failure=False)
        assert not result.converged

    def test_bisect(self):
        assert nm.bisect(math.cos, 0.0, 3.0).root == pytest.approx(math.pi / 2, rel=1e-12)
        with pytest.raises(ValueError, match="opposite signs"):
            nm.bisect(lambda x: x * x + 1, -1.0, 1.0)

    def test_newton_system(self):
        sol = nm.newton_system(lambda v: [v[0] ** 2 + v[1] ** 2 - 4, np.exp(v[0]) + v[1] - 1], [1.0, -1.0])
        x, y = sol.root
        assert x**2 + y**2 == pytest.approx(4.0)
        assert math.exp(x) + y == pytest.approx(1.0)
        assert sol.residual_norm < 1e-10

    def test_python_exceptions_propagate_from_callbacks(self):
        def broken(x):
            raise RuntimeError("boom")

        with pytest.raises(RuntimeError, match="boom"):
            nm.newton(broken, 1.0)


class TestOptimize:
    def test_rosenbrock(self):
        result = nm.newton_minimize(lambda v: (1 - v[0]) ** 2 + 100 * (v[1] - v[0] ** 2) ** 2, [-1.2, 1.0])
        np.testing.assert_allclose(result.x, [1.0, 1.0], atol=1e-6)
        assert result.value < 1e-12
        assert result.converged

    def test_quadratic_in_one_newton_step(self, rng):
        a = rng.normal(size=(3, 3))
        q = a @ a.T + 3 * np.eye(3)
        target = np.array([1.0, -2.0, 0.5])
        result = nm.newton_minimize(lambda x: (x - target) @ q @ (x - target), np.zeros(3))
        np.testing.assert_allclose(result.x, target, atol=1e-7)
        assert result.iterations <= 3

    def test_escapes_saddle_region(self):
        # Hessian is indefinite at the start; the fallback must still find a minimum.
        f = lambda v: v[0] ** 4 - 2 * v[0] ** 2 + v[1] ** 2  # noqa: E731
        result = nm.newton_minimize(f, [0.1, 0.5])
        assert abs(result.x[0]) == pytest.approx(1.0, abs=1e-6)
        assert result.x[1] == pytest.approx(0.0, abs=1e-6)


class TestQuadrature:
    def test_rules(self):
        x = np.linspace(0, math.pi, 201)
        assert nm.simpson(np.sin(x), x[1] - x[0]) == pytest.approx(2.0, rel=1e-9)
        assert nm.trapezoid(np.sin(x), x[1] - x[0]) == pytest.approx(2.0, rel=1e-4)

    @pytest.mark.parametrize("n", [2, 3, 4, 5, 6, 7, 10, 11])
    def test_simpson_exact_for_cubics(self, n):
        x = np.linspace(0.0, 2.0, n)
        expected = 2.0**4 / 4 + 2.0**2
        y = x**3 + 2 * x
        if n == 2:
            pytest.skip("two samples fall back to the trapezoid rule")
        assert nm.simpson(y, x[1] - x[0]) == pytest.approx(expected, rel=1e-12)

    def test_adaptive(self):
        assert nm.integrate(lambda x: math.exp(-x * x), -8, 8) == pytest.approx(math.sqrt(math.pi), rel=1e-10)
        assert nm.integrate(math.sqrt, 0.0, 1.0, tolerance=1e-12) == pytest.approx(2 / 3, rel=1e-9)


class TestTransforms:
    @pytest.mark.parametrize("n", [1, 2, 3, 5, 8, 12, 17, 64, 100, 1000, 1024])
    def test_fft_matches_numpy(self, n, rng):
        x = rng.normal(size=n) + 1j * rng.normal(size=n)
        np.testing.assert_allclose(nm.fft(x), np.fft.fft(x), rtol=1e-10, atol=1e-10 * n)
        np.testing.assert_allclose(nm.ifft(nm.fft(x)), x, atol=1e-12 * n)

    def test_fft_frequencies_match_numpy(self):
        for n in (1, 6, 7):
            np.testing.assert_allclose(nm.fft_frequencies(n, 0.1), np.fft.fftfreq(n, 0.1))

    def test_fourier_series_of_square_wave(self):
        n = 4096
        x = (np.arange(n) + 0.5) * 2 * np.pi / n
        series = nm.fourier_series(np.sign(np.sin(x)), n_terms=9)
        magnitudes = np.hypot(series.a, series.b)
        for k in (1, 3, 5, 7, 9):
            assert magnitudes[k] == pytest.approx(4 / (math.pi * k), rel=1e-4)
        np.testing.assert_allclose(magnitudes[[0, 2, 4, 6, 8]], 0.0, atol=1e-12)

    def test_fourier_series_round_trip(self):
        period = 3.0
        x = np.linspace(0, period, 40, endpoint=False)
        f = lambda t: 0.5 + np.cos(2 * np.pi * t / period) + 0.25 * np.sin(6 * np.pi * t / period)  # noqa: E731
        series = nm.fourier_series(f(x), n_terms=5, period=period)
        probe = np.linspace(-1, 4, 23)
        np.testing.assert_allclose(nm.fourier_series_evaluate(series, probe), f(probe), atol=1e-12)

    def test_fourier_series_nyquist_limit(self):
        with pytest.raises(ValueError, match="Nyquist"):
            nm.fourier_series(np.ones(10), n_terms=5)

    @pytest.mark.parametrize(
        ("f", "s", "expected"),
        [
            (lambda t: 1.0, 2.0, 0.5),
            (lambda t: t * t, 1.5, 2 / 1.5**3),
            (lambda t: math.exp(-3 * t), 1.0, 0.25),
            (lambda t: math.cos(2 * t), 3.0, 3 / 13),
        ],
    )
    def test_laplace_transform_table(self, f, s, expected):
        assert nm.laplace_transform(f, s) == pytest.approx(expected, rel=1e-8)

    @pytest.mark.parametrize(
        ("transform", "f"),
        [
            (lambda s: 1 / (s + 2), lambda t: math.exp(-2 * t)),
            (lambda s: 1 / s**2, lambda t: t),
            (lambda s: 1 / (s * (s + 1)), lambda t: 1 - math.exp(-t)),
        ],
    )
    def test_inverse_laplace(self, transform, f):
        # Gaver-Stehfest (N=14) in double precision is good to ~1e-3 relative on decaying exponentials.
        for t in (0.5, 1.0, 2.0):
            assert nm.inverse_laplace_transform(transform, t) == pytest.approx(f(t), rel=2e-3, abs=1e-6)

    def test_laplace_round_trip(self):
        f = lambda t: t * math.exp(-t)  # noqa: E731
        approx = nm.inverse_laplace_transform(lambda s: nm.laplace_transform(f, s, tolerance=1e-12), 1.3)
        assert approx == pytest.approx(f(1.3), rel=1e-3)

// tests/cpp/test_numerics.cpp
#include <algorithm>
#include <cmath>
#include <complex>
#include <stdexcept>
#include <vector>

#include "physicsforge/constants.h"
#include "physicsforge/numerics/calculus.h"
#include "physicsforge/numerics/integrate.h"
#include "physicsforge/numerics/ode.h"
#include "physicsforge/numerics/optimize.h"
#include "physicsforge/numerics/roots.h"
#include "physicsforge/numerics/transforms.h"
#include "test_framework.h"

namespace {
    using physicsforge::math::Vec3;
    namespace c = physicsforge::constants;

    // Naive O(n^2) DFT used as the reference for the FFT.
    std::vector<std::complex<double>> naive_dft(const std::vector<std::complex<double>> &x) {
        const std::size_t n = x.size();
        std::vector<std::complex<double>> out(n);
        for (std::size_t k = 0; k < n; ++k) {
            for (std::size_t j = 0; j < n; ++j) {
                const double angle = -2.0 * c::PI * static_cast<double>(k * j) / static_cast<double>(n);
                out[k] += x[j] * std::polar(1.0, angle);
            }
        }
        return out;
    }

    void check_fft_length(const std::size_t n) {
        std::vector<std::complex<double>> x(n);
        for (std::size_t i = 0; i < n; ++i) {
            x[i] = {std::sin(0.3 * static_cast<double>(i)) + 0.1 * static_cast<double>(i),
                    std::cos(1.7 * static_cast<double>(i))};
        }
        const auto fast = physicsforge::numerics::fft(x);
        const auto slow = naive_dft(x);
        const auto back = physicsforge::numerics::ifft(fast);
        double max_err = 0.0;
        double max_round_trip = 0.0;
        for (std::size_t i = 0; i < n; ++i) {
            max_err = std::max(max_err, std::abs(fast[i] - slow[i]));
            max_round_trip = std::max(max_round_trip, std::abs(back[i] - x[i]));
        }
        PF_CHECK(max_err < 1e-9 * static_cast<double>(n));
        PF_CHECK(max_round_trip < 1e-12 * static_cast<double>(n));
    }
} // namespace

void test_numerics() {
    using namespace physicsforge::numerics;

    // --- Differentiation ---
    const ScalarFunction sin_fn = [](const double x) { return std::sin(x); };
    PF_CHECK_CLOSE_TOL(derivative(sin_fn, 1.0), std::cos(1.0), 1e-11, 0.0);
    PF_CHECK_CLOSE_TOL(second_derivative(sin_fn, 1.0), -std::sin(1.0), 1e-7, 0.0);

    // f(x, y) = x^2 y + 3 y^2  ->  grad = (2xy, x^2 + 6y), H = [[2y, 2x], [2x, 6]]
    const MultivariateFunction f = [](const std::vector<double> &v) { return v[0] * v[0] * v[1] + 3.0 * v[1] * v[1]; };
    const std::vector<double> p{1.5, -2.0};
    const auto g = gradient(f, p);
    PF_CHECK_CLOSE_TOL(g[0], -6.0, 1e-10, 0.0);
    PF_CHECK_CLOSE_TOL(g[1], 2.25 - 12.0, 1e-10, 0.0);
    const auto h = hessian(f, p);
    PF_CHECK_CLOSE_TOL(h[0], -4.0, 1e-6, 0.0);
    PF_CHECK_CLOSE_TOL(h[1], 3.0, 1e-6, 0.0);
    PF_CHECK(h[1] == h[2]); // exactly symmetric
    PF_CHECK_CLOSE_TOL(h[3], 6.0, 1e-6, 0.0);
    PF_CHECK_CLOSE_TOL(laplacian(f, p), -4.0 + 6.0, 1e-6, 0.0);

    const VectorFunction F = [](const std::vector<double> &v) {
        return std::vector<double>{v[0] * v[1], std::sin(v[0]) + v[1] * v[1] * v[1]};
    };
    const auto J = jacobian(F, std::vector<double>{0.5, 2.0});
    PF_CHECK_CLOSE_TOL(J[0], 2.0, 1e-10, 0.0);
    PF_CHECK_CLOSE_TOL(J[1], 0.5, 1e-10, 0.0);
    PF_CHECK_CLOSE_TOL(J[2], std::cos(0.5), 1e-10, 0.0);
    PF_CHECK_CLOSE_TOL(J[3], 12.0, 1e-10, 0.0);

    // Vector calculus: F = (x y, y z, z x)  ->  div = y + z + x, curl = (-y, -z, -x).
    const VectorField field = [](const Vec3 &r) { return Vec3{r.x * r.y, r.y * r.z, r.z * r.x}; };
    const Vec3 pt{1.0, 2.0, 3.0};
    PF_CHECK_CLOSE_TOL(divergence(field, pt), 6.0, 1e-10, 0.0);
    const Vec3 cu = curl(field, pt);
    PF_CHECK_CLOSE_TOL(cu.x, -2.0, 1e-10, 0.0);
    PF_CHECK_CLOSE_TOL(cu.y, -3.0, 1e-10, 0.0);
    PF_CHECK_CLOSE_TOL(cu.z, -1.0, 1e-10, 0.0);
    // Gradient fields are curl-free.
    const VectorField conservative = [](const Vec3 &r) { return Vec3{2 * r.x, 2 * r.y, 2 * r.z}; };
    PF_CHECK_CLOSE_TOL(physicsforge::math::norm(curl(conservative, pt)), 0.0, 0.0, 1e-9);

    // --- Root finding ---
    const ScalarFunction cubic = [](const double x) { return x * x * x - 2.0 * x - 5.0; };
    const ScalarFunction dcubic = [](const double x) { return 3.0 * x * x - 2.0; };
    constexpr double cubic_root = 2.0945514815423265;
    const RootResult nr = newton_raphson(cubic, dcubic, 2.0);
    PF_CHECK(nr.converged);
    PF_CHECK_CLOSE(nr.root, cubic_root);
    PF_CHECK(nr.iterations < 10);
    const RootResult nr_numeric = newton_raphson(cubic, 2.0);
    PF_CHECK(nr_numeric.converged);
    PF_CHECK_CLOSE_TOL(nr_numeric.root, cubic_root, 1e-12, 0.0);
    const RootResult bis = bisection(cubic, 2.0, 3.0);
    PF_CHECK(bis.converged);
    PF_CHECK_CLOSE_TOL(bis.root, cubic_root, 1e-11, 0.0);
    PF_CHECK_THROWS(bisection(cubic, 3.0, 4.0), std::invalid_argument);
    // Newton fails gracefully at a stationary point.
    const RootResult stuck =
            newton_raphson([](double x) { return x * x + 1.0; }, [](double x) { return 2.0 * x; }, 0.0);
    PF_CHECK(!stuck.converged);

    // Multivariate Newton: intersection of x^2 + y^2 = 4 and y = x  ->  (sqrt 2, sqrt 2).
    const VectorFunction system = [](const std::vector<double> &v) {
        return std::vector<double>{v[0] * v[0] + v[1] * v[1] - 4.0, v[1] - v[0]};
    };
    const SystemRootResult sys = newton_system(system, std::vector<double>{1.0, 0.5});
    PF_CHECK(sys.converged);
    PF_CHECK_CLOSE_TOL(sys.root[0], std::sqrt(2.0), 1e-12, 0.0);
    PF_CHECK_CLOSE_TOL(sys.root[1], std::sqrt(2.0), 1e-12, 0.0);

    // --- Newton minimisation with the Hessian ---
    // Rosenbrock function: minimum at (1, 1).
    const MultivariateFunction rosenbrock = [](const std::vector<double> &v) {
        const double a = 1.0 - v[0];
        const double b = v[1] - v[0] * v[0];
        return a * a + 100.0 * b * b;
    };
    const OptimizeResult opt = newton_minimize(rosenbrock, std::vector<double>{-1.2, 1.0});
    PF_CHECK(opt.converged);
    PF_CHECK_CLOSE_TOL(opt.x[0], 1.0, 1e-6, 0.0);
    PF_CHECK_CLOSE_TOL(opt.x[1], 1.0, 1e-6, 0.0);
    // A quadratic is minimised in a single Newton step.
    const MultivariateFunction bowl = [](const std::vector<double> &v) {
        return (v[0] - 3.0) * (v[0] - 3.0) + 2.0 * (v[1] + 1.0) * (v[1] + 1.0) + v[0] * v[1];
    };
    const OptimizeResult q = newton_minimize(bowl, std::vector<double>{0.0, 0.0});
    PF_CHECK(q.converged);
    PF_CHECK(q.iterations <= 3);
    PF_CHECK_CLOSE_TOL(q.x[0], 28.0 / 7.0, 1e-6, 0.0); // solves 2(x-3) + y = 0, 4(y+1) + x = 0
    PF_CHECK_CLOSE_TOL(q.x[1], -2.0, 1e-6, 0.0);

    // --- Integration ---
    std::vector<double> samples(101);
    for (std::size_t i = 0; i < samples.size(); ++i) {
        samples[i] = std::sin(c::PI * static_cast<double>(i) / 100.0);
    }
    PF_CHECK_CLOSE_TOL(simpson(samples, c::PI / 100.0), 2.0, 1e-8, 0.0);
    PF_CHECK_CLOSE_TOL(trapezoid(samples, c::PI / 100.0), 2.0, 2e-4, 0.0);
    samples.pop_back(); // odd interval count exercises the 3/8 panel
    std::vector<double> cubic_samples(10);
    for (std::size_t i = 0; i < cubic_samples.size(); ++i) {
        const double x = static_cast<double>(i);
        cubic_samples[i] = x * x * x;
    }
    PF_CHECK_CLOSE(simpson(cubic_samples, 1.0), std::pow(9.0, 4) / 4.0); // exact for cubics
    PF_CHECK_CLOSE_TOL(integrate(sin_fn, 0.0, c::PI), 2.0, 1e-10, 0.0);
    PF_CHECK_CLOSE_TOL(integrate(sin_fn, c::PI, 0.0), -2.0, 1e-10, 0.0);
    PF_CHECK_CLOSE_TOL(integrate([](double x) { return std::exp(-x * x); }, -6.0, 6.0), std::sqrt(c::PI), 1e-10, 0.0);

    // --- FFT ---
    for (const std::size_t n: {1u, 2u, 8u, 64u, 3u, 7u, 12u, 100u, 127u}) {
        check_fft_length(n);
    }
    PF_CHECK(fft(std::vector<std::complex<double>>{}).empty());
    const auto freqs = fft_frequencies(5, 0.1);
    PF_CHECK_CLOSE(freqs[1], 2.0);
    PF_CHECK_CLOSE(freqs[3], -4.0);

    // --- Fourier series: square wave sign(sin x) has b_k = 4 / (pi k) for odd k ---
    constexpr int n_samples = 4096;
    std::vector<double> square(n_samples);
    for (int i = 0; i < n_samples; ++i) {
        const double x = c::TAU * (i + 0.5) / n_samples;
        square[static_cast<std::size_t>(i)] = std::sin(x) >= 0.0 ? 1.0 : -1.0;
    }
    // Shift by half a sample: samples are at cell centres, so use the matching evaluation grid.
    const FourierSeries series = fourier_series_coefficients(square, 7, c::TAU);
    PF_CHECK_CLOSE_TOL(std::hypot(series.a[1], series.b[1]), 4.0 / c::PI, 1e-5, 0.0);
    PF_CHECK_CLOSE_TOL(std::hypot(series.a[2], series.b[2]), 0.0, 0.0, 1e-12);
    PF_CHECK_CLOSE_TOL(std::hypot(series.a[3], series.b[3]), 4.0 / (3.0 * c::PI), 1e-5, 0.0);
    // Round trip on a band-limited signal is exact.
    std::vector<double> band(32);
    for (std::size_t i = 0; i < band.size(); ++i) {
        const double x = 2.0 * static_cast<double>(i) / 32.0;
        band[i] = 1.0 + 3.0 * std::cos(c::PI * x) - 2.0 * std::sin(3.0 * c::PI * x);
    }
    const FourierSeries bl = fourier_series_coefficients(band, 5, 2.0);
    PF_CHECK_CLOSE_TOL(bl.a[0], 2.0, 1e-12, 0.0);
    PF_CHECK_CLOSE_TOL(bl.a[1], 3.0, 1e-12, 0.0);
    PF_CHECK_CLOSE_TOL(bl.b[3], -2.0, 1e-12, 0.0);
    const std::vector<double> probe{0.37};
    PF_CHECK_CLOSE_TOL(fourier_series_evaluate(bl, probe)[0],
                       1.0 + 3.0 * std::cos(c::PI * 0.37) - 2.0 * std::sin(3.0 * c::PI * 0.37), 1e-12, 0.0);
    PF_CHECK_THROWS(fourier_series_coefficients(band, 16, 2.0), std::invalid_argument);

    // --- Laplace transforms ---
    // L{exp(-a t)} = 1 / (s + a);  L{t} = 1 / s^2;  L{sin(w t)} = w / (s^2 + w^2).
    PF_CHECK_CLOSE_TOL(laplace_transform([](double t) { return std::exp(-2.0 * t); }, 3.0), 0.2, 1e-9, 0.0);
    PF_CHECK_CLOSE_TOL(laplace_transform([](double t) { return t; }, 2.0), 0.25, 1e-9, 0.0);
    PF_CHECK_CLOSE_TOL(laplace_transform([](double t) { return std::sin(2.0 * t); }, 1.0), 0.4, 1e-8, 0.0);
    PF_CHECK_CLOSE_TOL(inverse_laplace_transform([](double s) { return 1.0 / (s + 1.0); }, 1.0), std::exp(-1.0), 1e-4,
                       0.0);
    PF_CHECK_CLOSE_TOL(inverse_laplace_transform([](double s) { return 1.0 / (s * s); }, 2.5), 2.5, 1e-4, 0.0);
    PF_CHECK_THROWS(inverse_laplace_transform([](double s) { return 1.0 / s; }, 1.0, 7), std::invalid_argument);
}

void test_ode() {
    using namespace physicsforge::numerics;
    using State = std::vector<double>;

    // y' = -y, y(0) = 1  =>  y(1) = e^-1. Global error ~ h^p, so halving h divides it by 2^p.
    const OdeFunction decay = [](double, const State &y) { return State{-y[0]}; };
    const auto final_error = [&](const auto &method, int steps) {
        return std::abs(method(decay, 0.0, State{1.0}, 1.0, steps).y.back() - std::exp(-1.0));
    };
    const auto euler_method = [](const OdeFunction &f, double t0, const State &y0, double t1, int n) {
        return euler(f, t0, y0, t1, n);
    };
    const auto rk4_method = [](const OdeFunction &f, double t0, const State &y0, double t1, int n) {
        return runge_kutta4(f, t0, y0, t1, n);
    };
    PF_CHECK_CLOSE_TOL(final_error(euler_method, 200) / final_error(euler_method, 400), 2.0, 0.02, 0.0);
    PF_CHECK_CLOSE_TOL(final_error(rk4_method, 20) / final_error(rk4_method, 40), 16.0, 0.05, 0.0);
    PF_CHECK(final_error(rk4_method, 100) < 1e-10);

    // Trajectory layout: steps + 1 rows, exact end time.
    const OdeSolution sol = runge_kutta4(decay, 0.0, State{1.0}, 2.0, 7);
    PF_CHECK(sol.t.size() == 8 && sol.y.size() == 8 && sol.dimension == 1);
    PF_CHECK(sol.t.front() == 0.0 && sol.t.back() == 2.0);

    // Harmonic oscillator x'' = -x from (1, 0): energy (x^2 + v^2) / 2 starts at 0.5.
    const OdeFunction oscillator = [](double, const State &y) { return State{y[1], -y[0]}; };
    const AccelerationFunction spring = [](double, const State &x) { return State{-x[0]}; };
    constexpr int steps = 10000;
    const double t_end = 100.0 * c::TAU; // 100 periods
    const auto energy = [](double x, double v) { return 0.5 * (x * x + v * v); };

    const OdeSolution explicit_run = euler(oscillator, 0.0, State{1.0, 0.0}, t_end, steps);
    const double explicit_energy = energy(explicit_run.y[2 * steps], explicit_run.y[2 * steps + 1]);
    PF_CHECK(explicit_energy > 10.0); // forward Euler pumps energy into the oscillator

    const SecondOrderSolution symplectic = semi_implicit_euler(spring, 0.0, State{1.0}, State{0.0}, t_end, steps);
    double max_drift = 0.0;
    for (int i = 0; i <= steps; ++i) {
        const auto row = static_cast<std::size_t>(i);
        max_drift = std::max(max_drift, std::abs(energy(symplectic.position[row], symplectic.velocity[row]) - 0.5));
    }
    PF_CHECK(max_drift < 0.05); // bounded: O(h), no secular growth

    const OdeSolution rk4_run = runge_kutta4(oscillator, 0.0, State{1.0, 0.0}, t_end, steps);
    PF_CHECK_CLOSE_TOL(rk4_run.y[2 * steps], std::cos(t_end), 0.0, 1e-4);

    // Invalid input.
    PF_CHECK_THROWS(euler(decay, 0.0, State{}, 1.0, 10), std::invalid_argument);
    PF_CHECK_THROWS(euler(decay, 1.0, State{1.0}, 0.0, 10), std::invalid_argument);
    PF_CHECK_THROWS(runge_kutta4(decay, 0.0, State{1.0}, 1.0, 0), std::invalid_argument);
    PF_CHECK_THROWS(euler(oscillator, 0.0, State{1.0}, 1.0, 10), std::invalid_argument); // wrong size
    PF_CHECK_THROWS(semi_implicit_euler(spring, 0.0, State{1.0}, State{0.0, 1.0}, 1.0, 10), std::invalid_argument);
}

// python/bindings/numerics.cpp - exposes physicsforge::numerics
//
// Routines that take user callables keep the GIL (they call back into Python);
// pure-array routines (FFT, Fourier-series evaluation) release it.
#include <complex>
#include <sstream>
#include <vector>

#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "binding_utils.h"
#include "physicsforge/numerics/calculus.h"
#include "physicsforge/numerics/integrate.h"
#include "physicsforge/numerics/optimize.h"
#include "physicsforge/numerics/roots.h"
#include "physicsforge/numerics/transforms.h"

namespace py = pybind11;

namespace {
    // Copies a std::vector into a fresh NumPy array (used for read-only struct members).
    py::array_t<double> copy_to_numpy(const std::vector<double> &values) {
        return py::array_t<double>(static_cast<py::ssize_t>(values.size()), values.data());
    }

    std::string bool_text(const bool b) { return b ? "True" : "False"; }
} // namespace

void bind_numerics(py::module_ &m) {
    namespace nm = physicsforge::numerics;
    using pfpy::as_span;
    using pfpy::ComplexArray;
    using pfpy::DoubleArray;
    using pfpy::to_numpy;

    // ---- Result types ------------------------------------------------------------------------
    py::class_<nm::RootResult>(m, "RootResult", "Outcome of a scalar root search.")
            .def_readonly("root", &nm::RootResult::root, "Final iterate.")
            .def_readonly("residual", &nm::RootResult::residual, "f(root).")
            .def_readonly("iterations", &nm::RootResult::iterations, "Iterations performed.")
            .def_readonly("converged", &nm::RootResult::converged, "Whether the tolerance was met.")
            .def("__repr__", [](const nm::RootResult &r) {
                std::ostringstream out;
                out.precision(17);
                out << "RootResult(root=" << r.root << ", residual=" << r.residual << ", iterations=" << r.iterations
                    << ", converged=" << bool_text(r.converged) << ")";
                return out.str();
            });

    py::class_<nm::SystemRootResult>(m, "SystemRootResult", "Outcome of a multivariate root search.")
            .def_property_readonly("root", [](const nm::SystemRootResult &r) { return copy_to_numpy(r.root); })
            .def_readonly("residual_norm", &nm::SystemRootResult::residual_norm)
            .def_readonly("iterations", &nm::SystemRootResult::iterations)
            .def_readonly("converged", &nm::SystemRootResult::converged)
            .def("__repr__", [](const nm::SystemRootResult &r) {
                std::ostringstream out;
                out << "SystemRootResult(residual_norm=" << r.residual_norm << ", iterations=" << r.iterations
                    << ", converged=" << bool_text(r.converged) << ")";
                return out.str();
            });

    py::class_<nm::OptimizeResult>(m, "OptimizeResult", "Outcome of a minimisation.")
            .def_property_readonly("x", [](const nm::OptimizeResult &r) { return copy_to_numpy(r.x); })
            .def_readonly("value", &nm::OptimizeResult::value)
            .def_property_readonly("gradient", [](const nm::OptimizeResult &r) { return copy_to_numpy(r.gradient); })
            .def_readonly("iterations", &nm::OptimizeResult::iterations)
            .def_readonly("converged", &nm::OptimizeResult::converged)
            .def("__repr__", [](const nm::OptimizeResult &r) {
                std::ostringstream out;
                out.precision(17);
                out << "OptimizeResult(value=" << r.value << ", iterations=" << r.iterations
                    << ", converged=" << bool_text(r.converged) << ")";
                return out.str();
            });

    py::class_<nm::FourierSeries>(m, "FourierSeries", "Truncated real Fourier series a0/2 + sum(a_k cos + b_k sin).")
            .def(py::init([](const DoubleArray &a, const DoubleArray &b, const double period) {
                     const auto sa = as_span(a, "a");
                     const auto sb = as_span(b, "b");
                     return nm::FourierSeries{{sa.begin(), sa.end()}, {sb.begin(), sb.end()}, period};
                 }),
                 py::arg("a"), py::arg("b"), py::arg("period"))
            .def_property_readonly("a", [](const nm::FourierSeries &s) { return copy_to_numpy(s.a); })
            .def_property_readonly("b", [](const nm::FourierSeries &s) { return copy_to_numpy(s.b); })
            .def_readonly("period", &nm::FourierSeries::period)
            .def("__repr__", [](const nm::FourierSeries &s) {
                std::ostringstream out;
                out << "FourierSeries(n_terms=" << (s.a.empty() ? 0 : s.a.size() - 1) << ", period=" << s.period << ")";
                return out.str();
            });

    // ---- Calculus ----------------------------------------------------------------------------
    m.def("derivative", &nm::derivative, py::arg("f"), py::arg("x"), py::arg("step") = 0.0);
    m.def("second_derivative", &nm::second_derivative, py::arg("f"), py::arg("x"), py::arg("step") = 0.0);
    m.def(
            "gradient",
            [](const nm::MultivariateFunction &f, const std::vector<double> &x) {
                return to_numpy(nm::gradient(f, x));
            },
            py::arg("f"), py::arg("x"));
    m.def(
            "hessian",
            [](const nm::MultivariateFunction &f, const std::vector<double> &x) {
                const auto n = static_cast<py::ssize_t>(x.size());
                return to_numpy(nm::hessian(f, x)).reshape({n, n});
            },
            py::arg("f"), py::arg("x"));
    m.def(
            "jacobian",
            [](const nm::VectorFunction &f, const std::vector<double> &x) {
                const auto n = static_cast<py::ssize_t>(x.size());
                py::array_t<double> flat = to_numpy(nm::jacobian(f, x));
                return flat.reshape({flat.size() / n, n});
            },
            py::arg("f"), py::arg("x"));
    m.def("laplacian", &nm::laplacian, py::arg("f"), py::arg("x"));
    m.def("divergence", &nm::divergence, py::arg("field"), py::arg("point"));
    m.def("curl", &nm::curl, py::arg("field"), py::arg("point"));

    // ---- Roots and optimisation --------------------------------------------------------------
    m.def("newton_raphson",
          py::overload_cast<const nm::ScalarFunction &, const nm::ScalarFunction &, double, double, int>(
                  &nm::newton_raphson),
          py::arg("f"), py::arg("derivative"), py::arg("x0"), py::arg("tolerance"), py::arg("max_iterations"));
    m.def("newton_raphson_numeric",
          py::overload_cast<const nm::ScalarFunction &, double, double, int>(&nm::newton_raphson), py::arg("f"),
          py::arg("x0"), py::arg("tolerance"), py::arg("max_iterations"));
    m.def("bisection", &nm::bisection, py::arg("f"), py::arg("a"), py::arg("b"), py::arg("tolerance"),
          py::arg("max_iterations"));
    m.def("newton_system", &nm::newton_system, py::arg("f"), py::arg("x0"), py::arg("tolerance"),
          py::arg("max_iterations"));
    m.def("newton_minimize", &nm::newton_minimize, py::arg("f"), py::arg("x0"), py::arg("tolerance"),
          py::arg("max_iterations"));

    // ---- Integration -------------------------------------------------------------------------
    m.def(
            "trapezoid", [](const DoubleArray &y, const double dx) { return nm::trapezoid(as_span(y, "samples"), dx); },
            py::arg("samples"), py::arg("spacing"));
    m.def(
            "simpson", [](const DoubleArray &y, const double dx) { return nm::simpson(as_span(y, "samples"), dx); },
            py::arg("samples"), py::arg("spacing"));
    m.def("integrate", &nm::integrate, py::arg("f"), py::arg("a"), py::arg("b"), py::arg("tolerance"),
          py::arg("max_depth"));

    // ---- Transforms --------------------------------------------------------------------------
    m.def(
            "fft",
            [](const ComplexArray &x) {
                const auto span = as_span(x, "x");
                std::vector<std::complex<double>> out;
                {
                    py::gil_scoped_release release;
                    out = nm::fft(span);
                }
                return to_numpy(std::move(out));
            },
            py::arg("x"));
    m.def(
            "ifft",
            [](const ComplexArray &x) {
                const auto span = as_span(x, "x");
                std::vector<std::complex<double>> out;
                {
                    py::gil_scoped_release release;
                    out = nm::ifft(span);
                }
                return to_numpy(std::move(out));
            },
            py::arg("x"));
    m.def(
            "fft_frequencies", [](const int n, const double d) { return to_numpy(nm::fft_frequencies(n, d)); },
            py::arg("n"), py::arg("sample_spacing"));
    m.def(
            "fourier_series_coefficients",
            [](const DoubleArray &samples, const int n_terms, const double period) {
                return nm::fourier_series_coefficients(as_span(samples, "samples"), n_terms, period);
            },
            py::arg("samples"), py::arg("n_terms"), py::arg("period"));
    m.def(
            "fourier_series_evaluate",
            [](const nm::FourierSeries &series, const DoubleArray &points) {
                const auto span = as_span(points, "points");
                std::vector<double> out;
                {
                    py::gil_scoped_release release;
                    out = nm::fourier_series_evaluate(series, span);
                }
                return to_numpy(std::move(out));
            },
            py::arg("series"), py::arg("points"));
    m.def("laplace_transform", &nm::laplace_transform, py::arg("f"), py::arg("s"), py::arg("tolerance"));
    m.def("inverse_laplace_transform", &nm::inverse_laplace_transform, py::arg("transform"), py::arg("t"),
          py::arg("n_terms"));
}

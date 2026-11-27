// src/bindings/vector.cpp - exposes physicsforge::math::Vec3 as Vector3

#include <pybind11/numpy.h>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "physicsforge/math/vector3.h"

namespace py = pybind11;
using physicsforge::math::Vec3;

namespace {
    double component(const Vec3 &v, py::ssize_t index) {
        if (index < 0) {
            index += 3;
        }
        switch (index) {
            case 0:
                return v.x;
            case 1:
                return v.y;
            case 2:
                return v.z;
            default:
                throw py::index_error("Vector3 index out of range");
        }
    }

    // Uses Python's float repr (shortest round-trip form), e.g. Vector3(1.0, 0.1, 3e-08).
    py::str repr(const Vec3 &v) { return py::str("Vector3({!r}, {!r}, {!r})").format(v.x, v.y, v.z); }
} // namespace

void bind_vector(py::module_ &m) {
    namespace pm = physicsforge::math;

    py::class_<Vec3>(m, "Vector3", R"pbdoc(
A vector in R^3 with Cartesian components (x, y, z).

Supports +, - (binary and unary), scalar * and /, ==, indexing, iteration,
len() == 3, pickling, and conversion with numpy.asarray().
)pbdoc")
            .def(py::init<>(), "Zero vector.")
            .def(py::init([](const double x, const double y, const double z) { return Vec3{x, y, z}; }), py::arg("x"),
                 py::arg("y"), py::arg("z"))
            .def_readwrite("x", &Vec3::x, "x component")
            .def_readwrite("y", &Vec3::y, "y component")
            .def_readwrite("z", &Vec3::z, "z component")
            .def(py::self + py::self)
            .def(py::self - py::self)
            .def(-py::self)
            .def(py::self * double())
            .def(double() * py::self)
            .def(py::self / double())
            .def(py::self == py::self)
            .def(py::self != py::self)
            .def("__getitem__", &component, py::arg("index"))
            .def("__len__", [](const Vec3 &) { return 3; })
            .def("__iter__", [](const Vec3 &v) { return py::iter(py::make_tuple(v.x, v.y, v.z)); })
            .def("__repr__", &repr)
            .def(
                    "__array__",
                    [](const Vec3 &v, const py::object &dtype, const py::object &) {
                        py::array_t<double> out(3);
                        auto r = out.mutable_unchecked<1>();
                        r(0) = v.x;
                        r(1) = v.y;
                        r(2) = v.z;
                        return dtype.is_none() ? py::object(out) : py::object(out.attr("astype")(dtype));
                    },
                    py::arg("dtype") = py::none(), py::arg("copy") = py::none())
            .def(
                    "to_tuple", [](const Vec3 &v) { return py::make_tuple(v.x, v.y, v.z); },
                    "Return the components as a tuple (x, y, z).")
            .def(
                    "dot", [](const Vec3 &a, const Vec3 &b) { return pm::dot(a, b); }, py::arg("other"),
                    "Dot product with another vector.")
            .def(
                    "cross", [](const Vec3 &a, const Vec3 &b) { return pm::cross(a, b); }, py::arg("other"),
                    "Cross product with another vector.")
            .def(
                    "norm", [](const Vec3 &a) { return pm::norm(a); }, "Euclidean length |v|.")
            .def(
                    "norm_squared", [](const Vec3 &a) { return pm::norm_squared(a); }, "Squared length |v|^2.")
            .def(
                    "normalized", [](const Vec3 &a) { return pm::normalized(a); },
                    "Unit vector in the same direction. Raises ValueError for the zero vector.")
            .def(py::pickle([](const Vec3 &v) { return py::make_tuple(v.x, v.y, v.z); },
                            [](const py::tuple &t) {
                                if (t.size() != 3) {
                                    throw std::runtime_error("invalid Vector3 pickle state");
                                }
                                return Vec3{t[0].cast<double>(), t[1].cast<double>(), t[2].cast<double>()};
                            }))
            .attr("__hash__") = py::none(); // mutable => unhashable

    m.def("vec_dot", [](const Vec3 &a, const Vec3 &b) { return pm::dot(a, b); }, py::arg("a"), py::arg("b"));
    m.def("vec_cross", [](const Vec3 &a, const Vec3 &b) { return pm::cross(a, b); }, py::arg("a"), py::arg("b"));
    m.def("vec_norm", [](const Vec3 &a) { return pm::norm(a); }, py::arg("a"));
    m.def("vec_normalized", &pm::normalized, py::arg("a"));
    m.def("vec_angle_between", &pm::angle_between, py::arg("a"), py::arg("b"));
    m.def("vec_project", &pm::project, py::arg("a"), py::arg("b"));
    m.def(
            "vec_triple_product",
            [](const Vec3 &a, const Vec3 &b, const Vec3 &c) { return pm::triple_product(a, b, c); }, py::arg("a"),
            py::arg("b"), py::arg("c"));
    m.def("vec_distance", &pm::distance, py::arg("a"), py::arg("b"));
}

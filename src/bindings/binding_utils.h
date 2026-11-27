// src/bindings/binding_utils.h
//
// Shared helpers for the pybind11 layer: zero-copy views of NumPy inputs and
// zero-copy hand-off of std::vector results to NumPy.
//
// Keeping these here (rather than in include/) keeps the C++ core free of any
// Python dependency, so it can be shipped and consumed as a plain C++ library.

#ifndef PHYSICSFORGE_BINDING_UTILS_H
#define PHYSICSFORGE_BINDING_UTILS_H

#include <complex>
#include <span>
#include <string>
#include <utility>
#include <vector>

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

namespace pfpy {

    namespace py = pybind11;

    // forcecast + c_style: any array-like of numbers is converted (copied only if needed)
    // to a contiguous float64 / complex128 array before reaching C++.
    using DoubleArray = py::array_t<double, py::array::c_style | py::array::forcecast>;
    using ComplexArray = py::array_t<std::complex<double>, py::array::c_style | py::array::forcecast>;

    // Read-only view of a 1-D NumPy array. Throws ValueError for other ranks.
    template<typename T, int Flags>
    std::span<const T> as_span(const py::array_t<T, Flags> &array, const char *name) {
        if (array.ndim() != 1) {
            throw py::value_error(std::string(name) + " must be 1-dimensional");
        }
        return {array.data(), static_cast<std::size_t>(array.shape(0))};
    }

    // Moves a std::vector into a NumPy array without copying; the array owns the buffer.
    template<typename T>
    py::array_t<T> to_numpy(std::vector<T> &&values) {
        auto *heap = new std::vector<T>(std::move(values));
        py::capsule owner(heap, [](void *p) { delete static_cast<std::vector<T> *>(p); });
        return py::array_t<T>(static_cast<py::ssize_t>(heap->size()), heap->data(), owner);
    }

} // namespace pfpy

#endif // PHYSICSFORGE_BINDING_UTILS_H

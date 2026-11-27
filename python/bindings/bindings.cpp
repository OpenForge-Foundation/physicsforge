// python/bindings/bindings.cpp
//
// Defines the compiled extension module physicsforge._physicsforge and wires in
// each domain's bindings. The public, documented Python API lives in
// python/physicsforge/; this module is an implementation detail.

#include <pybind11/pybind11.h>

#include "physicsforge/version.h"

namespace py = pybind11;

// Forward declarations of the per-domain binders (alphabetical, except vector which
// must run first because other domains use the Vector3 type in their signatures).
void bind_vector(py::module_ &m);
void bind_constants(py::module_ &m);
void bind_electromagnetism(py::module_ &m);
void bind_mechanics(py::module_ &m);
void bind_numerics(py::module_ &m);
void bind_quantum(py::module_ &m);
void bind_relativity(py::module_ &m);
void bind_thermodynamics(py::module_ &m);

PYBIND11_MODULE(_physicsforge, m) {
    m.doc() = "physicsforge compiled core (internal; use the physicsforge package API)";

    bind_vector(m);
    bind_constants(m);
    bind_electromagnetism(m);
    bind_mechanics(m);
    bind_numerics(m);
    bind_quantum(m);
    bind_relativity(m);
    bind_thermodynamics(m);

    m.attr("__version__") = PHYSICSFORGE_VERSION_STRING;
}

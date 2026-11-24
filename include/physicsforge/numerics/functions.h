// include/physicsforge/numerics/functions.h
//
// Callable types accepted by the numerical routines.
//
// std::function is used (rather than templates) so the routines can be compiled
// once into the library and bound to Python callables without code bloat. The
// type-erasure overhead is negligible relative to the user function itself.

#ifndef PHYSICSFORGE_NUMERICS_FUNCTIONS_H
#define PHYSICSFORGE_NUMERICS_FUNCTIONS_H

#include <functional>
#include <vector>

#include "physicsforge/math/vector3.h"

namespace physicsforge::numerics {

    // f : R -> R
    using ScalarFunction = std::function<double(double)>;
    // f : R^n -> R
    using MultivariateFunction = std::function<double(const std::vector<double> &)>;
    // F : R^n -> R^m
    using VectorFunction = std::function<std::vector<double>(const std::vector<double> &)>;
    // F : R^3 -> R^3 (a vector field, e.g. E(r) or B(r))
    using VectorField = std::function<math::Vec3(const math::Vec3 &)>;

} // namespace physicsforge::numerics

#endif // PHYSICSFORGE_NUMERICS_FUNCTIONS_H

// include/physicsforge/detail/validation.h
//
// Internal argument-validation helpers.
//
// Every public physicsforge function validates its physical preconditions (finite
// inputs, positive masses, non-zero distances, ...) and throws
// std::invalid_argument with a descriptive message when they are violated. The
// Python bindings translate this to ValueError; the C API catches it and returns
// NaN. The checks are a handful of branch-predicted comparisons on the hot path;
// message formatting only happens on the (cold) failure path.

#ifndef PHYSICSFORGE_DETAIL_VALIDATION_H
#define PHYSICSFORGE_DETAIL_VALIDATION_H

#include <cmath>
#include <sstream>
#include <stdexcept>
#include <string>

namespace physicsforge::detail {

    // Throws std::invalid_argument of the form "<name> must be <requirement> (got <value>)".
    [[noreturn]] inline void throw_invalid(const char *name, const char *requirement, double value) {
        std::ostringstream message;
        message << name << " must be " << requirement << " (got " << value << ")";
        throw std::invalid_argument(message.str());
    }

    // Throws std::invalid_argument with a free-form message.
    [[noreturn]] inline void throw_invalid(const std::string &message) { throw std::invalid_argument(message); }

    inline void require_finite(const double value, const char *name) {
        if (!std::isfinite(value)) [[unlikely]] {
            throw_invalid(name, "finite", value);
        }
    }

    inline void require_positive(const double value, const char *name) {
        if (!(std::isfinite(value) && value > 0.0)) [[unlikely]] {
            throw_invalid(name, "finite and positive", value);
        }
    }

    inline void require_non_negative(const double value, const char *name) {
        if (!(std::isfinite(value) && value >= 0.0)) [[unlikely]] {
            throw_invalid(name, "finite and non-negative", value);
        }
    }

    inline void require_non_zero(const double value, const char *name) {
        if (!(std::isfinite(value) && value != 0.0)) [[unlikely]] {
            throw_invalid(name, "finite and non-zero", value);
        }
    }

    // Requires lower <= value <= upper (closed interval).
    inline void require_in_range(const double value, const double lower, const double upper, const char *name) {
        if (!(std::isfinite(value) && value >= lower && value <= upper)) [[unlikely]] {
            std::ostringstream requirement;
            requirement << "in [" << lower << ", " << upper << "]";
            throw_invalid(name, requirement.str().c_str(), value);
        }
    }

} // namespace physicsforge::detail

#endif // PHYSICSFORGE_DETAIL_VALIDATION_H

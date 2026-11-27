/**
 * @file validation.h
 * @brief Internal argument validation and the C API's exception barrier.
 *
 * Every public function validates its physical preconditions (finite inputs, positive masses, non-zero
 * distances, ...) and throws @c std::invalid_argument with a message naming the parameter and its value.
 * The Python bindings translate this to @c ValueError; the C API catches it and returns NaN. The checks
 * are a few predictable comparisons on the hot path; messages are only formatted on failure.
 *
 * Not part of the public API: names live in physicsforge::detail and may change without notice.
 */

#ifndef PHYSICSFORGE_VALIDATION_H
#define PHYSICSFORGE_VALIDATION_H

#include <cmath>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>

namespace physicsforge::detail {

    /**
     * @brief Throws @c std::invalid_argument reading "<name> must be <requirement> (got <value>)".
     * @param name Parameter name.
     * @param requirement Human-readable requirement, e.g. "finite and positive".
     * @param value Offending value.
     */
    [[noreturn]] inline void throw_invalid(const char *name, const char *requirement, double value) {
        std::ostringstream message;
        message << name << " must be " << requirement << " (got " << value << ")";
        throw std::invalid_argument(message.str());
    }

    /**
     * @brief Throws @c std::invalid_argument with a free-form message.
     * @param message Error message.
     */
    [[noreturn]] inline void throw_invalid(const std::string &message) { throw std::invalid_argument(message); }

    /** @brief Requires @p value to be finite. */
    inline void require_finite(const double value, const char *name) {
        if (!std::isfinite(value)) [[unlikely]] {
            throw_invalid(name, "finite", value);
        }
    }

    /** @brief Requires @p value to be finite and strictly positive. */
    inline void require_positive(const double value, const char *name) {
        if (!(std::isfinite(value) && value > 0.0)) [[unlikely]] {
            throw_invalid(name, "finite and positive", value);
        }
    }

    /** @brief Requires @p value to be finite and non-negative. */
    inline void require_non_negative(const double value, const char *name) {
        if (!(std::isfinite(value) && value >= 0.0)) [[unlikely]] {
            throw_invalid(name, "finite and non-negative", value);
        }
    }

    /** @brief Requires @p value to be finite and non-zero. */
    inline void require_non_zero(const double value, const char *name) {
        if (!(std::isfinite(value) && value != 0.0)) [[unlikely]] {
            throw_invalid(name, "finite and non-zero", value);
        }
    }

    /** @brief Requires @p lower <= @p value <= @p upper (closed interval). */
    inline void require_in_range(const double value, const double lower, const double upper, const char *name) {
        if (!(std::isfinite(value) && value >= lower && value <= upper)) [[unlikely]] {
            std::ostringstream requirement;
            requirement << "in [" << lower << ", " << upper << "]";
            throw_invalid(name, requirement.str().c_str(), value);
        }
    }

    /**
     * @brief Exception barrier for the C ABI.
     *
     * C++ exceptions must never cross an @c extern "C" boundary, so every C entry point evaluates its body
     * through this function, which maps any exception to a quiet NaN.
     *
     * @param fn Callable returning @c double.
     * @return The result of @p fn, or NaN if it threw.
     */
    template<typename Fn>
    double guarded(Fn &&fn) noexcept {
        try {
            return std::forward<Fn>(fn)();
        } catch (...) {
            return std::numeric_limits<double>::quiet_NaN();
        }
    }

} // namespace physicsforge::detail

#endif // PHYSICSFORGE_VALIDATION_H

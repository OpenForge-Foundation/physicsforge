// src/api/api_guard.h
//
// Exception barrier for the C ABI. C++ exceptions must never cross an extern "C"
// boundary, so every C entry point evaluates its body through guarded(), which
// maps any exception to a quiet NaN.

#ifndef PHYSICSFORGE_API_GUARD_H
#define PHYSICSFORGE_API_GUARD_H

#include <limits>
#include <utility>

namespace physicsforge::api {

    template<typename Fn>
    double guarded(Fn &&fn) noexcept {
        try {
            return std::forward<Fn>(fn)();
        } catch (...) {
            return std::numeric_limits<double>::quiet_NaN();
        }
    }

} // namespace physicsforge::api

#endif // PHYSICSFORGE_API_GUARD_H

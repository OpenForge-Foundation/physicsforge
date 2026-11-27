// tests/cpp/test_framework.h
//
// Minimal, dependency-free test harness.
//
// Unlike <cassert>, these checks stay active in Release builds (NDEBUG), report
// the file/line and values of every failure, and keep running so a single run
// surfaces all failures. test_app.cpp returns a non-zero exit code if any check
// failed, which CTest and CI treat as a test failure.

#ifndef PHYSICSFORGE_TEST_FRAMEWORK_H
#define PHYSICSFORGE_TEST_FRAMEWORK_H

#include <algorithm>
#include <cmath>
#include <iostream>

namespace pftest {

    inline int &failure_count() {
        static int count = 0;
        return count;
    }

    inline int &check_count() {
        static int count = 0;
        return count;
    }

    inline void report(const char *file, const int line, const char *expression) {
        ++failure_count();
        std::cerr << "  FAILED " << file << ":" << line << ": " << expression << "\n";
    }

    // Relative-or-absolute closeness: |a - b| <= max(rel * max(|a|, |b|), abs).
    inline bool close(const double actual, const double expected, const double rel, const double abs_tol) {
        if (actual == expected) {
            return true;
        }
        const double scale = std::max(std::abs(actual), std::abs(expected));
        return std::abs(actual - expected) <= std::max(rel * scale, abs_tol);
    }

} // namespace pftest

#define PF_CHECK(condition)                                 \
    do {                                                    \
        ++pftest::check_count();                            \
        if (!(condition)) {                                 \
            pftest::report(__FILE__, __LINE__, #condition); \
        }                                                   \
    } while (false)

#define PF_CHECK_CLOSE_TOL(actual, expected, rel, abs_tol)                                    \
    do {                                                                                      \
        ++pftest::check_count();                                                              \
        const double pf_actual_ = (actual);                                                   \
        const double pf_expected_ = (expected);                                               \
        if (!pftest::close(pf_actual_, pf_expected_, (rel), (abs_tol))) {                     \
            pftest::report(__FILE__, __LINE__, #actual " ~= " #expected);                     \
            std::cerr.precision(17);                                                          \
            std::cerr << "    actual=" << pf_actual_ << " expected=" << pf_expected_ << "\n"; \
        }                                                                                     \
    } while (false)

#define PF_CHECK_CLOSE(actual, expected) PF_CHECK_CLOSE_TOL(actual, expected, 1e-12, 1e-300)

#define PF_CHECK_THROWS(expression, exception_type)                                     \
    do {                                                                                \
        ++pftest::check_count();                                                        \
        bool pf_threw_ = false;                                                         \
        try {                                                                           \
            static_cast<void>(expression);                                              \
        } catch (const exception_type &) {                                              \
            pf_threw_ = true;                                                           \
        } catch (...) {                                                                 \
        }                                                                               \
        if (!pf_threw_) {                                                               \
            pftest::report(__FILE__, __LINE__, #expression " throws " #exception_type); \
        }                                                                               \
    } while (false)

#endif // PHYSICSFORGE_TEST_FRAMEWORK_H

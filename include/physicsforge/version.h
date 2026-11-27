/**
 * @file version.h
 * @brief Library version and the symbol-visibility macro for the C API.
 *
 * The version follows Semantic Versioning (MAJOR.MINOR.PATCH). These macros are the single source of
 * truth: CMake, the C API and the Python package all read their version from here.
 */

#ifndef PHYSICSFORGE_VERSION_H
#define PHYSICSFORGE_VERSION_H

/**
 * @def PF_API
 * @brief Marks a C API symbol for export or import.
 *
 * When the C API is built as a shared library (@c PHYSICSFORGE_BUILD_C_SHARED=ON), CMake defines
 * @c PHYSICSFORGE_C_SHARED for the library and its consumers, and @c PHYSICSFORGE_C_EXPORTS only while
 * building the library itself. Otherwise the macro is empty.
 */
#if defined(PHYSICSFORGE_C_SHARED)
#if defined(_WIN32)
#if defined(PHYSICSFORGE_C_EXPORTS)
#define PF_API __declspec(dllexport)
#else
#define PF_API __declspec(dllimport)
#endif
#else
#define PF_API __attribute__((visibility("default")))
#endif
#else
#define PF_API
#endif

/** @brief Major version: incremented for incompatible API changes. */
#define PHYSICSFORGE_VERSION_MAJOR 0
/** @brief Minor version: incremented for backwards-compatible additions. */
#define PHYSICSFORGE_VERSION_MINOR 1
/** @brief Patch version: incremented for backwards-compatible fixes. */
#define PHYSICSFORGE_VERSION_PATCH 0

#define PHYSICSFORGE_STRINGIFY_IMPL(x) #x
#define PHYSICSFORGE_STRINGIFY(x) PHYSICSFORGE_STRINGIFY_IMPL(x)

/** @brief Full version string, e.g. "0.1.0". */
#define PHYSICSFORGE_VERSION_STRING                    \
    PHYSICSFORGE_STRINGIFY(PHYSICSFORGE_VERSION_MAJOR) \
    "." PHYSICSFORGE_STRINGIFY(PHYSICSFORGE_VERSION_MINOR) "." PHYSICSFORGE_STRINGIFY(PHYSICSFORGE_VERSION_PATCH)

#ifdef __cplusplus
extern "C" {
#endif

/** @brief Major version of the compiled library. */
PF_API int pf_version_major(void);
/** @brief Minor version of the compiled library. */
PF_API int pf_version_minor(void);
/** @brief Patch version of the compiled library. */
PF_API int pf_version_patch(void);
/** @brief Full version string of the compiled library (static storage, never NULL). */
PF_API const char *pf_version_string(void);

#ifdef __cplusplus
}
#endif

#endif // PHYSICSFORGE_VERSION_H

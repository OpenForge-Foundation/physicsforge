// include/physicsforge/version.h
//
// Version information for the physicsforge library.
//
// The version follows Semantic Versioning (https://semver.org): MAJOR.MINOR.PATCH.
// The macros below are the single source of truth; CMake, the C API and the Python
// package all derive their reported version from them.

#ifndef PHYSICSFORGE_VERSION_H
#define PHYSICSFORGE_VERSION_H

#include "physicsforge/export.h"

#define PHYSICSFORGE_VERSION_MAJOR 0
#define PHYSICSFORGE_VERSION_MINOR 1
#define PHYSICSFORGE_VERSION_PATCH 0

#define PHYSICSFORGE_STRINGIFY_IMPL(x) #x
#define PHYSICSFORGE_STRINGIFY(x) PHYSICSFORGE_STRINGIFY_IMPL(x)

// Full version string, e.g. "0.1.0".
#define PHYSICSFORGE_VERSION_STRING                    \
    PHYSICSFORGE_STRINGIFY(PHYSICSFORGE_VERSION_MAJOR) \
    "." PHYSICSFORGE_STRINGIFY(PHYSICSFORGE_VERSION_MINOR) "." PHYSICSFORGE_STRINGIFY(PHYSICSFORGE_VERSION_PATCH)

#ifdef __cplusplus
extern "C" {
#endif

// Returns the major, minor and patch version numbers of the compiled library.
PF_API int pf_version_major(void);
PF_API int pf_version_minor(void);
PF_API int pf_version_patch(void);
// Returns the full version string (static storage, never NULL).
PF_API const char *pf_version_string(void);

#ifdef __cplusplus
}
#endif

#endif // PHYSICSFORGE_VERSION_H

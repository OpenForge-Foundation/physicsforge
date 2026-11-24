// include/physicsforge/export.h
//
// Symbol visibility macro for the C API.
//
// When the C API is compiled into a shared library (PHYSICSFORGE_BUILD_C_SHARED=ON),
// CMake defines PHYSICSFORGE_C_SHARED for both the library and its consumers, and
// PHYSICSFORGE_C_EXPORTS only while building the library itself.

#ifndef PHYSICSFORGE_EXPORT_H
#define PHYSICSFORGE_EXPORT_H

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

#endif // PHYSICSFORGE_EXPORT_H

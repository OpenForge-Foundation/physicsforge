// src/api/version.cpp - physicsforge version C API
#include "physicsforge/version.h"

extern "C" int pf_version_major(void) { return PHYSICSFORGE_VERSION_MAJOR; }

extern "C" int pf_version_minor(void) { return PHYSICSFORGE_VERSION_MINOR; }

extern "C" int pf_version_patch(void) { return PHYSICSFORGE_VERSION_PATCH; }

extern "C" const char *pf_version_string(void) { return PHYSICSFORGE_VERSION_STRING; }

cmake_minimum_required(VERSION 3.20)

# -----------------------------------------------------------------------------
# Version: parsed from include/physicsforge/version.h (single source of truth)
# -----------------------------------------------------------------------------
file(STRINGS "${CMAKE_CURRENT_SOURCE_DIR}/include/physicsforge/version.h" _pf_version_lines
        REGEX "#define PHYSICSFORGE_VERSION_(MAJOR|MINOR|PATCH) ")
foreach (_line IN LISTS _pf_version_lines)
    string(REGEX MATCH "PHYSICSFORGE_VERSION_(MAJOR|MINOR|PATCH) ([0-9]+)" _ "${_line}")
    set(_PF_VERSION_${CMAKE_MATCH_1} "${CMAKE_MATCH_2}")
endforeach ()

project(physicsforge
        VERSION ${_PF_VERSION_MAJOR}.${_PF_VERSION_MINOR}.${_PF_VERSION_PATCH}
        DESCRIPTION "High-performance physics library with a C++20 core and Python bindings"
        LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if (NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif ()

# -----------------------------------------------------------------------------
# Options
# -----------------------------------------------------------------------------
option(PHYSICSFORGE_BUILD_TESTS "Build the C++ unit tests" ON)
option(PHYSICSFORGE_BUILD_PYTHON "Build the Python extension module" ON)
option(PHYSICSFORGE_BUILD_C_SHARED "Build the C API as a shared library (physicsforge_c)" OFF)
option(PHYSICSFORGE_WARNINGS_AS_ERRORS "Treat compiler warnings as errors" OFF)

# -----------------------------------------------------------------------------
# Compiler warnings (applied to every target we own)
# -----------------------------------------------------------------------------
add_library(physicsforge_warnings INTERFACE)
if (MSVC)
    target_compile_options(physicsforge_warnings INTERFACE /W4 /permissive- /utf-8)
    if (PHYSICSFORGE_WARNINGS_AS_ERRORS)
        target_compile_options(physicsforge_warnings INTERFACE /WX)
    endif ()
else ()
    target_compile_options(physicsforge_warnings INTERFACE -Wall -Wextra -Wpedantic -Wshadow -Wconversion
            -Wno-sign-conversion)
    if (PHYSICSFORGE_WARNINGS_AS_ERRORS)
        target_compile_options(physicsforge_warnings INTERFACE -Werror)
    endif ()
endif ()

# -----------------------------------------------------------------------------
# Core library (pure C++, no Python dependency)
# -----------------------------------------------------------------------------
set(PHYSICSFORGE_CORE_SOURCES
        src/math/linalg.cpp

        src/mechanics/kinematics.cpp
        src/mechanics/dynamics.cpp
        src/mechanics/gravitation.cpp
        src/mechanics/rotation.cpp
        src/mechanics/rocketry.cpp

        src/electromagnetism/circuits.cpp
        src/electromagnetism/electrostatics.cpp
        src/electromagnetism/maxwell.cpp

        src/thermodynamics/heat.cpp
        src/thermodynamics/heat_equation.cpp

        src/modern/quantum.cpp
        src/modern/relativity.cpp

        src/numerics/calculus.cpp
        src/numerics/integrate.cpp
        src/numerics/optimize.cpp
        src/numerics/roots.cpp
        src/numerics/transforms.cpp
)

set(PHYSICSFORGE_API_SOURCES
        src/api/version.cpp
        src/api/mechanics.cpp
        src/api/electromagnetism.cpp
        src/api/thermodynamics.cpp
        src/api/modern.cpp
)

add_library(physicsforge_core STATIC ${PHYSICSFORGE_CORE_SOURCES} ${PHYSICSFORGE_API_SOURCES})
add_library(physicsforge::core ALIAS physicsforge_core)
target_include_directories(physicsforge_core
        PUBLIC
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
)
target_compile_features(physicsforge_core PUBLIC cxx_std_20)
set_target_properties(physicsforge_core PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_link_libraries(physicsforge_core PRIVATE physicsforge_warnings)

# -----------------------------------------------------------------------------
# Optional shared C library exposing include/physicsforge/api.h
# -----------------------------------------------------------------------------
if (PHYSICSFORGE_BUILD_C_SHARED)
    add_library(physicsforge_c SHARED ${PHYSICSFORGE_CORE_SOURCES} ${PHYSICSFORGE_API_SOURCES})
    add_library(physicsforge::c ALIAS physicsforge_c)
    target_include_directories(physicsforge_c PUBLIC
            $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
            $<INSTALL_INTERFACE:include>)
    target_compile_definitions(physicsforge_c PUBLIC PHYSICSFORGE_C_SHARED PRIVATE PHYSICSFORGE_C_EXPORTS)
    set_target_properties(physicsforge_c PROPERTIES
            VERSION ${PROJECT_VERSION}
            SOVERSION ${PROJECT_VERSION_MAJOR}
            CXX_VISIBILITY_PRESET hidden
            VISIBILITY_INLINES_HIDDEN ON)
    target_link_libraries(physicsforge_c PRIVATE physicsforge_warnings)
endif ()

# -----------------------------------------------------------------------------
# Python extension module: physicsforge._physicsforge
# -----------------------------------------------------------------------------
if (PHYSICSFORGE_BUILD_PYTHON)
    find_package(Python 3.10 REQUIRED COMPONENTS Interpreter Development.Module)

    # Prefer the pybind11 installed alongside the build interpreter (pip install pybind11);
    # fall back to fetching a pinned release.
    execute_process(
            COMMAND "${Python_EXECUTABLE}" -m pybind11 --cmakedir
            OUTPUT_VARIABLE _pybind11_cmake_dir
            OUTPUT_STRIP_TRAILING_WHITESPACE
            ERROR_QUIET)
    if (_pybind11_cmake_dir)
        list(APPEND CMAKE_PREFIX_PATH "${_pybind11_cmake_dir}")
    endif ()
    find_package(pybind11 2.13 CONFIG QUIET)
    if (NOT pybind11_FOUND)
        message(STATUS "pybind11 not found locally; fetching v3.0.1")
        include(FetchContent)
        FetchContent_Declare(pybind11
                GIT_REPOSITORY https://github.com/pybind/pybind11.git
                GIT_TAG v3.0.1
                GIT_SHALLOW TRUE)
        FetchContent_MakeAvailable(pybind11)
    endif ()

    pybind11_add_module(_physicsforge
            python/bindings/bindings.cpp
            python/bindings/constants.cpp
            python/bindings/vector.cpp
            python/bindings/mechanics.cpp
            python/bindings/electromagnetism.cpp
            python/bindings/thermodynamics.cpp
            python/bindings/modern.cpp
            python/bindings/numerics.cpp
    )
    target_include_directories(_physicsforge PRIVATE python/bindings)
    target_link_libraries(_physicsforge PRIVATE physicsforge_core physicsforge_warnings)
    target_compile_definitions(_physicsforge PRIVATE PHYSICSFORGE_PYTHON_VERSION_INFO="${PROJECT_VERSION}")
endif ()

# -----------------------------------------------------------------------------
# Tests
# -----------------------------------------------------------------------------
if (PHYSICSFORGE_BUILD_TESTS)
    enable_testing()
    add_executable(physicsforge_tests
            tests/test_main.cpp
            tests/test_constants.cpp
            tests/test_vector.cpp
            tests/test_mechanics.cpp
            tests/test_electromagnetism.cpp
            tests/test_thermodynamics.cpp
            tests/test_modern.cpp
            tests/test_numerics.cpp
            tests/test_api.cpp
    )
    target_link_libraries(physicsforge_tests PRIVATE physicsforge_core physicsforge_warnings)
    add_test(NAME physicsforge_tests COMMAND physicsforge_tests)
endif ()

# -----------------------------------------------------------------------------
# Install (C++ headers + static core; optional shared C library)
# -----------------------------------------------------------------------------
include(GNUInstallDirs)
install(DIRECTORY include/physicsforge DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(TARGETS physicsforge_core ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
if (PHYSICSFORGE_BUILD_C_SHARED)
    install(TARGETS physicsforge_c
            LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
            ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
            RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
endif ()

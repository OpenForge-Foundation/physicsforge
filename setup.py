"""Build script for physicsforge.

The C++ core and the pybind11 extension are built by CMake; setuptools only
drives the build and packages the pure-Python layer around ``_physicsforge``.

Environment variables
---------------------
CMAKE_GENERATOR
    Override the CMake generator (e.g. ``Ninja``).
CMAKE_GENERATOR_INSTANCE
    Pin a specific Visual Studio installation on Windows (useful when several
    are installed and only one has the C++ workload).
CMAKE_ARGS
    Extra space-separated arguments forwarded to the CMake configure step.
CMAKE_BUILD_PARALLEL_LEVEL
    Parallel build jobs.
DEBUG
    ``1`` builds a Debug extension instead of Release.
"""

from __future__ import annotations

import os
import re
import subprocess
import sys
from pathlib import Path

from setuptools import Extension, find_packages, setup
from setuptools.command.build_ext import build_ext

ROOT = Path(__file__).parent.resolve()

# Convert distutils Windows platform specifiers to CMake -A arguments.
PLAT_TO_CMAKE = {
    "win32": "Win32",
    "win-amd64": "x64",
    "win-arm32": "ARM",
    "win-arm64": "ARM64",
}


def read_version() -> str:
    """Single source of truth: the version macros in ``include/physicsforge/version.h``."""
    header = (ROOT / "include" / "physicsforge" / "version.h").read_text(encoding="utf-8")
    parts = [
        re.search(rf"#define PHYSICSFORGE_VERSION_{part}\s+(\d+)", header).group(1)
        for part in ("MAJOR", "MINOR", "PATCH")
    ]
    return ".".join(parts)


class CMakeExtension(Extension):
    """An extension whose sources are built by CMake rather than setuptools."""

    def __init__(self, name: str, sourcedir: str = "") -> None:
        super().__init__(name, sources=[])
        self.sourcedir = os.fspath(Path(sourcedir).resolve())


class CMakeBuild(build_ext):
    """Configure and build ``CMakeExtension`` targets with CMake."""

    def build_extension(self, ext: CMakeExtension) -> None:
        ext_fullpath = Path.cwd() / self.get_ext_fullpath(ext.name)
        extdir = ext_fullpath.parent.resolve()

        debug = int(os.environ.get("DEBUG", 0)) if self.debug is None else self.debug
        cfg = "Debug" if debug else "Release"
        cmake_generator = os.environ.get("CMAKE_GENERATOR", "")

        cmake_args = [
            f"-DCMAKE_LIBRARY_OUTPUT_DIRECTORY={extdir}{os.sep}",
            f"-DPython_EXECUTABLE={sys.executable}",
            f"-DCMAKE_BUILD_TYPE={cfg}",
            "-DPHYSICSFORGE_BUILD_PYTHON=ON",
            "-DPHYSICSFORGE_BUILD_TESTS=OFF",
        ]
        build_args = ["--target", "_physicsforge"]

        if "CMAKE_ARGS" in os.environ:
            cmake_args += [item for item in os.environ["CMAKE_ARGS"].split(" ") if item]

        if self.compiler.compiler_type != "msvc":
            # Prefer Ninja when available: it parallelises automatically.
            if not cmake_generator or cmake_generator == "Ninja":
                try:
                    import ninja

                    ninja_executable_path = Path(ninja.BIN_DIR) / "ninja"
                    cmake_args += ["-GNinja", f"-DCMAKE_MAKE_PROGRAM:FILEPATH={ninja_executable_path}"]
                except ImportError:
                    pass
        else:
            single_config = any(x in cmake_generator for x in {"NMake", "Ninja"})
            contains_arch = any(x in cmake_generator for x in {"ARM", "Win64"})
            if not single_config and not contains_arch:
                cmake_args += ["-A", PLAT_TO_CMAKE[self.plat_name]]
            if not single_config:
                cmake_args += [f"-DCMAKE_LIBRARY_OUTPUT_DIRECTORY_{cfg.upper()}={extdir}"]
                build_args += ["--config", cfg]

        if sys.platform.startswith("darwin"):
            archs = re.findall(r"-arch (\S+)", os.environ.get("ARCHFLAGS", ""))
            if archs:
                cmake_args += ["-DCMAKE_OSX_ARCHITECTURES={}".format(";".join(archs))]

        if "CMAKE_BUILD_PARALLEL_LEVEL" not in os.environ and getattr(self, "parallel", None):
            build_args += [f"-j{self.parallel}"]

        build_temp = Path(self.build_temp) / ext.name
        build_temp.mkdir(parents=True, exist_ok=True)

        subprocess.run(["cmake", ext.sourcedir, *cmake_args], cwd=build_temp, check=True)
        subprocess.run(["cmake", "--build", ".", *build_args], cwd=build_temp, check=True)


setup(
    version=read_version(),
    package_dir={"": "python"},
    packages=find_packages(where="python", exclude=["tests", "tests.*"]),
    package_data={"physicsforge": ["py.typed"]},
    ext_modules=[CMakeExtension("physicsforge._physicsforge", sourcedir=str(ROOT))],
    cmdclass={"build_ext": CMakeBuild},
    zip_safe=False,
)

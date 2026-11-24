# physicsforge

Fast, validated physics for C++ and Python.

physicsforge is a C++20 library of physical laws and the numerical methods used with them. It has typed,
NumPy-vectorised Python bindings and a stable C ABI. It is the physics sibling of
[probforge](https://github.com/OpenForge-Foundation/probforge) and follows the same layout and conventions.

```python
>>> import physicsforge as pf
>>> pf.mechanics.rocket_delta_v(exhaust_velocity=4400.0, initial_mass=549_000.0, final_mass=25_600.0)
13488.2...
>>> pf.electromagnetism.voltage_divider(12.0, [1e3, 2e3], 2e3)        # broadcasts like NumPy
array([8., 6.])
>>> pf.cross((1, 0, 0), (0, 1, 0))
Vector3(0.0, 0.0, 1.0)
>>> rod = pf.thermodynamics.HeatRod(length=1.0, diffusivity=1e-4, left_temperature=100.0)
>>> rod.solve(initial=20.0, time=600.0)                                # Crank–Nicolson
array([100.        ,  98.1432395 , ...])
```

---

## Contents

- [Design principles](#design-principles)
- [Domain coverage](#domain-coverage)
- [Architecture](#architecture)
- [Building and installing](#building-and-installing)
- [Testing](#testing)
- [Adding a new law](#adding-a-new-law)
- [Conventions](#conventions)
- [Code formatting and pre-commit hooks](#code-formatting-and-pre-commit-hooks)
- [Release notes](#release-notes)
- [Roadmap](#roadmap)
- [Contributors](#contributors)

---

## Design principles

- **One implementation.** Every law is written once in C++. Python, the C ABI and native callers all call that same
  code.
- **Validate at the boundary.** Inputs that are physically meaningless throw `std::invalid_argument` with a message
  that names the parameter and shows the value, for example a negative mass, a NaN, or a final mass larger than the
  initial mass. Python receives this as `ValueError`; the C ABI returns `NaN`.
- **Numerically careful.** The formulas avoid cancellation and overflow:
  - `expm1` for propellant mass and Planck's law.
  - `hypot` for norms.
  - `atan2` for angles between vectors.
  - A cancellation-free relativistic kinetic energy.
  - A pivoted linear solve.
- **Fast by default.** Scalar laws are `py::vectorize`d, so arrays broadcast in C++ with no Python loop. The heat
  solvers release the GIL and hand results back to NumPy without copying.
- **SI units throughout.** Every docstring states each parameter's unit. CODATA 2022 constants.
- **Documented like a reference.** Every public Python function has a NumPy-style docstring with the formula,
  parameters, units, exceptions, and a doctest-verified example.

---

## Domain coverage

| Package                          | Contents                                                                                                                                                                                                                                                                                                                       |
|----------------------------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `physicsforge.constants`         | CODATA 2022: `c`, `h` (Planck), `hbar`, `e` (elementary charge), `k_B`, `N_A`, `GAS_CONSTANT`, `G`, `g0`, `epsilon_0`, `mu_0`, σ, Wien *b*, particle masses; plus mathematical π, τ and Euler's number                                                                                                                                       |
| `physicsforge.vector`            | `Vector3` (operators, NumPy interop, pickling), `dot`, `cross`, `norm`, `normalize`, `angle_between`, `project`, `triple_product`, `distance`                                                                                                                                                                                  |
| `physicsforge.mechanics`         | **Kinematics** (displacement, velocity, SUVAT, projectiles) · **Dynamics** (Newton's 2nd law, momentum, impulse, work, energy, power, Hooke, friction) · **Gravitation** (universal gravitation scalar and vector, field strength, escape/orbital velocity, Kepler's 3rd law) · **Rotation** (torque scalar and vector, angular momentum, moments of inertia, parallel axis) · **Rocketry** (Tsiolkovsky rocket equation, mass ratio, propellant mass, Isp, thrust, burn time) |
| `physicsforge.electromagnetism`  | **Circuits** (Ohm's law, power, voltage and current dividers, series/parallel networks, resistivity, RC transients) · **Electrostatics** (Coulomb's law, fields, potentials, capacitors) · **Maxwell** (Gauss's laws, Faraday–Lenz, Ampère–Maxwell with displacement current, Biot–Savart geometries, Lorentz force, EM wave speed and intensity) |
| `physicsforge.thermodynamics`    | **Heat laws** (calorimetry, latent heat, Fourier's law of conduction, thermal resistance, Newton's law of cooling, Stefan–Boltzmann, first law, Carnot, ideal gas, thermal expansion) · **Heat equation** (Fourier sine-series solution, FTCS, Crank–Nicolson, high-level `HeatRod`)                                                    |
| `physicsforge.modern`            | **Quantum** (Planck–Einstein, photon momentum, de Broglie, photoelectric effect, Planck's law, Wien, Bohr levels, Heisenberg) · **Relativity** (Lorentz factor, time dilation, length contraction, rest and kinetic energy, momentum)                                                                                           |
| `physicsforge.numerics`          | **Calculus** (derivatives, gradient, **Hessian**, Jacobian, Laplacian, divergence, curl) · **Roots** (Newton–Raphson, bisection, multivariate Newton) · **Optimisation** (damped Newton minimisation with the Hessian) · **Quadrature** (trapezoid, Simpson, adaptive Simpson) · **Transforms** (FFT/IFFT for any length, Fourier series, Laplace and inverse Laplace) |

### Fourier, Laplace and heat transfer

`thermodynamics.heat_equation` puts Fourier's two contributions side by side:

- **Fourier's law** of conduction is $q = -k\,\nabla T$.
- **Fourier series** solve the 1-D heat equation $\partial_t u = \alpha\,\partial_{xx} u$ exactly. Each sine mode
  $\sin(n\pi x/L)$ decays like $e^{-\alpha (n\pi/L)^2 t}$.

`HeatRod` wraps these with Dirichlet boundaries and three methods that can be checked against each other:

- `fourier_solution`: the analytic sine series.
- `solve(method="crank-nicolson")`: unconditionally stable and second order in time.
- `solve(method="explicit")`: FTCS, which rejects time steps beyond the stability limit $r \le 1/2$.

`numerics.laplace_transform` and `numerics.inverse_laplace_transform` (Gaver–Stehfest) cover the usual
transform-domain workflow for linear ODEs.

---

## Architecture

```
               ┌──────────────────────────────────────────────────────────┐
  Python user  │ physicsforge (typed, documented, py.typed)               │  python/physicsforge/
               │   mechanics/  electromagnetism/  thermodynamics/  ...    │
               └────────────────────────────┬─────────────────────────────┘
                                            │ thin wrappers (validation lives in C++)
               ┌────────────────────────────▼─────────────────────────────┐
               │ _physicsforge  (pybind11: py::vectorize, GIL release,    │  python/bindings/
               │                 zero-copy NumPy, exception translation)  │
               └────────────────────────────┬─────────────────────────────┘
                                            │
  C / FFI user ┌──────────────┐   ┌─────────▼─────────────────────────────┐
  ────────────►│ C ABI  pf_*  ├──►│ physicsforge_core (C++20, no deps)    │  include/physicsforge/
               │ (NaN on err) │   │   physicsforge::{math, mechanics,     │  src/
               └──────────────┘   │   electromagnetism, thermodynamics,   │
                src/api/          │   modern, numerics, constants}        │
                                  └───────────────────────────────────────┘
```

```
physicsforge/
├── CMakeLists.txt              # core lib, optional C shared lib, Python module, tests
├── setup.py / pyproject.toml   # CMake-driven wheel build (version read from version.h)
├── build_all.ps1               # multi-version wheel builder (3.12–3.14)
├── include/physicsforge/       # public C++ headers, one per topic, + api.h (C ABI)
│   ├── constants.h  config.h  version.h  export.h  physicsforge.h
│   ├── detail/validation.h
│   ├── math/        mechanics/  electromagnetism/  thermodynamics/  modern/  numerics/
├── src/                        # implementations mirroring include/, + src/api for the C ABI
├── tests/                      # native C++ test suite (ctest)
└── python/
    ├── bindings/               # pybind11 sources, one per domain
    ├── physicsforge/           # the Python package
    └── tests/                  # pytest suite
```

---

## Building and installing

### Requirements

- CMake ≥ 3.20 installed on your system. CLion's bundled CMake is not enough for `setup.py`.
- A C++20 compiler: MSVC 2022 (Desktop development with C++), GCC ≥ 11 or Clang ≥ 14.
- Python ≥ 3.10 with `pip`. pybind11 is pulled in automatically.

### Install from source (recommended)

```bash
git clone https://github.com/OpenForge-Foundation/physicsforge.git
cd physicsforge
pip install .
```

`pip` builds the C++ extension through CMake and installs the package. On Windows machines with several Visual Studio
installations, point CMake at the one that has the C++ workload:

```powershell
$env:CMAKE_GENERATOR_INSTANCE = "C:/Program Files (x86)/Microsoft Visual Studio/2022/BuildTools"
```

### Build wheels

```bash
python -m build                 # sdist + wheel in ./dist
```

For several Python versions at once (Windows):

```powershell
.\build_all.ps1                                   # 3.12, 3.13, 3.14
.\build_all.ps1 -PythonVersion 3.14 -Test -Clean  # one version, run pytest on the wheel, tidy up
```

| Flag             | Effect                                                                    |
|------------------|---------------------------------------------------------------------------|
| `-PythonVersion` | Versions to build (resolved with the `py` launcher).                      |
| `-Test`          | Installs each wheel into its build venv and runs the pytest suite on it. |
| `-Clean`         | Deletes the per-version venvs and intermediate folders afterwards.        |
| `-PipInstall`    | Installs the wheel for that version into `./.venv`.                       |

### Developer build (C++ + in-place module)

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --config Release --parallel
ctest --test-dir build -C Release --output-on-failure
```

| CMake option                       | Default | Purpose                                         |
|------------------------------------|---------|-------------------------------------------------|
| `PHYSICSFORGE_BUILD_TESTS`         | `ON`    | Native test executable registered with CTest.   |
| `PHYSICSFORGE_BUILD_PYTHON`        | `ON`    | The `_physicsforge` pybind11 module.            |
| `PHYSICSFORGE_BUILD_C_SHARED`      | `OFF`   | `physicsforge_c` shared library exposing `pf_*`. |
| `PHYSICSFORGE_WARNINGS_AS_ERRORS`  | `OFF`   | Treat compiler warnings as errors.              |

To work on the Python layer without reinstalling, copy the built `_physicsforge*.pyd`/`.so` into
`python/physicsforge/` and run tests from `python/`.

> **Windows path length:** MSVC's file tracker fails once a path passes 260 characters. If a build in a deeply
> nested folder reports `The directory name is invalid` or a FileTracker error, use a shorter build folder, such as
> `-B C:/b/pf`.

### Using the C++ library

```cmake
add_subdirectory(physicsforge)
target_link_libraries(my_app PRIVATE physicsforge_core)
```

```cpp
#include <physicsforge/physicsforge.h>
namespace pf = physicsforge;

double dv = pf::mechanics::rocket_delta_v(4400.0, 549'000.0, 25'600.0);
auto r = pf::numerics::newton_raphson([](double x) { return x * x - 2.0; }, 1.0);
```

---

## Testing

| Suite              | Command                                              | What it covers                                     |
|--------------------|------------------------------------------------------|----------------------------------------------------|
| C++ unit tests     | `ctest --test-dir build -C Release`                  | Every law, validation path and numerical method; C ABI |
| Python tests       | `pytest` (from the repo root)                        | Physical identities, broadcasting, errors, packaging |
| Doctests           | included in `pytest` via `--doctest-modules`         | Every docstring example is executed                |

Where possible, the tests check physics against itself rather than only against hard-coded numbers. For example:

- Integrating Planck's law gives Stefan–Boltzmann.
- The curl of a gradient is zero.
- Gauss's law and Coulomb's law agree.
- Crank–Nicolson converges at second order.

CI runs the following on every push and pull request:

- clang-format and ruff.
- A C++ build and test matrix on Linux, Windows and macOS.
- Wheel builds tested on Python 3.10, 3.12 and 3.14.

---

## Adding a new law

`rocket_delta_v` is the reference example; the four steps below follow it from C++ to Python. Names below are for a
hypothetical `mechanics.drag_force`.

### 1. C++ declaration and implementation

`include/physicsforge/mechanics/dynamics.h`:

```cpp
// Aerodynamic drag, F = 1/2 rho v^2 C_d A. Requires rho, C_d, A >= 0.
double drag_force(double density, double speed, double drag_coefficient, double area);
```

`src/mechanics/dynamics.cpp`. Validate every argument with the helpers in `detail/validation.h`, which produce
uniform error messages:

```cpp
double drag_force(const double density, const double speed, const double drag_coefficient, const double area) {
    require_non_negative(density, "density");
    require_finite(speed, "speed");
    require_non_negative(drag_coefficient, "drag_coefficient");
    require_non_negative(area, "area");
    return 0.5 * density * speed * speed * drag_coefficient * area;
}
```

A new source file must be added to the `physicsforge_core` sources in `CMakeLists.txt`.

### 2. Native test

Add checks to the matching suite in `tests/` (`tests/test_mechanics.cpp`), covering both values and rejected inputs:

```cpp
PF_CHECK_CLOSE(drag_force(1.2, 10.0, 0.5, 2.0), 60.0);
PF_CHECK_THROWS(drag_force(-1.0, 10.0, 0.5, 2.0));
```

### 3. Python binding

`python/bindings/mechanics.cpp`. Scalar laws are wrapped with `py::vectorize`, so they broadcast:

```cpp
m.def("drag_force", py::vectorize(&pm::drag_force), py::arg("density"), py::arg("speed"),
      py::arg("drag_coefficient"), py::arg("area"));
```

A new binding file needs a `bind_<domain>` declaration and call in `python/bindings/bindings.cpp` (alphabetical), and
an entry in the `pybind11_add_module` list in `CMakeLists.txt`.

### 4. Python API and docstring

`python/physicsforge/mechanics/dynamics.py`. Add a thin typed wrapper with a NumPy-style docstring and a doctest:

```python
def drag_force(density: ArrayLike, speed: ArrayLike, drag_coefficient: ArrayLike, area: ArrayLike) -> FloatOrArray:
    r"""Aerodynamic drag, :math:`F_D = \tfrac12 \rho v^2 C_d A`.

    Parameters
    ----------
    density : float or array_like
        Fluid density :math:`\rho` [kg/m³]. Must be non-negative.
    ...

    Examples
    --------
    >>> drag_force(1.2, 10.0, 0.5, 2.0)
    60.0
    """
    return _core.drag_force(density, speed, drag_coefficient, area)
```

Then export it in the module's `__all__` and in `python/physicsforge/mechanics/__init__.py`. Add a pytest case to
`python/tests/test_mechanics.py`. A physical identity is better than a magic number.

Optionally expose it through the C ABI: declare `pf_mechanics_drag_force` in `include/physicsforge/api.h` and define it
in `src/api/mechanics.cpp` with `guarded(...)`, so errors become `NaN`.

---

## Conventions

| Topic         | Convention                                                                                       |
|---------------|--------------------------------------------------------------------------------------------------|
| Units         | SI everywhere; angles in radians; temperatures in kelvin unless the name says Celsius.          |
| Errors        | C++ `std::invalid_argument` → Python `ValueError`; non-convergence → `ConvergenceError` (carries `.result`); C ABI → `NaN`. |
| Return types  | Scalars in → Python `float` out; arrays in → `numpy.ndarray` out (NumPy broadcasting rules).    |
| Vectors       | Any length-3 sequence or array is accepted wherever a vector is expected; results are `Vector3`. |
| Matrices      | Row-major in C++; `numpy.ndarray` of shape `(n, n)` in Python.                                  |
| Naming        | `snake_case` laws named for the quantity they return (`escape_velocity`, `rocket_delta_v`).    |
| C ABI prefix  | `pf_<domain>_<function>`.                                                                        |
| Docstrings    | NumPy style, formula in the summary line, units in brackets, at least one doctest example.      |

---

## Code formatting and pre-commit hooks

C++ is formatted with `clang-format` (Google-based; see `.clang-format`). Python is linted and formatted with `ruff`.

```bash
pip install pre-commit
pre-commit install          # run automatically on every commit
pre-commit run --all-files  # or run by hand
```

Please run pre-commit before opening a pull request. CI enforces it.

---

## Release notes

### v0.1.0 — Initial pre-release

This first release sets up the architecture, build system and API surface across C++, the C ABI and Python.

**Mechanics:** kinematics, projectile motion, Newton's laws, energy and power, Hooke's law, friction, universal
gravitation (scalar and vector), orbital mechanics, torque, angular momentum, moments of inertia and the Tsiolkovsky
rocket equation.

**Electromagnetism:** Ohm's law, power, voltage and current dividers, series and parallel networks, RC circuits,
Coulomb's law, electric fields and potentials, capacitors, and Maxwell's equations in integral form (Gauss,
Gauss for magnetism, Faraday–Lenz, Ampère–Maxwell), Biot–Savart geometries, the Lorentz force and EM waves.

**Thermodynamics:** calorimetry, Fourier conduction, thermal resistance networks, Newton's cooling, Stefan–Boltzmann
radiation, the first law, Carnot limits, the ideal gas law, and the 1-D heat equation (Fourier series, FTCS,
Crank–Nicolson, `HeatRod`).

**Modern physics:** photon relations, the photoelectric effect, de Broglie, Planck's law, Wien's law, Bohr levels,
Heisenberg uncertainty and special relativity.

**Numerics:** finite-difference derivatives, gradient, Hessian, Jacobian, Laplacian, divergence and curl;
Newton–Raphson, bisection and multivariate Newton; damped Newton minimisation; trapezoid, Simpson and adaptive
quadrature; FFT for any length (radix-2 plus Bluestein); Fourier series; Laplace and inverse Laplace transforms.

**Vectors and constants:** `Vector3` with dot and cross products and geometry helpers; CODATA 2022 constants.

**Quality:**
- Native C++ tests, a pytest suite and executed doctests.
- Warning-free under MSVC `/W4`.
- CI for formatting, a C++ matrix and wheels.

---

## Roadmap

- CUDA / batched back-ends for the vectorised laws and the heat solvers, as in probforge.
- 2-D/3-D heat equation and Neumann/Robin boundary conditions.
- ODE integrators (RK4, adaptive RK45, symplectic Verlet) for N-body and oscillator problems.
- Fluid mechanics (Bernoulli, Reynolds number, drag), waves and optics, and statistical mechanics.
- AC circuit analysis with complex impedance.
- Unit-aware wrappers, possibly with optional `pint` interop.
- Enable warnings-as-errors on GCC/Clang in CI once their `-Wconversion` output is verified clean.
- Published wheels on PyPI (cibuildwheel) and hosted API documentation (Sphinx + numpydoc).
- Benchmarks against NumPy/SciPy equivalents.

---

## Contributors

Special thanks to:

- Manny McGrail
- Zach Pipes

## License

Apache License 2.0. See [LICENSE](LICENSE).

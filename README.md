# physicsforge

## General Information

- This library provides high-performance implementations of foundational physics laws and the numerical methods used
  alongside them, with a C++20 core and NumPy-vectorised Python bindings.
- All quantities are SI (angles in radians, temperatures in kelvin). Physically invalid input raises `ValueError` in
  Python, throws `std::invalid_argument` in C++, and returns `NaN` from the C API.
- Derivations of the laws are widely available and are therefore not duplicated here; each function's docstring (Python)
  and Doxygen comment (C++) gives its formula, units and preconditions.

```python
>>> import physicsforge as pf
>>> pf.mechanics.rocket_delta_v(exhaust_velocity=4400.0, initial_mass=549_000.0, final_mass=25_600.0)
13488.2...
>>> pf.electromagnetism.voltage_divider(12.0, [1e3, 2e3], 2e3)
array([8., 6.])
```

---

## Building and Installing `physicsforge`

### Requirements

- CMake 3.20+ installed on your system (CLion's bundled CMake does not work).
- A C++20 compiler: MSVC 2022 with the C++ workload, GCC 11+ or Clang 14+.
- Python 3.10+.

### Building and Installing the Python Package

From the **project root**:

```bash
pip install .
```

This builds the C++ extension through CMake and installs the package. To build a wheel instead:

```bash
python -m build
pip install .\dist\physicsforge-0.1.0-cpXXX-cpXXX-win_amd64.whl --force-reinstall
```

**Note:** on Windows, MSVC fails on paths longer than 260 characters. If a build reports
`The directory name is invalid`, build from a shorter folder.

### Building the C++ Library and Tests

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --config Release
ctest --test-dir build -C Release --output-on-failure
```

Add `-DPHYSICSFORGE_BUILD_C_SHARED=ON` to also build the C API as a shared library.

---

## Creating A Physics Function

`rocket_delta_v` is a good end-to-end reference. The steps below add a hypothetical `drag_force` to mechanics.

### Writing the C++ Code

1. Declare the function with a Doxygen comment in the domain header, e.g.
   `include/physicsforge/mechanics/dynamics.h`:

```cpp
/**
 * @brief Aerodynamic drag, \f$ F = \tfrac12 \rho v^2 C_d A \f$.
 * @param density Fluid density [kg/m³]; must be non-negative.
 * ...
 * @return Drag force [N].
 * @throws std::invalid_argument If an input is negative.
 */
double drag_force(double density, double speed, double drag_coefficient, double area);
```

2. Define it in the matching source file (`src/mechanics/dynamics.cpp`), validating inputs with the helpers in
   `include/physicsforge/validation.h`.
3. Add checks to the matching test file (`tests/cpp/test_mechanics.cpp`).
4. New source or test files must be added to `CMakeLists.txt`.
5. Optionally expose it through the C API: declare `pf_mechanics_drag_force` in `include/physicsforge/api.h` and define
   it in `src/api/mechanics.cpp`.

### Creating the Python Bindings

1. Add the function to the domain's binding file, `src/bindings/mechanics.cpp`:

```cpp
m.def("drag_force", py::vectorize(&pm::drag_force), py::arg("density"), py::arg("speed"),
      py::arg("drag_coefficient"), py::arg("area"));
```

2. A new binding file needs `void bind_<domain>(py::module_ &m);` and `bind_<domain>(m);` added to
   `src/bindings/bindings.cpp` (in alphabetical order).
3. Regenerate the type stub, which CI checks against the bindings:

```bash
pybind11-stubgen physicsforge._physicsforge -o stubs
cp stubs/physicsforge/_physicsforge.pyi python/physicsforge/
```

### Creating the Python Function

1. Add a typed wrapper with a NumPy-style docstring and a doctest example to
   `python/physicsforge/mechanics/dynamics.py`, and list it in that module's `__all__`.
2. Add a test to `tests/python/test_mechanics.py`.

### Notes

`include/physicsforge/mechanics/rocketry.h` contains the C++ declarations. \
`src/mechanics/rocketry.cpp` contains the C++ definitions. \
`src/bindings/mechanics.cpp` contains the pybind11 bindings. \
`src/bindings/bindings.cpp` contains the module bindings. \
`python/physicsforge/mechanics/rocketry.py` contains the Python functions.

---

## Building Wheels for Multiple Python Versions (3.12–3.14)

1. Install each required Python version: https://www.python.org/downloads/
2. Install CMake: https://cmake.org/download/
3. From the project root, run:

```powershell
.\build_all.ps1
```

4. Optional flags:
    - `-PythonVersion`: build for specific versions only (e.g. `-PythonVersion 3.13`).
    - `-Test`: run the pytest suite against each built wheel.
    - `-Clean`: remove temporary build artifacts and virtual environments afterwards.
    - `-PipInstall`: install the built wheel for that version into `.venv` (e.g. `-PipInstall 3.14`).

---

## Code Formatting and Pre-Commit Hooks

This repository enforces consistent formatting using `clang-format` for C++ and `ruff` for Python.

```bash
pip install pre-commit
pre-commit run --all-files
```

**Note:** running pre-commit is **required** before submitting commits or pull requests.

---

## Release Notes

### v0.1.0 — Initial Pre-Release

The first pre-release of physicsforge, establishing the architecture, API surface and build system across C++, the C
API and Python.

Mechanics:

- Kinematics and projectile motion
- Newton's laws, momentum, work, energy, power, Hooke's law and friction
- Universal gravitation (scalar and vector), escape and orbital velocity, Kepler's third law
- Torque (scalar and vector), angular momentum, moments of inertia, parallel-axis theorem
- Tsiolkovsky rocket equation, specific impulse and thrust

Electromagnetism:

- Ohm's law, power, voltage and current dividers, series/parallel networks, RC transients
- Coulomb's law, electric fields and potentials, capacitors
- Maxwell's equations (Gauss, Faraday, Ampère–Maxwell), Biot–Savart geometries, Lorentz force, EM waves

Thermodynamics:

- Calorimetry, Fourier conduction, thermal resistance, Newton's cooling, Stefan–Boltzmann radiation
- First law, Carnot limits, ideal gas law, thermal expansion
- 1-D heat equation: Fourier-series solution, FTCS and Crank–Nicolson, wrapped by `HeatRod`

Quantum:

- Planck–Einstein relation, photon momentum, de Broglie wavelength, photoelectric effect
- Planck's law, Wien's law, Bohr energy levels, Heisenberg uncertainty

Relativity:

- Lorentz factor, time dilation, length contraction, rest energy, relativistic momentum and kinetic energy

Numerics:

- Derivatives, gradient, Hessian, Jacobian, Laplacian, divergence and curl
- Newton–Raphson, bisection, multivariate Newton, Newton minimisation with the Hessian
- ODE integrators: Euler, semi-implicit (symplectic) Euler and Runge–Kutta 4
- Trapezoid, Simpson and adaptive quadrature
- FFT, Fourier series, Laplace and inverse Laplace transforms

Vectors and Constants:

- `Vector3` with dot and cross products; CODATA 2022 constants including Planck's constant, e and π

Testing and CI:

- C++ unit tests, a pytest suite and executed doctests
- GitHub Actions for formatting, a per-push build and test check, type-stub verification and wheel builds on `main`

---

Long-term plans:

- Add CUDA/batched back-ends for the vectorised laws and heat solvers, as in probforge
- Add adaptive ODE integration (RK45) and velocity Verlet
- Extend the heat equation to 2-D/3-D and Neumann/Robin boundaries
- Add fluid mechanics, waves and optics, and AC circuit analysis
- Publish wheels to PyPI and generate hosted API documentation from the Doxygen comments and docstrings
- Create benchmarks against NumPy/SciPy

---

## Contributors

Special thanks to:

- Manny McGrail
- Zach Pipes

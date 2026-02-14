# Copilot Instructions

## Overview

Astrodynamics library for the Circular Restricted Three-Body Problem (CR3BP) and related orbital mechanics models. Author: Abram Aguilar (NASA JSC).

## Project structure

```text
astro.py              # Core: CelestialObject class, CR3BP ODEs, Lagrange points, Jacobi constant,
                      # differential correctors (solveLyapunov, solveHalo), continuation methods
general_utils.py      # Unit conversions, constants (EARTH_MOON_CR3BP_*), KDTree proximity search,
                      # osculating elements, JPL 3-Body API client (Jpl3BodyData)
plot_utils.py         # Plotly-based 3D trajectory/Poincare plotting
asset_utils.py        # ASSET optimal control solver wrappers (Jet parallel jobs)
data_utils.py         # Pickle serialization (writeDataToPkl / readDataFromPkl)
Astro.jl              # Julia port of astro.py (standalone module, uses DifferentialEquations.jl)

Models/
  CR3BP.py            # ASSET-format CR3BP ODE (CR3BPFrame, CR3BP classes)
  BCR4BP.py           # Bicircular Restricted 4-Body Problem (Sun perturbation)
  CR3BP_LT.py         # CR3BP with low-thrust (mass + control variables)
  CR3BP_ST.py         # CR3BP with Sundman Transform (regularization near P2)
  EPPR_ST.py          # Ephemeris Pulsating Rotating frame with Sundman Transform
  LVLH_CR3BP.py       # Local Vertical Local Horizontal frame for RPO

CR3BP_Utils/
  stm.py              # State Transition Matrix computation (serial & parallel)
  manifolds.py        # Stable/unstable manifold computation along periodic orbits
  events.py           # ASSET integrator event functions (axis crossings, distance triggers)
  solve_periodic_orbit.py  # Generic periodic orbit solver using ASSET phase optimization
  constraints.py      # ASSET equality/inequality constraints (position matching, keep-out zones)
  objectives.py       # ASSET objective functions (delta-V minimization, low-thrust optimal)
  control_laws.py     # VNB-frame thrust direction control laws
  ftle.py             # Finite-Time Lyapunov Exponent computation

examples/             # Example scripts (Python + Julia NRHO solvers)
data/                 # Pickle orbit files, periodicLagrangeOrbits.csv
```

## Languages and dependencies

**Python** (primary): numpy, scipy, pandas, matplotlib, plotly, requests, asset_asrl. No package manager or virtual environment configured.

**Julia** (secondary): `Astro.jl` is a standalone module. Dependencies: LinearAlgebra, StaticArrays, DifferentialEquations, Roots, CSV, DataFrames, Plots.

The `asset_asrl` package is required by all files in `Models/`, `CR3BP_Utils/`, `asset_utils.py`, and parts of `general_utils.py`. The core `astro.py` does NOT depend on ASSET (uses scipy directly).

## Running code

No test suite, linter, or CI/CD exists. Scripts are run directly:

```bash
python examples/Solve_NRHOs_CR3BP.py
julia examples/Solve_NRHOs_CR3BP.jl
```

Example scripts use `sys.path.append` to find the project root. Imports reference modules by relative path from project root (e.g., `from Models.CR3BP import CR3BP`).

## Naming conventions

Python:
- Classes: `UpperCamelCase` (CelestialObject, CR3BPFrame, Jpl3BodyData)
- Class methods: `lower_case_with_underscores` (equatorial_radius, orbit_semi_major_axis)
- Standalone functions: `camelCase` (characteristicQuantities, lagrangePoints, solveHalo)
- Module-level constants: `UPPER_SNAKE_CASE` (EARTH_MOON_CR3BP_MU, M_TO_KM)
- ASSET alias convention at top of files: `vf = ast.VectorFunctions`, `oc = ast.OptimalControl`, `Args = vf.Arguments`

Julia:
- Types: `CamelCase` (CelestialObject, BodyData)
- Functions: `snake_case` (characteristic_quantities, solve_halo)
- In-place functions: `snake_case!` (cr3bp_ode!, cr3bp_ode_stm!)

## Key domain concepts

All dynamics use **non-dimensional** coordinates. The CR3BP state vector is `[x, y, z, vx, vy, vz]` (6 elements). ASSET-format states append time as element 7: `[x, y, z, vx, vy, vz, t]`.

- **mu (mass ratio)**: mu = m2 / (m1 + m2), always between 0 and 0.5
- **lstar**: characteristic length (distance between primaries, km)
- **tstar**: characteristic time (seconds), derived as sqrt(lstar^3 / (G * m_star))
- **vstar**: characteristic velocity = lstar / tstar
- **P1**: primary body at [-mu, 0, 0]; **P2**: secondary body at [1-mu, 0, 0]
- **L1-L5**: Lagrange equilibrium points computed by each Frame class
- **STM**: 6x6 State Transition Matrix; when appended to state gives a 42-element vector
- **Jacobi constant**: energy-like integral of motion; computed by `jacobi()` / `jacobi_constant()`
- **Sundman Transform (ST)**: time regularization that scales integration variable by distance to P2

## Docstring style

NumPy-style docstrings with Parameters/Returns sections. Each file starts with a module-level docstring listing all functions and classes. Follow this pattern when adding new code.

## Common patterns

**Creating a CR3BP system:**

```python
from astro import CelestialObject, characteristicQuantities
from Models.CR3BP import CR3BP

earth = CelestialObject("Earth")
moon = CelestialObject("Moon")
cq = characteristicQuantities(earth, moon)
ode = CR3BP(cq["mu"], cq["lstar"], cq["tstar"])
```

**ASSET ODE class pattern** (Models/ directory): each model defines a Frame class with equations of motion, then an ODE class that inherits from both `ODEBase` and the Frame class. The `__init__` calls `oc.ODEArguments(n_states, n_controls)`, builds a VectorFunction, and passes it to `ODEBase.__init__`.

**Differential correction**: `solveLyapunov` and `solveHalo` in `astro.py` use scipy's `solve_ivp` with STM propagation. The ASSET-based `solvePeriodic` in `CR3BP_Utils/solve_periodic_orbit.py` uses ASSET's phase optimization.

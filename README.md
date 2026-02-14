# Astro Tools

Astrodynamics library for the Circular Restricted Three-Body Problem (CR3BP) and related orbital mechanics models.

**Author:** Abram Aguilar (NASA Johnson Space Center)

## Overview

The **Circular Restricted Three-Body Problem (CR3BP)** is a classical problem in celestial mechanics that models the motion of a massless particle (spacecraft) under the gravitational influence of two primary bodies (e.g., Earth and Moon). This formulation assumes the two primaries orbit their common barycenter in circular orbits, and the spacecraft's mass is negligible such that it doesn't affect the primaries' motion.

Astro Tools provides a comprehensive suite of numerical methods and models for:
- **Cislunar space trajectory design** (Earth-Moon system dynamics)
- **Orbital mechanics research** (periodic orbits, manifolds, stability analysis)
- **Optimal control problems** (low-thrust trajectory optimization)
- **Multi-body dynamics** (CR3BP extensions including solar perturbations)

This library offers a **dual-language implementation**:
- **Python** (primary): Full-featured implementation with integration to NASA's ASSET optimal control library
- **Julia** (secondary): Standalone module (`Astro.jl`) for high-performance differential equation solving

## Key Features

- **CR3BP Dynamics and Lagrange Points**: Compute equations of motion in the rotating reference frame, find L1-L5 equilibrium points, and analyze their stability properties

- **Multiple Orbital Mechanics Models**:
  - **BCR4BP** (Bicircular Restricted 4-Body Problem): Includes solar perturbations on cislunar dynamics
  - **CR3BP_LT** (Low-Thrust CR3BP): Extended state with mass and control variables for electric propulsion
  - **CR3BP_ST** (Sundman Transform): Time regularization for trajectories near the secondary body
  - **EPPR_ST** (Ephemeris Pulsating Rotating with Sundman Transform): High-fidelity ephemeris-based dynamics
  - **LVLH_CR3BP** (Local Vertical Local Horizontal): Rendezvous and proximity operations frame

- **Periodic Orbit Solving and Differential Correction**: Algorithms to find and refine periodic orbits (Lyapunov, Halo, NRHO families) using State Transition Matrix-based corrections

- **State Transition Matrix (STM) Computation**: Serial and parallel computation of linearized dynamics for sensitivity analysis and targeting

- **Stable/Unstable Manifold Analysis**: Compute invariant manifolds along periodic orbits for designing low-energy transfers and understanding transport dynamics

- **Optimal Control Integration**: Seamless integration with NASA's ASSET (Astrogator Software for Exploring Trajectories) library for trajectory optimization problems

- **Finite-Time Lyapunov Exponent (FTLE) Analysis**: Identify Lagrangian Coherent Structures (LCS) that govern transport in phase space

- **3D Trajectory Visualization**: Plotly-based interactive 3D plotting for visualizing orbits, manifolds, and phase space structures

- **Poincaré Section Plotting**: Generate Poincaré maps to analyze orbit families and bifurcations

- **Continuation Methods**: Trace families of periodic orbits as system parameters vary

## Quick Start

### Creating a CR3BP System

```python
from astro import CelestialObject, characteristicQuantities
from Models.CR3BP import CR3BP

# Define the two primary bodies
earth = CelestialObject("Earth")
moon = CelestialObject("Moon")

# Compute characteristic quantities (mu, lstar, tstar, vstar)
cq = characteristicQuantities(earth, moon)

# Create the CR3BP equations of motion
ode = CR3BP(cq["mu"], cq["lstar"], cq["tstar"])
```

### Computing Lagrange Points

```python
from astro import lagrangePoints

# Compute all five Lagrange points (L1-L5)
L_points = lagrangePoints(ode.mu)  # Returns positions in rotating frame

print(f"L1 position: {L_points[0]}")
print(f"L2 position: {L_points[1]}")
```

For complete working examples, see:
- `examples/Solve_NRHOs_CR3BP.py` - Solves 9:2 and 4:1 L2 Southern Near Rectilinear Halo Orbits
- `examples/Solve_NRHOs_CR3BP.jl` - Julia implementation of NRHO solver

## Project Structure

```
astro.py              # Core: CelestialObject class, CR3BP ODEs, Lagrange points, Jacobi constant,
                      #   differential correctors (solveLyapunov, solveHalo), continuation methods
general_utils.py      # Unit conversions, constants (EARTH_MOON_CR3BP_*), KDTree proximity search,
                      #   osculating elements, JPL 3-Body API client (Jpl3BodyData)
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

## Installation and Dependencies

### Python Dependencies

**Core scientific computing:**
- numpy
- scipy
- pandas
- matplotlib
- plotly
- requests

**Optimal control (required for Models/ and CR3BP_Utils/):**
- asset_asrl (NASA ASSET library)

**Note:** The core module `astro.py` does NOT depend on ASSET and uses scipy directly for integration and root-finding.

### Julia Dependencies (for Astro.jl)

- LinearAlgebra
- StaticArrays
- DifferentialEquations
- Roots
- CSV
- DataFrames
- Plots

### Installation Notes

- **No package manager or virtual environment is currently configured**
- Scripts are run directly: `python examples/Solve_NRHOs_CR3BP.py`
- Example scripts use `sys.path.append` to locate the project root
- Imports reference modules by relative path from project root (e.g., `from Models.CR3BP import CR3BP`)

## Key Technical Concepts

### Non-dimensional Coordinates

All dynamics in this library use **non-dimensional coordinates** to improve numerical stability and enable dimensionless analysis across different three-body systems.

**State Vector Formats:**
- Standard CR3BP: `[x, y, z, vx, vy, vz]` (6 elements)
- ASSET format: `[x, y, z, vx, vy, vz, t]` (7 elements, with time appended)
- STM-augmented: `[x, y, z, vx, vy, vz, φ₁₁, φ₁₂, ..., φ₆₆]` (42 elements)

### Mass Ratio (μ)

The **mass ratio** μ defines the three-body system:

```
μ = m₂ / (m₁ + m₂)
```

where:
- m₁ = mass of the primary body (larger mass)
- m₂ = mass of the secondary body (smaller mass)
- μ is always between 0 and 0.5

For the Earth-Moon system: μ ≈ 0.01215

### Characteristic Quantities

Characteristic quantities provide the scaling factors between dimensional and non-dimensional coordinates:

- **lstar** (characteristic length): Distance between the two primaries, measured in kilometers. For Earth-Moon, lstar = 384,400 km.

- **tstar** (characteristic time): Computed as:
  ```
  tstar = √(lstar³ / (G · mstar))
  ```
  where G is the gravitational constant and mstar = m₁ + m₂. Measured in seconds.

- **vstar** (characteristic velocity):
  ```
  vstar = lstar / tstar
  ```
  Measured in km/s.

These allow conversion between systems:
- **Non-dimensional → Dimensional**: multiply by characteristic quantity
- **Dimensional → Non-dimensional**: divide by characteristic quantity

### Primaries and Reference Frame

The CR3BP uses a **rotating reference frame** centered at the barycenter of the two primary bodies:

- **P1** (primary body): Located at `[-μ, 0, 0]`
- **P2** (secondary body): Located at `[1-μ, 0, 0]`

**Frame properties:**
- Origin at the system barycenter
- x-axis connects P1 and P2 (through both primaries)
- z-axis is perpendicular to the orbital plane (angular momentum direction)
- y-axis completes the right-handed coordinate system
- Frame rotates with angular velocity ω = 1 (non-dimensional)

### Lagrange Points (L1-L5)

The **Lagrange points** (also called libration points) are five equilibrium positions in the rotating frame where the gravitational forces of the two primaries and the centrifugal force balance exactly.

**Collinear Points (L1, L2, L3):**
- **L1**: Between P1 and P2, approximately 326,000 km from Earth (for Earth-Moon system)
- **L2**: Beyond P2, approximately 449,000 km from Earth
- **L3**: Beyond P1, on the opposite side from P2
- All collinear points are **unstable** (saddle × center × center dynamics)

**Triangular Points (L4, L5):**
- **L4**: Forms an equilateral triangle with P1 and P2 (leading the secondary body by 60°)
- **L5**: Forms an equilateral triangle with P1 and P2 (trailing the secondary body by 60°)
- These points are **stable** for sufficiently small mass ratios (μ < μcrit ≈ 0.0385)
- For Earth-Moon (μ ≈ 0.01215), L4 and L5 are linearly stable

**Applications:**
- Gateway lunar outpost is planned for a Near Rectilinear Halo Orbit (NRHO) around L2
- James Webb Space Telescope orbits L2 of the Sun-Earth system
- Trojan asteroids reside near L4 and L5 of the Sun-Jupiter system

### State Transition Matrix (STM)

The **State Transition Matrix** Φ(t, t₀) is a 6×6 matrix that describes how small perturbations in initial conditions evolve over time in the linearized dynamics:

```
δx(t) = Φ(t, t₀) · δx(t₀)
```

where δx is a small perturbation in the state vector.

**Properties:**
- Computed by integrating the variational equations alongside the trajectory
- Initial condition: Φ(t₀, t₀) = I₆ (6×6 identity matrix)
- When appended to the state vector, creates a 42-element augmented state: `[x, y, z, vx, vy, vz, φ₁₁, φ₁₂, ..., φ₆₆]`

**Applications:**
- **Differential correction**: Computing Newton-Raphson updates to find periodic orbits
- **Targeting**: Designing maneuvers to reach specific target states
- **Stability analysis**: Eigenvalues of Φ(T, 0) determine orbit stability
- **Manifold computation**: Eigenvectors define stable/unstable directions

### Jacobi Constant

The **Jacobi constant** C is an energy-like integral of motion in the CR3BP. It remains constant along any trajectory and is defined as:

```
C = -2(v² - U(x, y, z))
```

where:
- v² = vx² + vy² + vz² (velocity magnitude squared)
- U is the effective potential (pseudopotential):
  ```
  U = (x² + y²)/2 + (1-μ)/r₁ + μ/r₂ + μ(1-μ)/2
  ```
- r₁ = distance to P1
- r₂ = distance to P2

**Physical interpretation:**
- Lower C → higher energy (spacecraft can access more regions)
- Higher C → lower energy (more restricted motion)
- Zero-velocity curves define **forbidden regions** where motion is not possible
- Periodic orbits have constant C values

**Applications:**
- Quick check for numerical integration accuracy (C should remain constant)
- Classifying orbit families
- Understanding reachable regions of space

### Differential Correction

**Differential correction** is an iterative Newton-Raphson algorithm used to find periodic orbits with high precision. The method leverages the State Transition Matrix to compute how changes in initial conditions affect the periodicity constraints.

**Basic algorithm for symmetric periodic orbits:**
1. Start with an initial guess x₀ that is approximately periodic
2. Propagate x₀ with STM for half-period T/2 to reach x(T/2)
3. For symmetric orbits crossing x-z plane, enforce: y(T/2) = 0, vx(T/2) = 0, vz(T/2) = 0
4. Use STM to compute corrections: Δx₀ = -[∂f/∂x₀]⁻¹ · f(x₀)
5. Update: x₀ ← x₀ + Δx₀
6. Repeat until convergence (|f(x₀)| < tolerance)

**Implementations in this library:**
- `solveLyapunov(x0, T, mu)` in `astro.py`: Scipy-based solver for planar Lyapunov orbits
- `solveHalo(x0, T, mu)` in `astro.py`: Scipy-based solver for 3D Halo orbits
- `solvePeriodic(x0, T, ode, fix_init)` in `CR3BP_Utils/solve_periodic_orbit.py`: ASSET-based general periodic orbit solver with flexible constraints

**Advantages:**
- Extremely fast convergence (quadratic near solution)
- Can refine approximate solutions to machine precision
- Enables tracing orbit families via continuation

### Sundman Transform (ST)

The **Sundman Transform** is a time regularization technique that addresses numerical difficulties when integrating trajectories that pass close to a primary body (especially P2).

**Concept:**
Replace physical time t with a fictitious time τ defined by:
```
dτ = dt / r
```
where r is the distance to P2.

**Effect:**
- Near P2 (small r): Physical time advances slowly relative to fictitious time
- Far from P2 (large r): Physical time advances quickly relative to fictitious time
- This "stretches out" the integration near close approaches, improving numerical stability

**Trade-offs:**
- **Advantage**: Prevents step size issues and improves accuracy near P2
- **Disadvantage**: Equations become more complex; additional state variable needed
- **Use case**: Essential for low-altitude lunar orbits or trajectories with close lunar flybys

**Implementation:** See `Models/CR3BP_ST.py` for the CR3BP formulation with Sundman Transform.

### Manifolds

**Invariant manifolds** are sets of trajectories that asymptotically approach (stable manifold) or depart from (unstable manifold) a periodic orbit.

**Mathematical definition:**
- **Stable manifold** Wˢ: All trajectories that approach the periodic orbit as t → +∞
- **Unstable manifold** Wᵘ: All trajectories that approach the periodic orbit as t → -∞

**Computation:**
1. Compute the **monodromy matrix** M = Φ(T, 0) for one complete orbit period T
2. Find eigenvalues λ and eigenvectors v of M
3. For unstable manifolds: Use eigenvectors corresponding to |λ| > 1
4. For stable manifolds: Use eigenvectors corresponding to |λ| < 1
5. Perturb initial conditions along eigenvector directions: x₀ ± ε·v
6. Propagate forward (unstable) or backward (stable) in time

**Applications:**
- **Low-energy transfers**: Design trajectories that "ride" manifolds between periodic orbits
- **Interplanetary transport**: Manifolds of L1/L2 orbits create "tubes" connecting Earth and Moon
- **Mission design**: Gateway station transfers, ballistic lunar transfers, asteroid tour missions
- **Chaos and transport**: Understanding how trajectories move between different regions of phase space

**Implementation:** See `CR3BP_Utils/manifolds.py` for manifold computation along periodic orbits.

### Finite-Time Lyapunov Exponent (FTLE)

The **Finite-Time Lyapunov Exponent** measures the maximum rate of separation of initially nearby trajectories over a finite time interval T.

**Definition:**
```
FTLE(x₀, T) = (1/|T|) · ln(σmax(Φ(T, 0)))
```
where σmax is the largest singular value of the flow map Φ.

**Interpretation:**
- **High FTLE** (ridges in FTLE field): Indicates strong stretching → identifies Lagrangian Coherent Structures (LCS)
- **Low FTLE**: Indicates weak stretching → more predictable motion
- **LCS**: Act as "transport barriers" that organize phase space into distinct dynamical regions

**Applications:**
- **Mission planning**: Identify regions with sensitive dynamics that require precise navigation
- **Trajectory design**: Exploit or avoid transport barriers depending on mission objectives
- **Space situational awareness**: Predict long-term evolution of debris or natural objects
- **Scientific discovery**: Reveal hidden structure in complex dynamical systems

**Computation:**
1. Create a grid of initial conditions in phase space
2. For each initial condition, propagate forward (or backward) for time T with STM
3. Extract the flow map Φ(T, 0) and compute largest singular value
4. Apply FTLE formula
5. Visualize as a scalar field (color map or contour plot)

**Implementation:** See `CR3BP_Utils/ftle.py` for FTLE computation utilities.

## Code Organization and Conventions

### Python Naming Conventions

- **Classes**: `UpperCamelCase`
  - Examples: `CelestialObject`, `CR3BPFrame`, `Jpl3BodyData`

- **Class methods**: `lower_case_with_underscores`
  - Examples: `equatorial_radius()`, `orbit_semi_major_axis()`

- **Standalone functions**: `camelCase`
  - Examples: `characteristicQuantities()`, `lagrangePoints()`, `solveHalo()`

- **Module-level constants**: `UPPER_SNAKE_CASE`
  - Examples: `EARTH_MOON_CR3BP_MU`, `M_TO_KM`, `EARTH_MOON_CR3BP_DAY`

- **ASSET alias convention** (top of files in Models/ and CR3BP_Utils/):
  ```python
  import asset_asrl as ast
  vf = ast.VectorFunctions
  oc = ast.OptimalControl
  Args = vf.Arguments
  ```

### Julia Naming Conventions

- **Types**: `CamelCase`
  - Examples: `CelestialObject`, `BodyData`

- **Functions**: `snake_case`
  - Examples: `characteristic_quantities()`, `solve_halo()`, `lagrange_points()`

- **In-place functions**: `snake_case!` (with exclamation mark)
  - Examples: `cr3bp_ode!()`, `cr3bp_ode_stm!()`

### Docstring Style

- **NumPy-style docstrings** with `Parameters` and `Returns` sections
- Each file begins with a **module-level docstring** listing all functions and classes
- Follow this pattern when adding new code

**Example:**
```python
def solveHalo(x0, T, mu, tol=1e-12, max_iter=100):
    """
    Solve for a periodic Halo orbit using differential correction.

    Parameters
    ----------
    x0 : array_like, shape (6,)
        Initial guess for state vector [x, y, z, vx, vy, vz]
    T : float
        Initial guess for orbital period (non-dimensional time)
    mu : float
        CR3BP mass ratio
    tol : float, optional
        Convergence tolerance (default: 1e-12)
    max_iter : int, optional
        Maximum number of iterations (default: 100)

    Returns
    -------
    x0_corrected : ndarray, shape (6,)
        Corrected initial state for periodic orbit
    T_corrected : float
        Corrected orbital period
    converged : bool
        True if differential correction converged
    """
```

## Usage Examples

This library includes working example scripts that demonstrate key capabilities:

### Solving Near Rectilinear Halo Orbits (NRHOs)

**Python:** `examples/Solve_NRHOs_CR3BP.py`

Solves two lunar NRHOs (9:2 and 4:1 L2 Southern families) using ASSET-based periodic orbit solver. Initial conditions are taken from Dr. Emily Spreen's PhD dissertation. The script:
1. Defines the Earth-Moon CR3BP system
2. Sets up initial guesses from literature
3. Applies differential correction with ASSET optimizer
4. Validates results against published Jacobi constant and period values
5. Saves converged trajectories to pickle files

**Julia:** `examples/Solve_NRHOs_CR3BP.jl`

Julia implementation of the same NRHO solver, demonstrating the `Astro.jl` standalone module.

### Running Examples

```bash
# Python version
python examples/Solve_NRHOs_CR3BP.py

# Julia version
julia examples/Solve_NRHOs_CR3BP.jl
```

**Note:** Python examples use `sys.path.append(str(pathlib.Path(__file__).parents[1].absolute()))` to locate the project root, so they can be run from any directory.

### Computed Data

The `data/` directory contains:
- **Pickle files** (`.pkl`): Saved trajectory data from example scripts
- **`periodicLagrangeOrbits.csv`**: Tabulated periodic orbit families with initial conditions and properties

## Additional Documentation

For comprehensive developer documentation, see **`CLAUDE.md`**, which includes:
- **Detailed common patterns**: Step-by-step recipes for typical tasks
- **ASSET ODE class design patterns**: How to create new dynamical models
- **Complete function and class listings**: Full API reference for each module
- **Integration instructions for coding agents**: Guidelines for AI assistants working with this codebase

The `CLAUDE.md` file serves as the detailed technical reference, while this README provides an accessible entry point for new users.

## Disclaimer

**Note:** Some code in this repository may be AI-generated. This repository serves as an educational resource for learning how to work with coding agents for astrodynamics research and development.

**Original Author:** Abram Aguilar (NASA Johnson Space Center)

**Purpose:** This repository is intended for learning, research, and educational purposes. Users should validate results independently for production applications or mission-critical systems.

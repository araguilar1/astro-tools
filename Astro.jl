"""
    module Astro

Julia module for astrodynamic applications, focusing on the Circular Restricted
Three-Body Problem (CR3BP). Provides celestial body data, CR3BP equations of motion,
Lagrange point computation, periodic orbit solvers (Lyapunov/Halo), manifold analysis
utilities, and zero-velocity curve plotting.

Converted from astro.py by Abram Aguilar (original: 12/29/2019).

Naming conventions:
    Types:     CamelCase        (e.g. CelestialObject)
    Functions: snake_case       (e.g. characteristic_quantities)
    Mutating:  snake_case!      (e.g. cr3bp_ode!)
"""
module Astro

using LinearAlgebra
using Printf
using StaticArrays
using DifferentialEquations
using Roots
using CSV, DataFrames
using Plots

export CelestialObject,
       axial_rotation_period, equatorial_radius, mu, mass,
       orbit_semi_major_axis, orbital_period, orbit_eccentricity,
       orbit_inclination, body_name,
       characteristic_quantities,
       cr3bp_ode!, cr3bp_ode_stm!,
       u_star, jacobi_constant,
       lagrange_points, lagrange_points_roots,
       non_dim, dimensionalize,
       cr3bp_to_j2k,
       periodic_orbit_df, plot_zvc,
       two_body_ode!,
       lagrange_point_eigendecomp,
       solve_lyapunov, lyapunov_continuation,
       solve_halo, halo_continuation

# =============================================================================
# Constants
# =============================================================================

"""Gravitational constant in km³/(kg·s²)"""
const G_KM3 = 6.6743e-11 / 1e9  # 6.6743e-20

# =============================================================================
# CelestialObject
# =============================================================================

"""
Data for a single celestial body.

Fields with `Nothing` indicate unavailable data (e.g. the Sun has no orbital elements).
Bodies with synchronous rotation store `nothing` for `axial_rotation_period`.
"""
struct BodyData
    axial_rotation_period::Union{Float64, Nothing}  # Rev/Day ("Synch" → nothing)
    equatorial_radius::Float64                       # km
    mu_grav::Float64                                 # km³/s²
    orbit_semi_major_axis::Union{Float64, Nothing}   # km ("--" → nothing)
    orbital_period_days::Union{Float64, Nothing}     # days
    orbit_eccentricity::Union{Float64, Nothing}
    orbit_inclination::Union{Float64, Nothing}       # deg
    name::String
end

const BODY_DATABASE = BodyData[
    #                         rot_period    eq_radius   mu               sma              period     ecc        inc        name
    BodyData(                 0.039401,     695700.0,   132712440041.94, nothing,         nothing,   nothing,   nothing,   "Sun"),
    BodyData(                 0.0366,       1738.0,     4902.8,          384400.0,        27.32,     0.0549,    5.145,     "Moon"),
    BodyData(                 0.017051,     2440.0,     22031.78,        4902.8,          87.97,     0.205647,  7.0039,    "Mercury"),
    BodyData(                 0.004115,     6051.89,    324858.59,       108208441.28,    224.7,     0.006794,  3.39449,   "Venus"),
    BodyData(                 1.002737,     6378.14,    398600.44,       149657721.28,    365.26,    0.016192,  0.0045,    "Earth"),
    BodyData(                 0.974700,     3394.0,     42828.38,        227937168.37,    686.98,    0.093343,  1.84814,   "Mars"),
    BodyData(                 2.418111,     71492.0,    126686534.91,    778350260.22,    4332.59,   0.048708,  1.30394,   "Jupiter"),
    BodyData(                 2.252205,     60268.0,    37931207.80,     1433203470.67,   10755.70,  0.050663,  2.48560,   "Saturn"),
    BodyData(                -1.3921089,    25559.0,    5793951.32,      2863429878.70,   30685.40,  0.048551,  0.77151,   "Uranus"),
    BodyData(                 1.489754,     24766.0,    6835099.50,      4501859020.15,   60189.0,   0.007183,  1.77740,   "Neptune"),
    BodyData(                 0.156563,     1188.30,    869.34,          6018076570.89,   91101.50,  0.258313,  17.22524,  "Pluto"),
    BodyData(                 0.15625,      605.0,      102.27,          19596.84,        6.39,      0.00005,   112.89596, "Charon"),
    BodyData(                 0.5291,       43.33,      0.0003,          48690.0,         24.85,     0.238214,  112.88839, "Nix"),
    BodyData(                 2.328,        65.0,       0.000320,        64738.0,         38.2,      0.0058652, 0.24200,   "Hydra"),
    BodyData(                 nothing,      2634.0,     9891.0,          1070042.8,       7.15,      0.0006,    0.186,     "Ganymede"),
    BodyData(                 nothing,      2575.5,     8978.13,         1221870.0,       15.95,     0.0288,    0.28,      "Titan"),
    BodyData(                 nothing,      788.9,      235.4,           435800.0,        8.71,      0.0022,    0.1,       "Titania"),
    BodyData(                 2.644860,     469.7,      62.63,           413968739.37,    1680.22,   0.076103,  10.6007,   "Ceres"),
    BodyData(                 nothing,      252.30,     1.21135,         238040.0,        1.370218,  0.0047,    0.009,     "Enceladus"),
    BodyData(                 nothing,      13.10,      0.000721,        9377.20,         0.32,      0.0151,    1.082,     "Phobos"),
    BodyData(                 nothing,      1352.60,    1432.93,         354760.0,        5.876854,  0.000016,  156.834,   "Triton"),
    BodyData(                 0.05992,      2403.0,     7181.32,         1883000.0,       16.69,     0.007,     0.281,     "Callisto"),
    BodyData(                 4.2625,       151.959,    3203.56,         671100.0,        3.552,     0.0094,    7.483,     "Europa"),
]

const NAME_TO_ID = Dict{String, Int}(
    "sun"       => 1,  "moon"      => 2,  "mercury"   => 3,  "venus"    => 4,
    "earth"     => 5,  "mars"      => 6,  "jupiter"   => 7,  "saturn"   => 8,
    "uranus"    => 9,  "neptune"   => 10, "pluto"     => 11, "charon"   => 12,
    "nix"       => 13, "hydra"     => 14, "ganymede"  => 15, "titan"    => 16,
    "titania"   => 17, "ceres"     => 18, "enceladus" => 19, "phobos"   => 20,
    "triton"    => 21, "callisto"  => 22, "europa"    => 23,
)

"""
    CelestialObject(name::String)

Celestial body with physical and orbital data. Case-insensitive name lookup
for 23 bodies (Sun, Moon, 8 planets, Pluto, and assorted moons).
"""
struct CelestialObject
    body::String
    id::Int
    data::BodyData

    function CelestialObject(name::String)
        key = lowercase(name)
        if !haskey(NAME_TO_ID, key)
            error("Celestial object not found: $name")
        end
        id = NAME_TO_ID[key]
        new(name, id, BODY_DATABASE[id])
    end
end

axial_rotation_period(obj::CelestialObject) = obj.data.axial_rotation_period
equatorial_radius(obj::CelestialObject)     = obj.data.equatorial_radius
mu(obj::CelestialObject)                    = obj.data.mu_grav
mass(obj::CelestialObject)                  = obj.data.mu_grav / G_KM3
orbit_semi_major_axis(obj::CelestialObject) = obj.data.orbit_semi_major_axis
orbit_eccentricity(obj::CelestialObject)    = obj.data.orbit_eccentricity
orbit_inclination(obj::CelestialObject)     = obj.data.orbit_inclination
body_name(obj::CelestialObject)             = obj.data.name

function orbital_period(obj::CelestialObject; unit::String="day")
    p = obj.data.orbital_period_days
    if p === nothing
        return nothing
    end
    if unit == "day"
        return p
    elseif unit == "year"
        return p / 365.25
    else
        @warn "Unit not recognized, defaulting to days."
        return p
    end
end

# =============================================================================
# CR3BP Characteristic Quantities
# =============================================================================

"""
    characteristic_quantities(P1::CelestialObject, P2::CelestialObject)

Compute the CR3BP non-dimensional parameters for a primary–secondary system.
Returns a NamedTuple `(mu, lstar, tstar)`.
"""
function characteristic_quantities(P1::CelestialObject, P2::CelestialObject)
    m_star = mass(P1) + mass(P2)
    mu_ratio = mass(P2) / m_star
    lstar = orbit_semi_major_axis(P2)
    tstar = sqrt(lstar^3 / (G_KM3 * m_star))
    return (mu=mu_ratio, lstar=lstar, tstar=tstar)
end

# =============================================================================
# Equations of Motion
# =============================================================================

"""
    cr3bp_ode!(du, u, p, t)

CR3BP equations of motion (in-place). State `u = [x, y, z, vx, vy, vz]`, parameter `p = μ`.
"""
function cr3bp_ode!(du, u, p, t)
    μ = p
    x, y, z, vx, vy, vz = u

    d = sqrt((x + μ)^2 + y^2 + z^2)
    r = sqrt((x + μ - 1)^2 + y^2 + z^2)
    d3 = d^3
    r3 = r^3

    du[1] = vx
    du[2] = vy
    du[3] = vz
    du[4] = -(1 - μ) * (x + μ) / d3 - μ * (x - 1 + μ) / r3 + 2vy + x
    du[5] = -(1 - μ) * y / d3 - μ * y / r3 - 2vx + y
    du[6] = -(1 - μ) * z / d3 - μ * z / r3
    return nothing
end

"""
    cr3bp_ode_stm!(du, u, p, t)

CR3BP equations of motion with State Transition Matrix propagation (in-place).
State `u` is 42-element: `[x, y, z, vx, vy, vz, Φ₁₁, Φ₂₁, ..., Φ₆₆]`.
The STM is stored in column-major order. Parameter `p = μ`.
"""
function cr3bp_ode_stm!(du, u, p, t)
    μ = p
    x, y, z, vx, vy, vz = u[1], u[2], u[3], u[4], u[5], u[6]

    # STM as 6×6 matrix (column-major view)
    phi = reshape(@view(u[7:42]), 6, 6)

    d = sqrt((x + μ)^2 + y^2 + z^2)
    r = sqrt((x + μ - 1)^2 + y^2 + z^2)
    d3, d5 = d^3, d^5
    r3, r5 = r^3, r^5

    sigXX = 1 - (1-μ)/d3 - μ/r3 + 3(1-μ)*(x+μ)^2/d5 + 3μ*(x-1+μ)^2/r5
    sigYY = 1 - (1-μ)/d3 - μ/r3 + 3(1-μ)*y^2/d5 + 3μ*y^2/r5
    sigZZ = -(1-μ)/d3 - μ/r3 + 3(1-μ)*z^2/d5 + 3μ*z^2/r5
    sigXY = 3(1-μ)*(x+μ)*y/d5 + 3μ*(x-1+μ)*y/r5
    sigXZ = 3(1-μ)*(x+μ)*z/d5 + 3μ*(x-1+μ)*z/r5
    sigYZ = 3(1-μ)*y*z/d5 + 3μ*y*z/r5

    A = @SMatrix [0.0   0.0   0.0   1.0  0.0  0.0;
                  0.0   0.0   0.0   0.0  1.0  0.0;
                  0.0   0.0   0.0   0.0  0.0  1.0;
                  sigXX sigXY sigXZ 0.0  2.0  0.0;
                  sigXY sigYY sigYZ -2.0 0.0  0.0;
                  sigXZ sigYZ sigZZ 0.0  0.0  0.0]

    # State derivatives
    ax = -(1 - μ) * (x + μ) / d3 - μ * (x - 1 + μ) / r3 + 2vy + x
    ay = -(1 - μ) * y / d3 - μ * y / r3 - 2vx + y
    az = -(1 - μ) * z / d3 - μ * z / r3

    du[1] = vx
    du[2] = vy
    du[3] = vz
    du[4] = ax
    du[5] = ay
    du[6] = az

    # STM derivative: Φ̇ = A * Φ
    phi_dot = A * phi
    du_phi = reshape(@view(du[7:42]), 6, 6)
    du_phi .= phi_dot
    return nothing
end

"""
    two_body_ode!(du, u, p, t)

Two-body problem equations of motion (in-place). State `u = [x,y,z,vx,vy,vz]`, `p = μ`.
"""
function two_body_ode!(du, u, p, t)
    μ = p
    r = @view u[1:3]
    r3 = norm(r)^3
    du[1] = u[4]
    du[2] = u[5]
    du[3] = u[6]
    du[4] = -μ * u[1] / r3
    du[5] = -μ * u[2] / r3
    du[6] = -μ * u[3] / r3
    return nothing
end

# =============================================================================
# Pseudo-Potential and Jacobi Constant
# =============================================================================

"""
    u_star(ndx, μ)

CR3BP pseudo-potential U* at position `ndx[1:3]`.
"""
function u_star(ndx, μ)
    x, y, z = ndx[1], ndx[2], ndx[3]
    d = sqrt((x + μ)^2 + y^2 + z^2)
    r = sqrt((x + μ - 1)^2 + y^2 + z^2)
    return (1 - μ) / d + μ / r + (x^2 + y^2) / 2
end

"""
    jacobi_constant(ndx, μ)

Jacobi constant C for state `ndx = [x, y, z, vx, vy, vz]` and mass ratio `μ`.
"""
function jacobi_constant(ndx, μ)
    ustar = u_star(ndx, μ)
    v = @view ndx[4:6]
    return 2ustar - dot(v, v)
end

# =============================================================================
# Lagrange Points
# =============================================================================

"""
    lagrange_points(μ)

Compute the 5 Lagrange points using hand-rolled Newton-Raphson iteration.
Returns a 5×2 matrix of `[x, y]` coordinates.
"""
function lagrange_points(μ)
    tol = 1e-12

    function colin(x)
        d = abs(x + μ)
        r = abs(x - 1 + μ)
        return x - (1 - μ) * (x + μ) / d^3 - μ * (x - 1 + μ) / r^3
    end

    function fd_colin(x)
        xm = x - tol
        xp = x + tol
        ym = colin(xm)
        y  = colin(x)
        yp = colin(xp)
        dm = (y - ym) / (x - xm)
        dp = (yp - y) / (xp - x)
        return (dm + dp) / 2
    end

    # L1
    f1(gam)  = colin(1 - μ - gam)
    f1p(gam) = -fd_colin(1 - μ - gam)
    gam1 = μ
    delta = 100.0
    while abs(delta) > tol
        delta = f1(gam1) / f1p(gam1)
        gam1 -= delta
    end
    x1 = 1 - μ - gam1

    # L2
    f2(gam)  = colin(1 - μ + gam)
    f2p(gam) = fd_colin(1 - μ + gam)
    gam2 = μ
    delta = 100.0
    while abs(delta) > tol
        delta = f2(gam2) / f2p(gam2)
        gam2 -= delta
    end
    x2 = 1 - μ + gam2

    # L3
    f3(gam)  = colin(-μ - gam)
    f3p(gam) = -fd_colin(-μ - gam)
    gam3 = 1.0
    delta = 100.0
    while abs(delta) > tol
        delta = f3(gam3) / f3p(gam3)
        gam3 -= delta
    end
    x3 = -μ - gam3

    # L4/L5
    x45 = 0.5 - μ
    y4  = sqrt(3) / 2
    y5  = -sqrt(3) / 2

    return [x1 0.0; x2 0.0; x3 0.0; x45 y4; x45 y5]
end

"""
    lagrange_points_roots(μ)

Compute the 5 Lagrange points using `Roots.find_zero`. Returns a 5×2 matrix.
"""
function lagrange_points_roots(μ)
    function collinear(x)
        d = abs(x + μ)
        r = abs(x - 1 + μ)
        return x - (1 - μ) * (x + μ) / d^3 - μ * (x - 1 + μ) / r^3
    end

    x1 = find_zero(collinear, 0.0)
    x2 = find_zero(collinear, 1.0)
    x3 = find_zero(collinear, -1.0)

    x45 = 0.5 - μ
    y4  = sqrt(3) / 2
    y5  = -sqrt(3) / 2

    return [x1 0.0; x2 0.0; x3 0.0; x45 y4; x45 y5]
end

# =============================================================================
# Coordinate Transformations
# =============================================================================

"""
    non_dim(x, lstar, tstar, μ; shift=true)

Non-dimensionalize a dimensional state `[x,y,z,vx,vy,vz]` (km, km/s) using
CR3BP characteristic quantities. If `shift=true`, shifts origin to barycenter.
"""
function non_dim(x, lstar, tstar, μ; shift=true)
    ndx = copy(x)
    ndx[1:3] ./= lstar
    ndx[4:6] .*= tstar / lstar
    if shift
        ndx[1] -= μ
    end
    return ndx
end

"""
    dimensionalize(ndx, lstar, tstar, μ; shift=true)

Convert a non-dimensional CR3BP state to dimensional (km, km/s).
If `shift=true`, shifts origin back to primary body.
"""
function dimensionalize(ndx, lstar, tstar, μ; shift=true)
    x = copy(ndx)
    if shift
        x[1] += μ
    end
    x[1:3] .*= lstar
    x[4:6] .*= lstar / tstar
    return x
end

"""
    cr3bp_to_j2k(rp2, vp2; direction=0)

Rotation matrix for CR3BP → J2000 frame transformation (6×6).
States must be dimensionalized. Epoch-dependent (requires secondary body
position/velocity from ephemeris).

- `direction=0`: CR3BP → J2000 (returns R)
- `direction=1`: J2000 → CR3BP (returns R⁻¹)

Based on Pavlak's PhD dissertation, section 2.3.2.
"""
function cr3bp_to_j2k(rp2, vp2; direction=0)
    lstar = norm(rp2)
    h_ = cross(rp2, vp2)
    h = h_ / lstar^2
    xhat = rp2 / lstar
    zhat = h / norm(h)
    yhat = cross(zhat, xhat)
    hn = norm(h)

    R0 = hcat(xhat, yhat, zhat)  # 3×3
    R2 = hcat(hn .* yhat, -hn .* xhat, zeros(3))  # 3×3
    R = [R0       zeros(3, 3);
         R2       R0]

    if direction == 0
        return R
    else
        return inv(R)
    end
end

# =============================================================================
# Eigendecomposition
# =============================================================================

"""
    lagrange_point_eigendecomp(ndx, μ)

Eigenvalue decomposition of the CR3BP state-space Jacobian evaluated at position `ndx`.
Returns `(eigenvalues, eigenvectors)`.
"""
function lagrange_point_eigendecomp(ndx, μ)
    x, y, z = ndx[1], ndx[2], ndx[3]

    d = sqrt((x + μ)^2 + y^2 + z^2)
    r = sqrt((x + μ - 1)^2 + y^2 + z^2)
    d3, d5 = d^3, d^5
    r3, r5 = r^3, r^5

    Uxx = 1 - (1-μ)/d3 - μ/r3 + 3(1-μ)*(x+μ)^2/d5 + 3μ*(x-1+μ)^2/r5
    Uyy = 1 - (1-μ)/d3 - μ/r3 + 3(1-μ)*y^2/d5 + 3μ*y^2/r5
    Uzz = -(1-μ)/d3 - μ/r3 + 3(1-μ)*z^2/d5 + 3μ*z^2/r5
    Uxy = 3(1-μ)*(x+μ)*y/d5 + 3μ*(x-1+μ)*y/r5
    Uxz = 3(1-μ)*(x+μ)*z/d5 + 3μ*(x-1+μ)*z/r5
    Uyz = 3(1-μ)*y*z/d5 + 3μ*y*z/r5

    A = [0.0  0.0  0.0  1.0  0.0  0.0;
         0.0  0.0  0.0  0.0  1.0  0.0;
         0.0  0.0  0.0  0.0  0.0  1.0;
         Uxx  Uxy  Uxz  0.0  2.0  0.0;
         Uxy  Uyy  Uyz -2.0  0.0  0.0;
         Uxz  Uyz  Uzz  0.0  0.0  0.0]

    F = eigen(A)
    return F.values, F.vectors
end

# =============================================================================
# I/O
# =============================================================================

"""
    periodic_orbit_df()

Load periodic Lagrange orbit data from CSV. Assumes `periodicLagrangeOrbits.csv`
is in the same directory as this module file.
"""
function periodic_orbit_df()
    data_path = joinpath(@__DIR__, "periodicLagrangeOrbits.csv")
    return CSV.read(data_path, DataFrame)
end

# =============================================================================
# Plotting
# =============================================================================

"""
    plot_zvc(μ, C; fill=false)

Plot the zero-velocity curve for mass ratio `μ` and Jacobi constant `C`.
Uses Plots.jl. If `fill=true`, produces filled contours.
"""
function plot_zvc(μ, C; fill=false)
    npoints = 1000
    xs = range(-1.5, 1.5, length=npoints)
    ys = range(0.0, 1.5, length=npoints)
    eps = 0.1

    ZZ = Matrix{Float64}(undef, npoints, npoints)
    fill!(ZZ, NaN)

    for i in 1:npoints, j in 1:npoints
        ndx = [xs[i], ys[j], 0.0, 0.0, 0.0, 0.0]
        temp = jacobi_constant(ndx, μ)
        if C >= temp && temp > C * (1 - eps)
            ZZ[j, i] = temp
        end
    end

    if fill
        contourf(xs, ys, ZZ)
        contourf!(xs, -ys, ZZ)
    else
        contour(xs, ys, ZZ)
        contour!(xs, -ys, ZZ)
    end
end

# =============================================================================
# Periodic Orbit Solvers
# =============================================================================

"""
    print_corrector_status(iteration, Fnorm, dx; label="")

Print iteration info for differential corrector solvers (solve_lyapunov, solve_halo).
"""
function print_corrector_status(iteration::Int, Fnorm::Float64, dx::AbstractVector; label::String="")
    header = isempty(label) ? "" : "[$label] "
    dx_str = join((@sprintf("% .6e", v) for v in dx), ", ")
    println("  $(header)Iter $(lpad(iteration, 3))  |F| = $(@sprintf("%.6e", Fnorm))  dx = [$dx_str]")
end

# Constraint functions
lyapunov_constraints(ndx) = [ndx[2], ndx[4]]       # y = 0, vx = 0
halo_constraints(ndx)     = [ndx[2], ndx[4], ndx[6]] # y = 0, vx = 0, vz = 0

# Event condition: y-axis crossing
x_crossing_condition(u, t, integrator) = u[2]

function _make_crossing_callback(direction; tol=1e-13, tmin=1e-9)
    affect!(integrator) = terminate!(integrator)

    # Guard against the trivial root at t=0 for states initialized on y=0.
    function guarded_condition(u, t, integrator)
        if t <= tmin
            return one(eltype(u))
        end
        return x_crossing_condition(u, t, integrator)
    end

    if direction > 0
        return ContinuousCallback(guarded_condition, affect!;
                                  affect_neg! = nothing, abstol=tol)
    elseif direction < 0
        return ContinuousCallback(guarded_condition, nothing, affect!;
                                  abstol=tol)
    else
        return ContinuousCallback(guarded_condition, affect!;
                                  abstol=tol)
    end
end

"""
Build the initial 42-element state vector for STM integration:
6 state variables + 36 entries of the identity STM (column-major).
"""
function _stm_initial_state(ndx6::AbstractVector)
    u0 = zeros(42)
    u0[1:6] .= ndx6
    # Identity matrix in column-major flat form
    for i in 1:6
        u0[6 + (i-1)*6 + i] = 1.0
    end
    return u0
end

"""
    solve_lyapunov(ndx, t; mu=0.0121, tol=1e-12, direction=1)

Solve for a periodic Lyapunov orbit using differential correction.

# Arguments
- `ndx`: 6-element initial state guess `[x, y, z, vx, vy, vz]`
- `t`: initial guess for half-period
- `mu`: CR3BP mass ratio (default Earth-Moon)
- `tol`: convergence tolerance
- `direction`: +1 for positive y-crossing, -1 for negative, 0 for either

# Returns
- `sol`: DifferentialEquations solution object (full period)
- `period`: full orbital period
"""
function solve_lyapunov(ndx, t; mu=0.0121, tol=1e-12, direction=1)
    effective_direction = direction
    if abs(ndx[2]) <= sqrt(tol)
        expected_direction = ndx[5] < 0 ? 1 : -1
        effective_direction = direction == 0 ? expected_direction : direction
        if direction != 0 && direction != expected_direction
            @warn "Initial state is on y=0; requested crossing direction may miss first nontrivial crossing" direction expected_direction vy0=ndx[5]
        end
    end

    cb = _make_crossing_callback(effective_direction; tol=tol/10)

    # Build 42-element initial state with identity STM
    x0 = _stm_initial_state(ndx)

    # Initial integration to find x-axis crossing
    prob = ODEProblem(cr3bp_ode_stm!, x0, (0.0, t), mu)
    sol = solve(prob, Vern9(); abstol=tol, reltol=tol, callback=cb)

    t0 = sol.t[end]

    F = lyapunov_constraints(sol.u[end])
    Fnorm = norm(F)
    iteration = 0

    while Fnorm > tol
        uf = sol.u[end]
        x, y = uf[1], uf[2]
        y_dot = uf[5]
        d = sqrt((x + mu)^2 + y^2)
        r = sqrt((x - 1 + mu)^2 + y^2)
        x_ddot = -(1 - mu) * (x + mu) / d^3 - mu * (x - 1 + mu) / r^3 + 2y_dot + x

        phi = reshape(uf[7:42], 6, 6)

        F = lyapunov_constraints(uf)
        Fnorm = norm(F)

        # Jacobian of constraints w.r.t. [vy₀, t]
        # Python phi[1,4] → Julia phi[2,5], Python phi[3,4] → Julia phi[4,5]
        DF = [phi[2,5] y_dot;
              phi[4,5] x_ddot]
        dx = DF \ F

        iteration += 1
        print_corrector_status(iteration, Fnorm, dx; label="Lyapunov")

        x0[5] -= dx[1]
        t0    -= dx[2]

        prob = ODEProblem(cr3bp_ode_stm!, x0, (0.0, t0), mu)
        sol = solve(prob, Vern9(); abstol=tol, reltol=tol, callback=cb)
    end

    # Full period integration
    prob = ODEProblem(cr3bp_ode_stm!, x0, (0.0, 2t0), mu)
    sol = solve(prob, Vern9(); abstol=tol, reltol=tol)
    return sol, sol.t[end]
end

"""
    lyapunov_continuation(ndx, t; mu=0.0121, dx=0.001, lim=0.9, direction=1)

Generate a family of Lyapunov orbits by natural parameter continuation in x.

Returns a vector of DifferentialEquations solution objects.
"""
function lyapunov_continuation(ndx, t; mu=0.0121, dx=0.001, lim=0.9, direction=1)
    traj_list = []
    x0 = copy(ndx)
    t0 = copy(t)

    # First orbit
    sol1, period = solve_lyapunov(x0, t0; mu=mu, tol=1e-12, direction=direction)
    push!(traj_list, sol1)

    sign_val = sign(lim - sol1.u[1][1])
    sign_last = sign_val

    while sign_val == sign_last
        # Initial guess from last converged orbit
        ig = copy(traj_list[end].u[1][1:6])
        t_ig = period / 2

        ig[1] += dx
        try
            sol, period = solve_lyapunov(ig, t_ig; mu=mu, tol=1e-12, direction=direction)
            push!(traj_list, sol)
            sign_last = sign_val
            sign_val = sign(lim - sol.u[1][1])
        catch
            continue
        end
    end

    return traj_list
end

"""
    solve_halo(ndx, t; mu=0.0121, direction=-1)

Solve for a periodic Halo orbit using differential correction.

# Arguments
- `ndx`: 6-element initial state guess `[x, y, z, vx, vy, vz]`
- `t`: initial guess for half-period (or integration time to first crossing)
- `mu`: CR3BP mass ratio
- `direction`: crossing direction for event detection (+1/-1), or 0 for either

# Returns
- `sol`: DifferentialEquations solution object (full period)
- `period`: full orbital period
"""
function solve_halo(ndx, t; mu=0.0121, direction=-1)
    tol = 1e-12

    effective_direction = direction
    if abs(ndx[2]) <= sqrt(tol)
        # For y0≈0, first nontrivial crossing is opposite the sign of vy0.
        expected_direction = ndx[5] < 0 ? 1 : -1
        effective_direction = direction == 0 ? expected_direction : direction
        if direction != 0 && direction != expected_direction
            @warn "Initial state is on y=0; requested crossing direction may miss first nontrivial crossing" direction expected_direction vy0=ndx[5]
        end
    end

    cb = _make_crossing_callback(effective_direction; tol=tol/10)

    x0 = _stm_initial_state(ndx)

    # Initial integration to x-axis crossing
    prob = ODEProblem(cr3bp_ode_stm!, x0, (0.0, t), mu)
    sol = solve(prob, Vern9(); abstol=tol, reltol=tol, callback=cb)

    t0 = sol.t[end]

    F = halo_constraints(sol.u[end])
    Fnorm = norm(F)
    iteration = 0

    while Fnorm > tol
        uf = sol.u[end]
        x, y, z = uf[1], uf[2], uf[3]
        y_dot = uf[5]
        d = sqrt((x + mu)^2 + y^2 + z^2)
        r = sqrt((x - 1 + mu)^2 + y^2 + z^2)
        x_ddot = -(1 - mu) * (x + mu) / d^3 - mu * (x - 1 + mu) / r^3 + 2y_dot + x
        z_ddot = -(1 - mu) * z / d^3 - mu * z / r^3

        phi = reshape(uf[7:42], 6, 6)

        F = halo_constraints(uf)
        Fnorm = norm(F)

        # Jacobian of constraints w.r.t. [x₀, vy₀, t]
        # Python phi[1,0]→phi[2,1], phi[3,0]→phi[4,1], phi[5,0]→phi[6,1]
        # Python phi[1,4]→phi[2,5], phi[3,4]→phi[4,5], phi[5,4]→phi[6,5]
        DF = [phi[2,1]  phi[2,5]  y_dot;
              phi[4,1]  phi[4,5]  x_ddot;
              phi[6,1]  phi[6,5]  z_ddot]
        dx = DF \ F

        iteration += 1
        print_corrector_status(iteration, Fnorm, dx; label="Halo")

        x0[1] -= dx[1]
        x0[5] -= dx[2]
        t0    -= dx[3]

        prob = ODEProblem(cr3bp_ode_stm!, x0, (0.0, t0), mu)
        sol = solve(prob, Vern9(); abstol=tol, reltol=tol, callback=cb)
    end

    # Full period integration
    prob = ODEProblem(cr3bp_ode_stm!, x0, (0.0, 2t0), mu)
    sol = solve(prob, Vern9(); abstol=tol, reltol=tol)
    return sol, sol.t[end]
end

"""
    halo_continuation(ndx, t; mu=0.0121, dx=-0.001, lim=0.04, direction=-1)

Generate a family of Halo orbits by natural parameter continuation in z.

Returns a vector of DifferentialEquations solution objects.
"""
function halo_continuation(ndx, t; mu=0.0121, dx=-0.001, lim=0.04, direction=-1)
    traj_list = []

    sol1, period = solve_halo(ndx, t; mu=mu, direction=direction)
    push!(traj_list, sol1)

    sign_val = sign(lim - sol1.u[1][3])
    sign_last = sign_val

    while sign_val == sign_last
        ig = copy(traj_list[end].u[1][1:6])
        t_ig = period

        ig[3] += dx  # adjust along z direction
        try
            sol, period = solve_halo(ig, t_ig; mu=mu, direction=direction)
            push!(traj_list, sol)
            sign_last = sign_val
            sign_val = sign(lim - sol.u[1][3])
        catch
            continue
        end
    end

    return traj_list
end

end # module Astro

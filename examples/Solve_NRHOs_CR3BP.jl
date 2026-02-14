"""
Propagate_NRHO_Astro.jl

Minimal example showing how to use `Astro.jl` to:
1) Build Earth-Moon CR3BP parameters,
2) Differentially correct an NRHO-like initial guess with `solve_halo`,
3) Plot and save the corrected periodic trajectory.

Initial conditions come from `examples/Solve_NRHOs_CR3BP.py`
(9:2 and 4:1 L2 Southern NRHO initial guesses).
"""

if !isdefined(Main, :Astro)
    include(joinpath(@__DIR__, "..", "Astro.jl"))
end
import .Astro
using Plots

# Earth-Moon CR3BP characteristic quantities
earth = Astro.CelestialObject("Earth")
moon = Astro.CelestialObject("Moon")
cq = Astro.characteristic_quantities(earth, moon)

μ = cq.mu
tstar = cq.tstar
day_nd = 86400.0 / tstar

# 9:2 L2 Southern NRHO initial guess (non-dimensional)
x0_92 = [
    1.02134,
    0.0,
    -0.18162,
    0.0,
    -0.10176,
    9.76561e-07,
]

# 4:1 L2 Southern NRHO initial guess (non-dimensional)
x0_41 = [
    1.03545,
    0.0,
    -0.19003,
    0.0,
    -0.13071,
    5.62991e-07,
]

# Half-period guesses from the Python example
tf_half_92 = 1.50206
tf_half_41 = 1.68981

# Differentially correct and propagate one full corrected period
sol_corr_92, period_nd_92 = Astro.solve_halo(x0_92, tf_half_92; mu=μ, direction=0)
sol_corr_41, period_nd_41 = Astro.solve_halo(x0_41, tf_half_41; mu=μ, direction=0)

# Diagnostics
C0_92 = Astro.jacobi_constant(sol_corr_92.u[1], μ)
Cf_92 = Astro.jacobi_constant(sol_corr_92.u[end], μ)
period_days_92 = period_nd_92 / day_nd

C0_41 = Astro.jacobi_constant(sol_corr_41.u[1], μ)
Cf_41 = Astro.jacobi_constant(sol_corr_41.u[end], μ)
period_days_41 = period_nd_41 / day_nd

println("=== 9:2 L2 Southern NRHO (corrected) ===")
println("Corrected initial state [x y z vx vy vz]:")
println(sol_corr_92.u[1])
println("Corrected period [ND]:    ", period_nd_92)
println("Initial Jacobi constant: ", C0_92)
println("Final Jacobi constant:   ", Cf_92)
println("|ΔC|:                    ", abs(Cf_92 - C0_92))
println("Corrected period [days]: ", period_days_92)

println()
println("=== 4:1 L2 Southern NRHO (corrected) ===")
println("Corrected initial state [x y z vx vy vz]:")
println(sol_corr_41.u[1])
println("Corrected period [ND]:    ", period_nd_41)
println("Initial Jacobi constant: ", C0_41)
println("Final Jacobi constant:   ", Cf_41)
println("|ΔC|:                    ", abs(Cf_41 - C0_41))
println("Corrected period [days]: ", period_days_41)

# Plot corrected 9:2 trajectory in rotating frame
xs_92 = [u[1] for u in sol_corr_92.u]
ys_92 = [u[2] for u in sol_corr_92.u]
zs_92 = [u[3] for u in sol_corr_92.u]

plt_92 = plot(
    xs_92,
    ys_92,
    zs_92,
    label="NRHO corrected periodic orbit",
    linewidth=2,
    xlabel="x [ND]",
    ylabel="y [ND]",
    zlabel="z [ND]",
    aspect_ratio=:equal,
    title="Earth-Moon CR3BP Corrected 9:2 NRHO (Astro.jl)",
    legend=:topright,
)

# Primary/secondary locations in CR3BP rotating frame
scatter!(plt_92, [-μ], [0.0], [0.0], label="Earth", markersize=5)
scatter!(plt_92, [1 - μ], [0.0], [0.0], label="Moon", markersize=4)

png_path_92 = joinpath(@__DIR__, "output", "nrho_92_corrected_trajectory.png")
savefig(plt_92, png_path_92)
println("Saved plot: ", png_path_92)

# Plot corrected 4:1 trajectory in rotating frame
xs_41 = [u[1] for u in sol_corr_41.u]
ys_41 = [u[2] for u in sol_corr_41.u]
zs_41 = [u[3] for u in sol_corr_41.u]

plt_41 = plot(
    xs_41,
    ys_41,
    zs_41,
    label="NRHO corrected periodic orbit",
    linewidth=2,
    xlabel="x [ND]",
    ylabel="y [ND]",
    zlabel="z [ND]",
    aspect_ratio=:equal,
    title="Earth-Moon CR3BP Corrected 4:1 NRHO (Astro.jl)",
    legend=:topright,
)

scatter!(plt_41, [-μ], [0.0], [0.0], label="Earth", markersize=5)
scatter!(plt_41, [1 - μ], [0.0], [0.0], label="Moon", markersize=4)

png_path_41 = joinpath(@__DIR__, "output", "nrho_41_corrected_trajectory.png")
savefig(plt_41, png_path_41)
println("Saved plot: ", png_path_41)

# Plot both corrected trajectories together in one figure
plt_both = plot(
    xs_92,
    ys_92,
    zs_92,
    label="9:2 corrected NRHO",
    linewidth=2,
    xlabel="x [ND]",
    ylabel="y [ND]",
    zlabel="z [ND]",
    aspect_ratio=:equal,
    title="Earth-Moon CR3BP Corrected NRHOs (9:2 and 4:1)",
    legend=:topright,
)

plot!(plt_both, xs_41, ys_41, zs_41, label="4:1 corrected NRHO", linewidth=2)
scatter!(plt_both, [-μ], [0.0], [0.0], label="Earth", markersize=5)
scatter!(plt_both, [1 - μ], [0.0], [0.0], label="Moon", markersize=4)

png_path_both = joinpath(@__DIR__, "output", "nrho_92_41_corrected_trajectory.png")
savefig(plt_both, png_path_both)
println("Saved plot: ", png_path_both)

display(plt_92)
display(plt_41)
display(plt_both)

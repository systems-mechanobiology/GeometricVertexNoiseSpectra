#!/usr/bin/env julia

"""
Figure 6: Sensitivity of Witness Motifs to Angular Perturbations

Generates a 3×2 panel figure showing, for each motif degree d=1..6, the
maximum weighted-degree change ΔWD as a function of hub deviation angle δ₀
from the bisector of the unoccupied arc.

Two curves per panel:
  - Deterministic baseline (dashed): hub moved at angle bisector+δ₀, spokes
    pushed radially with zero angular noise.
  - Envelope (solid): maximum ΔWD over N=50 random spoke-angle perturbations,
    using common random numbers (fixed seed).

For d ≤ 5 the extremum is unique at δ₀=0; for d=6, hexagonal symmetry
induces multiple equivalent maxima.

Parameters: k=100, r=50, δmax_deg=60, ngrid=121, N=50, seed=1

Author: Ben Cardoen
Source: scripts/patterns_new.jl (functions extracted verbatim)
"""

import Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using DynamicGeometricGraphs
using CairoMakie
using LinearAlgebra
using StaticArrays
using Random

const DDG = DynamicGeometricGraphs

# ============================================================================
# Core geometry helpers (from scripts/patterns_new.jl)
# ============================================================================

@inline function rot2(v::SVector{2,Float64}, ε::Float64)
    c, s = cos(ε), sin(ε)
    return SVector{2,Float64}(c*v[1] - s*v[2], s*v[1] + c*v[2])
end

function bisector_angle(d::Int)
    last = -(d - 1) * (π/3)
    occ  = abs(last)
    unocc = 2π - occ
    return last - unocc/2
end

function config_hub_and_spoke_noise(g::DynamicGeometricGraph{2,Float64}, d::Int, r::Float64, θ::Float64,
                                    δspokes::Vector{Float64}; hub_vid::Int=1)
    g_new = DDG.freeze(g)

    h  = DDG.get_vertex_coords(g, hub_vid)
    h′ = h + SVector{2,Float64}(r*cos(θ), r*sin(θ))
    DDG.update_coord!(g_new, h, h′)

    spoke_idx = 0
    for vid in sort(collect(DDG.vertices(g)))
        vid == hub_vid && continue
        spoke_idx += 1

        k  = DDG.get_vertex_coords(g, vid)
        v  = k - h′
        nv = norm(v)
        nv == 0 && continue

        dir = v / nv
        dir = rot2(dir, δspokes[spoke_idx])
        k′  = k + r * dir
        DDG.update_coord!(g_new, k, k′)
    end

    return g_new
end

# ============================================================================
# Envelope computation (from scripts/patterns_new.jl — final version)
# ============================================================================

function envelope_optimality(g, d; r::Float64, hub_vid::Int=1, δmax_deg::Real=30,
                             ngrid::Int=121, N::Int=50, seed::Int=1)

    rng = MersenneTwister(seed)
    bis = bisector_angle(d)
    wd0 = DDG.weighted_degree(g, hub_vid)

    δmax = deg2rad(δmax_deg)
    spoke_samples = [(rand(rng, d) .* 2 .- 1) .* δmax for _ in 1:N]

    δ0_grid = range(-deg2rad(δmax_deg), deg2rad(δmax_deg), length=ngrid)
    best = fill(-Inf, ngrid)
    det  = fill(-Inf, ngrid)

    δzero = zeros(Float64, d)

    for (j, δ0) in enumerate(δ0_grid)
        θ = bis + δ0

        gp0 = config_hub_and_spoke_noise(g, d, r, θ, δzero; hub_vid=hub_vid)
        Δ0  = DDG.weighted_degree(gp0, hub_vid) - wd0
        det[j] = Δ0

        best_j = Δ0
        for δsp in spoke_samples
            gp  = config_hub_and_spoke_noise(g, d, r, θ, collect(δsp); hub_vid=hub_vid)
            Δ   = DDG.weighted_degree(gp, hub_vid) - wd0
            best_j = max(best_j, Δ)
        end
        best[j] = best_j
    end

    return rad2deg.(collect(δ0_grid)), det, best
end

# ============================================================================
# Main: generate 3×2 panel figure
# ============================================================================

function generate_figure(; k=100.0, r=k/2, δmax_deg=60, ngrid=121, N=50, seed=1)
    fig = Figure(size=(1400, 900), fontsize=20)

    for d in 1:6
        g = DDG.generate_hub_spoke_graph(k, d)
        xdeg, det, env = envelope_optimality(g, d; r=r, δmax_deg=δmax_deg, ngrid=ngrid, N=N, seed=seed)

        row = (d - 1) % 3 + 1   # 1..3
        col = (d - 1) ÷ 3 + 1   # 1..2

        ax = Axis(fig[row, col]; title="d=$d")
        lines!(ax, xdeg, env; linewidth=3)
        lines!(ax, xdeg, det; linewidth=2, linestyle=:dash)
        vlines!(ax, [0.0]; linestyle=:dash)
        ax.xlabel = "hub deviation δ₀ (deg)"
        ax.ylabel = "max ΔWD"

        # Print peak values for verification
        mid = (ngrid + 1) ÷ 2
        println("  d=$d: ΔWD(δ₀=0) det=$(round(det[mid], digits=2)), env=$(round(env[mid], digits=2)), " *
                "det_max=$(round(maximum(det), digits=2)), env_max=$(round(maximum(env), digits=2))")
    end

    return fig
end

# ============================================================================
# Run
# ============================================================================

println("=" ^ 70)
println("Generating Figure 6: Motif Angular Sensitivity (Optimality)")
println("=" ^ 70)

fig = generate_figure(; k=100.0, δmax_deg=60, ngrid=121, N=50, seed=1)

output_pdf = joinpath(@__DIR__, "Figure6.pdf")
output_png = joinpath(@__DIR__, "Figure6.png")
save(output_pdf, fig)
save(output_png, fig, px_per_unit=2)

println("\nSaved: $output_pdf")
println("Saved: $output_png")
println("=" ^ 70)

#!/usr/bin/env julia

"""
Figure 5: Motif-informed graph-level sensitivity and spanner effects (S3I)

2×2 layout:
  - Panel A (1,1): S3I c₁=1 vs c₂=[2,3,4,5]  — noise amplitude separation
  - Panel B (1,2): S3I L₁=10 vs L₂=[20,30,40,50] — noise radius separation
  - Panel C (2,1): Motif frequency differences (1,2,3,4,5,6 vs 4×6)
  - Panel D (2,2): ε-thick vs unconstrained graph

Each panel: jittered scatter + bootstrap CI lines + significance checkmarks.
Saves CSV tables for all raw S3I values.

Author: Ben Cardoen
Source: archive/fig10_s3i_tests/generate.jl
"""

import Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using CalibratedTemperedPowerLaw
using DynamicGeometricGraphs
using CairoMakie
using Statistics
using Random
using Graphs
using LinearAlgebra
using StaticArrays
using Printf

const DDG = DynamicGeometricGraphs
const SEED = 42

# ============================================================================
# Graph construction helpers
# ============================================================================

function split_refgraph(g)
    gp = freeze(g)
    rem_vertex!(gp, 1)
    add_edge!(gp, 24, 9)
    add_edge!(gp, 13, 3)
    add_edge!(gp, 6, 21)
    gpp = freeze(gp)
    h = DDG.add_vertex!(gpp, SVector{2,Float64}(0.0, 0.0))
    add_edge!(gpp, h, 9)
    add_edge!(gpp, h, 3)
    add_edge!(gpp, h, 21)
    add_edge!(gpp, h, 24)
    add_edge!(gpp, h, 13)
    add_edge!(gpp, h, 6)
    @assert nv(gp) != nv(gpp)
    return gp, gpp
end

# ============================================================================
# Parameters
# ============================================================================

n_samples = 200
n_runs = 100
c1 = 1
limit = 40
p = 0.05

# ============================================================================
# Build graphs
# ============================================================================

G = refgraph()
gp, gpp = split_refgraph(G)
nu = refgraph(; pattern=[4,4,4,4,4,4])

l1 = autotune_lambda(c1, limit, p; mu=0.0)

println("=" ^ 70)
println("Generating Figure 5: S3I graph-level sensitivity")
println("=" ^ 70)
println("  n_samples=$n_samples, n_runs=$n_runs, c1=$c1, limit=$limit, p=$p")
println("  λ₁ ≈ $(@sprintf("%.6f", l1))")
println("  G: nv=$(nv(G)), ne=$(ne(G))")
println("  Gp (unconstrained): nv=$(nv(gp)), ne=$(ne(gp))")
println("  Gpp (ε-thick): nv=$(nv(gpp)), ne=$(ne(gpp))")
println("  ν (uniform 4×6): nv=$(nv(nu)), ne=$(ne(nu))")

# ============================================================================
# Panel A: S3I c₁=1 vs c₂=[2,3,4,5]
# ============================================================================

println("\n--- Panel A: varying c ---")
Random.seed!(SEED)

c2_levels = [2, 3, 4, 5]
s3i_c = Matrix{Float64}(undef, n_runs, length(c2_levels))

for (col_idx, c2_val) in enumerate(c2_levels)
    l2_temp = autotune_lambda(c2_val, limit, p; mu=0.0)
    println("  c₂=$c2_val, λ₂≈$(@sprintf("%.6f", l2_temp))")
    for r in 1:n_runs
        Da = strong_distances(G; n_samples=n_samples, λ=l1, c=c1, limit=limit, seed=r)
        Db = strong_distances(G; n_samples=n_samples, λ=l2_temp, c=c2_val, limit=limit, seed=r)
        SCa = SC(Da; ks=[1,2])
        SCb = SC(Db; ks=[1,2])
        s3i_c[r, col_idx] = S3I_from_SC(SCa[1], SCa[2], SCb[1], SCb[2])
    end
end

# ============================================================================
# Panel B: S3I L₁=10 vs L₂=[20,30,40,50]
# ============================================================================

println("\n--- Panel B: varying L ---")

limit_levels = [20, 30, 40, 50]
limit_base = 10
s3i_L = Matrix{Float64}(undef, n_runs, length(limit_levels))

for (col_idx, limit_val) in enumerate(limit_levels)
    l_base = autotune_lambda(c1, limit_base, p; mu=0.0)
    l_temp = autotune_lambda(c1, limit_val, p; mu=0.0)
    println("  L₂=$limit_val, λ_base≈$(@sprintf("%.6f", l_base)), λ₂≈$(@sprintf("%.6f", l_temp))")
    for r in 1:n_runs
        Da = strong_distances(G; n_samples=n_samples, λ=l_base, c=c1, limit=limit_base, seed=r)
        Db = strong_distances(G; n_samples=n_samples, λ=l_temp, c=c1, limit=limit_val, seed=r)
        SCa = SC(Da; ks=[1,2])
        SCb = SC(Db; ks=[1,2])
        s3i_L[r, col_idx] = S3I_from_SC(SCa[1], SCa[2], SCb[1], SCb[2])
    end
end

# ============================================================================
# Panel C: Motif frequency (G vs ν) at c=[1,2,3,4,5]
# ============================================================================

println("\n--- Panel C: motif frequency ---")

c_motif_levels = [1, 2, 3, 4, 5]
limit_motif = 40
s3i_motif = Matrix{Float64}(undef, n_runs, length(c_motif_levels))

for (col_idx, c_val) in enumerate(c_motif_levels)
    l_temp = autotune_lambda(c_val, limit_motif, p; mu=0.0)
    println("  c=$c_val, λ≈$(@sprintf("%.6f", l_temp))")
    for r in 1:n_runs
        D_G = strong_distances(G; n_samples=n_samples, λ=l_temp, c=c_val, limit=limit_motif, seed=r)
        SC_G = SC(D_G; ks=[1,2])
        D_ν = strong_distances(nu; n_samples=n_samples, λ=l_temp, c=c_val, limit=limit_motif, seed=r)
        SC_ν = SC(D_ν; ks=[1,2])
        s3i_motif[r, col_idx] = S3I_from_SC(SC_G[1], SC_G[2], SC_ν[1], SC_ν[2])
    end
end

# ============================================================================
# Panel D: ε-thick (Gpp) vs unconstrained (Gp) at c=[1,2,3,4,5]
# ============================================================================

println("\n--- Panel D: ε-thick vs unconstrained ---")

c_span_levels = [1, 2, 3, 4, 5]
limit_span = 40
s3i_span = Matrix{Float64}(undef, n_runs, length(c_span_levels))

for (col_idx, c_val) in enumerate(c_span_levels)
    l_temp = autotune_lambda(c_val, limit_span, p; mu=0.0)
    println("  c=$c_val, λ≈$(@sprintf("%.6f", l_temp))")
    for r in 1:n_runs
        D_p = strong_distances(gp; n_samples=n_samples, λ=l_temp, c=c_val, limit=limit_span, seed=r)
        SC_p = SC(D_p; ks=[1,2])
        D_pp = strong_distances(gpp; n_samples=n_samples, λ=l_temp, c=c_val, limit=limit_span, seed=r)
        SC_pp = SC(D_pp; ks=[1,2])
        s3i_span[r, col_idx] = S3I_from_SC(SC_p[1], SC_p[2], SC_pp[1], SC_pp[2])
    end
end

# ============================================================================
# Save CSV tables
# ============================================================================

println("\n--- Saving CSV data ---")

function save_s3i_csv(path, data, col_labels, col_header)
    open(path, "w") do io
        println(io, "run,$(col_header),s3i")
        for col_idx in 1:size(data, 2)
            for r in 1:size(data, 1)
                @printf(io, "%d,%s,%.8f\n", r, col_labels[col_idx], data[r, col_idx])
            end
        end
    end
    println("  Saved: $path")
end

save_s3i_csv(joinpath(@__DIR__, "s3i_panel_a_vary_c.csv"),
             s3i_c, ["c2=$c" for c in c2_levels], "comparison")

save_s3i_csv(joinpath(@__DIR__, "s3i_panel_b_vary_L.csv"),
             s3i_L, ["L2=$l" for l in limit_levels], "comparison")

save_s3i_csv(joinpath(@__DIR__, "s3i_panel_c_motif_freq.csv"),
             s3i_motif, ["c=$c" for c in c_motif_levels], "noise_level")

save_s3i_csv(joinpath(@__DIR__, "s3i_panel_d_spanner.csv"),
             s3i_span, ["c=$c" for c in c_span_levels], "noise_level")

# ============================================================================
# Figure
# ============================================================================

println("\n--- Plotting ---")

f = Figure(size=(1400, 1200), fontsize=18)

# Okabe-Ito scatter colors per group
const OI_COLORS = [colorant"#0072B2", colorant"#E69F00", colorant"#009E73",
                   colorant"#D55E00", colorant"#CC79A7"]

function plot_s3i_panel!(ax, data, labels; ylims_range=nothing)
    hlines!(ax, [0.0]; color=:gray, linestyle=:dash, linewidth=1)
    results = []
    for col_idx in 1:size(data, 2)
        vals = data[:, col_idx]
        result = bootstrap_s3i_test(vals; α=0.05)
        push!(results, (col_idx, result))

        # Jittered scatter
        x = col_idx .+ 0.15 .* (rand(size(data, 1)) .- 0.5)
        scatter!(ax, x, vals; markersize=4, alpha=0.5,
                 color=(OI_COLORS[mod1(col_idx, length(OI_COLORS))], 0.6))

        # CI lines
        linesegments!(ax, [col_idx - 0.3, col_idx + 0.3],
                       [result.ci_lower, result.ci_lower]; color=:black, linewidth=2)
        linesegments!(ax, [col_idx - 0.3, col_idx + 0.3],
                       [result.ci_upper, result.ci_upper]; color=:black, linewidth=2)
    end
    ax.xticks = (1:length(labels), labels)

    if ylims_range !== nothing
        ylims!(ax, ylims_range...)
    end

end

# Panel A
ax_a = Axis(f[1,1], title="S3I: c₁=1 vs c₂=[2,3,4,5]",
            xlabel="Comparison", ylabel="S3I")
Random.seed!(SEED + 1000)
plot_s3i_panel!(ax_a, s3i_c,
    ["c₁ vs c₂=$c" for c in c2_levels]; ylims_range=(0.0, 1.8))

# Panel B
ax_b = Axis(f[1,2], title="S3I: L₁=10 vs L₂=[20,30,40,50]",
            xlabel="Comparison", ylabel="S3I")
plot_s3i_panel!(ax_b, s3i_L,
    ["L₁ vs L₂=$l" for l in limit_levels]; ylims_range=(1.0, 3.0))

# Panel C
ax_c = Axis(f[2,1], title="Motif frequency differences (1,2,3,4,5,6 vs 4×6)",
            xlabel="Noise level", ylabel="S3I")
plot_s3i_panel!(ax_c, s3i_motif,
    ["c=$c" for c in c_motif_levels]; ylims_range=(0.0, 0.6))

# Panel D
ax_d = Axis(f[2,2], title="ε-thick vs unconstrained graph",
            xlabel="Noise level", ylabel="S3I")
plot_s3i_panel!(ax_d, s3i_span,
    ["c=$c" for c in c_span_levels]; ylims_range=(0.0, 0.6))

# ============================================================================
# Save figure
# ============================================================================

output_pdf = joinpath(@__DIR__, "Figure5.pdf")
output_png = joinpath(@__DIR__, "Figure5.png")
save(output_pdf, f)
save(output_png, f, px_per_unit=2)

println("\nSaved: $output_pdf")
println("Saved: $output_png")

# ============================================================================
# Summary statistics
# ============================================================================

println("\n" * "=" ^ 70)
println("Summary statistics")
println("=" ^ 70)

for (label, data, cols) in [
    ("Panel A (vary c)", s3i_c, ["c₂=$c" for c in c2_levels]),
    ("Panel B (vary L)", s3i_L, ["L₂=$l" for l in limit_levels]),
    ("Panel C (motif freq)", s3i_motif, ["c=$c" for c in c_motif_levels]),
    ("Panel D (spanner)", s3i_span, ["c=$c" for c in c_span_levels]),
]
    println("\n$label:")
    for (i, col) in enumerate(cols)
        vals = data[:, i]
        result = bootstrap_s3i_test(vals; α=0.05)
        @printf("  %s: mean=%.4f, CI=[%.4f, %.4f], separated=%s\n",
                col, mean(vals), result.ci_lower, result.ci_upper,
                result.is_separated ? "✓" : "✗")
    end
end

println("\n" * "=" ^ 70)

#!/usr/bin/env julia

"""
Publication-Quality CTPL Noise Model Figure

This script generates a comprehensive 3-panel figure illustrating:
- Panel A: PDF comparison (Lévy, tempered Lévy, calibrated CTPL)
- Panel B: Tail dominance (CTPL vs Gaussian vs Gamma)
- Panel C: Auto-tuned λ(c) for different tolerance levels

Layout: 2 rows — A spans full width on top, B and C side by side below.
Color palette: Okabe-Ito (colorblind-safe).

Author: Ben Cardoen
Date: 2025-12-16
"""

# Activate shared reproduction environment
import Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using CalibratedTemperedPowerLaw
using CairoMakie
using Statistics
using Distributions
using Random
using Printf

# Okabe-Ito colorblind-safe palette (Bang Wong, Nature Methods 2011)
const OI = [
    colorant"#E69F00",  # orange
    colorant"#56B4E9",  # sky blue
    colorant"#009E73",  # bluish green
    colorant"#F0E442",  # yellow
    colorant"#0072B2",  # blue
    colorant"#D55E00",  # vermillion
    colorant"#CC79A7",  # reddish purple
    colorant"#000000",  # black
]

"""
    create_ctpl_figure(; r_max=40.0, δ_target=0.05, c_range=[1.0, 2.0, 3.0, 4.0, 5.0])

Create comprehensive CTPL noise model figure with 2-row layout.
Row 1: Panel A (full width). Row 2: Panels B and C side by side.
"""
function create_ctpl_figure(; r_max=40.0, δ_target=0.05, c_range=[1.0, 2.0, 3.0, 4.0, 5.0])

    Random.seed!(42)

    fig = Figure(size=(1800, 1200), fontsize=18)

    # ========================================================================
    # Panel A: PDF Comparison — full width, top row
    # ========================================================================

    ax_pdf = Axis(fig[1, 1:2],
        xlabel = "Displacement r",
        ylabel = "Probability Density f(r)",
        title = "A: PDF Comparison",
        titlesize = 20,
        xscale = log10,
        yscale = log10,
        xticks = [0.1, 1, 10, 100],
        limits = (0.1, 200, 1e-6, 10)
    )

    r_vals = exp10.(range(log10(0.1), log10(200), length=500))

    # Lévy (λ=0)
    c_ref = 1.0
    levy_vals = [levy_pdf(r, c_ref, 0.0) for r in r_vals]
    lines!(ax_pdf, r_vals, levy_vals,
           label="Lévy (λ=0)",
           color=OI[5],       # blue
           linewidth=2.5,
           linestyle=:dash)

    # Tempered Lévy with uncalibrated λ
    λ_uncal = [0.05, 0.1, 0.2]
    tempered_colors = [OI[1], OI[6], OI[7]]  # orange, vermillion, reddish purple
    for (i, λ) in enumerate(λ_uncal)
        temp_vals = [tempered_levy_pdf(r, c_ref, λ, 0.0) for r in r_vals]
        lines!(ax_pdf, r_vals, temp_vals,
               label="Tempered (λ=$(λ))",
               color=tempered_colors[i],
               linewidth=2.0,
               linestyle=:dot)
    end

    # Calibrated CTPL for different c values
    ctpl_colors = [OI[3], OI[2], OI[4]]  # bluish green, sky blue, yellow
    for (i, c) in enumerate(c_range[1:3])
        λ_cal = autotune_lambda(c, r_max, δ_target)
        ctpl_vals = [tempered_levy_pdf(r, c, λ_cal, 0.0) for r in r_vals]
        lines!(ax_pdf, r_vals, ctpl_vals,
               label="CTPL (c=$(c), λ≈$(@sprintf("%.3f", λ_cal)))",
               color=ctpl_colors[i],
               linewidth=2.5)
    end

    # r_max marker
    vlines!(ax_pdf, [r_max], color=OI[8], linestyle=:dashdot, linewidth=1.5, label="r_max = $(r_max)")

    # Legend inside plot, top-right (data drops off there on log-log scale)
    axislegend(ax_pdf, position=:rt, framevisible=true, labelsize=14, rowgap=1)

    # Annotations
    text!(ax_pdf, 1.5, 0.3, text="Power-law\nregime ∝ r⁻³ᐟ²",
          fontsize=16, color=:gray30)
    text!(ax_pdf, 50, 1e-5, text="Exponential\ntempering",
          fontsize=16, color=:gray30)

    # ========================================================================
    # Panel B: Tail Dominance — bottom left
    # ========================================================================

    ax_tail = Axis(fig[2, 1],
        xlabel = "Displacement r",
        ylabel = "Survival Function P(X > r)",
        title = "B: Tail Dominance",
        titlesize = 20,
        yscale = log10,
        limits = (0, 100, 1e-5, 1)
    )

    c_demo = 1.0
    λ_demo = autotune_lambda(c_demo, r_max, δ_target)

    n_samples = 5000
    ctpl_samples = [sample_tempered_levy(c_demo, λ_demo; μ=0.0, limit=r_max*3) for _ in 1:n_samples]
    ctpl_mean = mean(ctpl_samples)
    ctpl_var = var(ctpl_samples)

    gauss_dist = Normal(ctpl_mean, sqrt(ctpl_var))
    gamma_shape = ctpl_mean^2 / ctpl_var
    gamma_scale = ctpl_var / ctpl_mean
    gamma_dist = Gamma(gamma_shape, gamma_scale)

    r_tail = range(0, 100, length=300)

    # Use ccdf for numerical stability; clamp empirical survival away from 0 for log scale.
    @inline safe_survival(p::Real) = p > 0 ? float(p) : eps(Float64)

    ctpl_survival  = [safe_survival(mean(ctpl_samples .> r)) for r in r_tail]
    gauss_survival = [safe_survival(ccdf(gauss_dist, r)) for r in r_tail]
    gamma_survival = [safe_survival(ccdf(gamma_dist, r)) for r in r_tail]

    lines!(ax_tail, r_tail, ctpl_survival,
           label="CTPL (calibrated)",
           color=OI[5],       # blue
           linewidth=2.5)

    lines!(ax_tail, r_tail, gauss_survival,
           label="Gaussian (matched moments)",
           color=OI[1],       # orange
           linewidth=2.0,
           linestyle=:dash)

    lines!(ax_tail, r_tail, gamma_survival,
           label="Gamma (matched moments)",
           color=OI[3],       # bluish green
           linewidth=2.0,
           linestyle=:dot)

    hlines!(ax_tail, [δ_target], color=OI[8], linestyle=:dashdot, linewidth=1.5)
    text!(ax_tail, 65, δ_target*2.0, text="Target δ = $(δ_target)", fontsize=16)

    vlines!(ax_tail, [r_max], color=OI[8], linestyle=:dashdot, linewidth=1.5)

    axislegend(ax_tail, position=:rt, framevisible=true, labelsize=16)

    text!(ax_tail, 50, 5e-3, text="CTPL maintains\nheavy tail",
          fontsize=16, color=OI[5])

    # ========================================================================
    # Panel C: Auto-Tuned λ(c) — bottom right
    # ========================================================================

    ax_lambda = Axis(fig[2, 2],
        xlabel = "Scale parameter c",
        ylabel = "Calibrated tempering λ",
        title = "C: Auto-Tuned λ(c)",
        titlesize = 20,
    )

    δ_levels = [0.01, 0.05, 0.10]
    δ_colors = [OI[5], OI[1], OI[3]]  # blue, orange, bluish green

    c_vals = range(1.0, 5.0, length=9)
    n_runs_per_point = 20

    all_upper_bounds = Float64[]

    for (idx, δ) in enumerate(δ_levels)
        λ_means = Float64[]
        λ_stds = Float64[]

        println("Computing λ(c) for δ = $(δ) with $(n_runs_per_point) runs per point...")
        for (i, c) in enumerate(c_vals)
            λ_samples = Float64[]
            for run in 1:n_runs_per_point
                seed_val = 1000*idx + 100*i + run
                Random.seed!(seed_val)
                λ = autotune_lambda(c, r_max, δ)
                push!(λ_samples, λ)
            end

            μ = mean(λ_samples)
            σ = std(λ_samples)
            push!(λ_means, μ)
            push!(λ_stds, σ)
            push!(all_upper_bounds, μ + σ)
        end

        lines!(ax_lambda, c_vals, λ_means,
               label="δ = $(δ)",
               color=δ_colors[idx],
               linewidth=2.0)

        errorbars!(ax_lambda, c_vals, λ_means, λ_stds,
                   color=δ_colors[idx],
                   linewidth=1.5,
                   whiskerwidth=8)

        scatter!(ax_lambda, c_vals, λ_means,
                color=δ_colors[idx],
                markersize=10,
                marker=:circle,
                strokewidth=1.5,
                strokecolor=:white)
    end

    y_max = maximum(all_upper_bounds) * 1.1
    ylims!(ax_lambda, 0, y_max)

    axislegend(ax_lambda, position=:lt, framevisible=true, labelsize=16)

    text!(ax_lambda, 1.2, y_max * 0.15,
          text="P(X > $(r_max)) ≤ δ",
          fontsize=16,
          color=:gray30)

    # ========================================================================
    # Final adjustments
    # ========================================================================

    colgap!(fig.layout, 20)
    rowgap!(fig.layout, 20)

    return fig
end

# ============================================================================
# Generate and save figure
# ============================================================================

println("=" ^ 70)
println("Generating CTPL Noise Model Figure")
println("=" ^ 70)

fig = create_ctpl_figure(
    r_max = 40.0,
    δ_target = 0.05,
    c_range = [1.0, 2.0, 3.0, 4.0, 5.0]
)

# Save as PDF
output_path = joinpath(@__DIR__, "Figure2.pdf")
println("\nSaving figure to: $(output_path)")
save(output_path, fig)

println("\nFigure saved successfully!")
println("=" ^ 70)

# Also save PNG for quick preview
png_path = joinpath(@__DIR__, "Figure2.png")
save(png_path, fig, px_per_unit=2)
println("Preview PNG saved to: $(png_path)")

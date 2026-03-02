#!/usr/bin/env julia

"""
Figure 4: Vertex Noise Induces Correlated Edge Perturbations

Demonstrates that i.i.d. vertex noise → non-i.i.d. edge noise
using an asymmetric 3-spoke hub motif at k = 2r.

Layout: 3 rows
  - Row 1: Panel A (strong oracle) + Panel B (weak oracle)
  - Row 2: Correlation heatmaps varying c (1, 3, 5) with ε=0
  - Row 3: Correlation heatmaps varying ε (0, 5, 10) with c=1

Color palette: Okabe-Ito (colorblind-safe).
Saves CSV with all correlation matrix values for auditability.

Author: Ben Cardoen
Source: scripts/graph_noise_model_figure.jl
"""

import Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using CalibratedTemperedPowerLaw
using CairoMakie
using Statistics
using Random
using LinearAlgebra
using StaticArrays
using Printf

# Okabe-Ito palette
const OI_BLUE    = colorant"#0072B2"
const OI_ORANGE  = colorant"#E69F00"
const OI_GREEN   = colorant"#009E73"
const OI_VERMIL  = colorant"#D55E00"
const OI_PURPLE  = colorant"#CC79A7"
const OI_BLACK   = colorant"#000000"

function create_graph_noise_figure(; k=100.0, r=50.0, n_samples=500, δ=0.05)

    Random.seed!(42)

    fig = Figure(size=(1800, 1650), fontsize=16)

    # Hub at origin, 3 spokes at angles: 0°, 60°, 210° (asymmetric)
    hub_pos = SVector{2,Float64}(0.0, 0.0)
    spoke_angles = [0.0, π/3, π/3 + 5π/6]  # 0°, 60°, 210°
    spoke_positions = [SVector{2,Float64}(k * cos(θ), k * sin(θ)) for θ in spoke_angles]

    spoke_colors = [OI_ORANGE, OI_GREEN, OI_VERMIL]

    println("Graph parameters:")
    println("  Hub position: $(hub_pos)")
    println("  3 spokes at angles: 0°, 60°, 210° (asymmetric)")
    println("  Edge length k = $(k) = 2r")
    println("  Noise radius r = $(r)")

    # ========================================================================
    # Row 1: Panel A (strong oracle) + Panel B (weak oracle)
    # ========================================================================

    # Use nested grid for row 1 to get equal-width panels
    row1 = fig[1, 1:3] = GridLayout()

    ax_strong = Axis(row1[1, 1],
        xlabel = "x", ylabel = "y",
        title = "A: Strong Oracle (True Positions)",
        titlesize = 18,
        aspect = DataAspect()
    )

    # True edges
    for spoke in spoke_positions
        lines!(ax_strong, [hub_pos[1], spoke[1]], [hub_pos[2], spoke[2]],
               color=OI_BLACK, linewidth=3)
    end
    lines!(ax_strong, [NaN], [NaN], color=OI_BLACK, linewidth=3, label="True edges")

    # Hub
    scatter!(ax_strong, [hub_pos[1]], [hub_pos[2]],
            color=OI_BLUE, markersize=25, marker=:star5, label="Hub")

    # Spokes
    for (i, spoke) in enumerate(spoke_positions)
        scatter!(ax_strong, [spoke[1]], [spoke[2]],
                color=spoke_colors[i], markersize=15, marker=:circle,
                label="Spoke $(i)")
    end

    # Uncertainty circles (radius r)
    θ_circle = range(0, 2π, length=100)
    x_c = hub_pos[1] .+ r .* cos.(θ_circle)
    y_c = hub_pos[2] .+ r .* sin.(θ_circle)
    lines!(ax_strong, x_c, y_c, color=:gray, linestyle=:dash, linewidth=1.5)
    for spoke in spoke_positions
        x_c = spoke[1] .+ r .* cos.(θ_circle)
        y_c = spoke[2] .+ r .* sin.(θ_circle)
        lines!(ax_strong, x_c, y_c, color=:gray, linestyle=:dash, linewidth=1.5)
    end

    text!(ax_strong, k/2, k*0.1,
          text="k = $(k) = 2r", fontsize=14, align=(:center, :bottom))
    text!(ax_strong, hub_pos[1] + r/√2, hub_pos[2] + r/√2,
          text="r = $(r)", fontsize=12, color=:gray30)

    axislegend(ax_strong, position=:lt, framevisible=true, labelsize=12)

    # Panel B: Weak Oracle
    ax_weak = Axis(row1[1, 2],
        xlabel = "x", ylabel = "y",
        title = "B: Weak Oracle (Noisy Observations)",
        titlesize = 18,
        aspect = DataAspect()
    )

    # True positions reference
    for spoke in spoke_positions
        lines!(ax_weak, [hub_pos[1], spoke[1]], [hub_pos[2], spoke[2]],
               color=OI_BLACK, linewidth=2, alpha=0.5)
    end
    scatter!(ax_weak, [hub_pos[1]], [hub_pos[2]],
            color=OI_BLACK, markersize=15, marker=:star5, alpha=0.5, label="True positions")

    # Noisy samples (same seed logic as original for first 30, extended to 50)
    c_viz = 1.0
    λ_viz = autotune_lambda(c_viz, r, δ; mu=0.0, tol=1e-3, max_iter=50, n_samples=5000)
    n_overlay = 50

    # Add legend entries for noisy elements (plot invisible markers first)
    spoke_labels = ["E1 (noisy)", "E2 (noisy)", "E3 (noisy)"]
    for (j, lbl) in enumerate(spoke_labels)
        lines!(ax_weak, [NaN], [NaN], color=spoke_colors[j], linewidth=2, label=lbl)
    end
    scatter!(ax_weak, [NaN], [NaN],
            color=OI_BLUE, markersize=8, marker=:star5, label="Hub (noisy)")

    for i in 1:n_overlay
        hub_noisy_raw = levy_perturb(hub_pos, λ_viz, c_viz; limit=r, seed=42+i*100)
        hub_noisy = SVector{2,Float64}(hub_noisy_raw)

        for (j, spoke) in enumerate(spoke_positions)
            spoke_noisy_raw = levy_perturb(spoke, λ_viz, c_viz; limit=r, seed=42+i*100+j)
            spoke_noisy = SVector{2,Float64}(spoke_noisy_raw)

            lines!(ax_weak, [hub_noisy[1], spoke_noisy[1]],
                           [hub_noisy[2], spoke_noisy[2]],
                   color=(spoke_colors[j], 0.3), linewidth=1.5)
        end

        scatter!(ax_weak, [hub_noisy[1]], [hub_noisy[2]],
                color=(OI_BLUE, 0.35), markersize=8, marker=:star5)
    end

    text!(ax_weak, hub_pos[1], hub_pos[2] - k*0.7,
          text="50 noisy samples",
          fontsize=12, align=(:center, :top), color=:gray30)

    axislegend(ax_weak, position=:lt, framevisible=true, labelsize=11)

    # ========================================================================
    # Generate correlation data
    # ========================================================================

    # Collect all correlation results for CSV
    csv_rows = Vector{NamedTuple}()

    # Row 2: Varying c (c=1,3,5) with ε=0
    c_values = [1.0, 3.0, 5.0]
    edge_data_c = Dict{Float64, Matrix{Float64}}()
    lambda_values_c = Dict{Float64, Float64}()

    for c_val in c_values
        λ = autotune_lambda(c_val, r, δ; mu=0.0, tol=1e-3, max_iter=50, n_samples=5000)
        lambda_values_c[c_val] = λ
        println("\nRow 2: c=$(c_val), λ≈$(@sprintf("%.4f", λ)), ε=0")

        edge_lengths = zeros(n_samples, 3)
        for i in 1:n_samples
            hub_noisy_raw = levy_perturb(hub_pos, λ, c_val; limit=r, seed=42+i*100+Int(c_val)*10000)
            hub_noisy = SVector{2,Float64}(hub_noisy_raw)

            for (j, spoke) in enumerate(spoke_positions)
                spoke_noisy_raw = levy_perturb(spoke, λ, c_val; limit=r, seed=42+i*100+j+Int(c_val)*10000)
                spoke_noisy = SVector{2,Float64}(spoke_noisy_raw)
                edge_lengths[i, j] = norm(spoke_noisy - hub_noisy)
            end
        end
        edge_data_c[c_val] = edge_lengths
    end

    # Row 3: Varying ε (ε=0,5,10) with c=1
    epsilon_values = [0.0, 5.0, 10.0]
    c_fixed = 1.0
    λ_fixed = autotune_lambda(c_fixed, r, δ; mu=0.0, tol=1e-3, max_iter=50, n_samples=5000)
    edge_data_eps = Dict{Float64, Matrix{Float64}}()

    for ε in epsilon_values
        k_current = 2 * r + ε
        spoke_positions_current = [SVector{2,Float64}(k_current * cos(θ), k_current * sin(θ)) for θ in spoke_angles]
        println("\nRow 3: ε=$(ε), k=$(k_current), c=$(c_fixed), λ≈$(@sprintf("%.4f", λ_fixed))")

        edge_lengths = zeros(n_samples, 3)
        for i in 1:n_samples
            hub_noisy_raw = levy_perturb(hub_pos, λ_fixed, c_fixed; limit=r, seed=42+i*100+Int(ε)*1000)
            hub_noisy = SVector{2,Float64}(hub_noisy_raw)

            for (j, spoke) in enumerate(spoke_positions_current)
                spoke_noisy_raw = levy_perturb(spoke, λ_fixed, c_fixed; limit=r, seed=42+i*100+j+Int(ε)*1000)
                spoke_noisy = SVector{2,Float64}(spoke_noisy_raw)
                edge_lengths[i, j] = norm(spoke_noisy - hub_noisy)
            end
        end
        edge_data_eps[ε] = edge_lengths
    end

    # ========================================================================
    # Row 2: Heatmaps varying c
    # ========================================================================

    Label(fig[2, 1:3],
          text="Varying c (c=1, 3, 5) with ε=0, k=2r",
          fontsize=17, font=:bold, valign=:bottom)

    for (idx, c_val) in enumerate(c_values)
        λ_val = lambda_values_c[c_val]
        ax_corr = Axis(fig[3, idx],
            xlabel = "Edge", ylabel = "Edge",
            title = "c=$(c_val), λ≈$(@sprintf("%.4f", λ_val))",
            titlesize = 16,
            aspect = DataAspect()
        )

        corr_matrix = cor(edge_data_c[c_val])
        heatmap!(ax_corr, corr_matrix, colormap=:RdBu, colorrange=(-1, 1))

        ax_corr.xticks = (1:3, ["E1", "E2", "E3"])
        ax_corr.yticks = (1:3, ["E1", "E2", "E3"])

        for i in 1:3, j in 1:3
            text!(ax_corr, j, i,
                  text=@sprintf("%.2f", corr_matrix[i, j]),
                  align=(:center, :center),
                  color=(abs(corr_matrix[i, j]) > 0.5 ? :white : :black),
                  fontsize=14)
            push!(csv_rows, (condition="vary_c", param=c_val, epsilon=0.0,
                             row_edge="E$i", col_edge="E$j",
                             correlation=corr_matrix[i, j]))
        end
    end

    # ========================================================================
    # Row 3: Heatmaps varying ε
    # ========================================================================

    Label(fig[4, 1:3],
          text="Varying ε (ε=0, 5, 10) with c=$(c_fixed), λ≈$(@sprintf("%.4f", λ_fixed))",
          fontsize=17, font=:bold, valign=:bottom)

    for (idx, ε) in enumerate(epsilon_values)
        k_test = 2 * r + ε
        ax_corr = Axis(fig[5, idx],
            xlabel = "Edge", ylabel = "Edge",
            title = "ε=$(ε), k=$(k_test)",
            titlesize = 16,
            aspect = DataAspect()
        )

        corr_matrix = cor(edge_data_eps[ε])
        heatmap!(ax_corr, corr_matrix, colormap=:RdBu, colorrange=(-1, 1))

        ax_corr.xticks = (1:3, ["E1", "E2", "E3"])
        ax_corr.yticks = (1:3, ["E1", "E2", "E3"])

        for i in 1:3, j in 1:3
            text!(ax_corr, j, i,
                  text=@sprintf("%.2f", corr_matrix[i, j]),
                  align=(:center, :center),
                  color=(abs(corr_matrix[i, j]) > 0.5 ? :white : :black),
                  fontsize=14)
            push!(csv_rows, (condition="vary_eps", param=c_fixed, epsilon=ε,
                             row_edge="E$i", col_edge="E$j",
                             correlation=corr_matrix[i, j]))
        end
    end

    # ========================================================================
    # Final adjustments
    # ========================================================================

    rowsize!(fig.layout, 2, Relative(0.03))
    rowsize!(fig.layout, 4, Relative(0.03))
    rowgap!(fig.layout, 10)
    colgap!(fig.layout, 15)

    # ========================================================================
    # Save CSV
    # ========================================================================

    csv_path = joinpath(@__DIR__, "correlations.csv")
    open(csv_path, "w") do io
        println(io, "condition,param,epsilon,row_edge,col_edge,correlation")
        for r in csv_rows
            @printf(io, "%s,%.1f,%.1f,%s,%s,%.6f\n",
                    r.condition, r.param, r.epsilon, r.row_edge, r.col_edge, r.correlation)
        end
    end
    println("\nCorrelation data saved to: $csv_path")

    return fig, edge_data_c, edge_data_eps
end

# ============================================================================
# Run
# ============================================================================

println("=" ^ 70)
println("Generating Figure 4: Vertex Noise → Correlated Edge Perturbations")
println("=" ^ 70)

fig, edge_data_c, edge_data_eps = create_graph_noise_figure(
    k = 100.0, r = 50.0, n_samples = 500, δ = 0.05
)

# Print correlation analysis
println("\n" * "=" ^ 70)
println("Edge Correlation Analysis")
println("=" ^ 70)

println("\nRow 2: Varying c (ε=0):")
for c_val in [1.0, 3.0, 5.0]
    corr = cor(edge_data_c[c_val])
    off = [corr[i, j] for i in 1:3, j in 1:3 if i != j]
    println("  c=$(c_val): Max |ρ|=$(@sprintf("%.3f", maximum(abs.(off)))), Mean |ρ|=$(@sprintf("%.3f", mean(abs.(off))))")
end

println("\nRow 3: Varying ε (c=1.0):")
for ε in [0.0, 5.0, 10.0]
    corr = cor(edge_data_eps[ε])
    off = [corr[i, j] for i in 1:3, j in 1:3 if i != j]
    println("  ε=$(ε) (k=$(2*50.0+ε)): Max |ρ|=$(@sprintf("%.3f", maximum(abs.(off)))), Mean |ρ|=$(@sprintf("%.3f", mean(abs.(off))))")
end

output_pdf = joinpath(@__DIR__, "Figure4.pdf")
output_png = joinpath(@__DIR__, "Figure4.png")
save(output_pdf, fig)
save(output_png, fig, px_per_unit=2)

println("\nSaved: $output_pdf")
println("Saved: $output_png")
println("=" ^ 70)

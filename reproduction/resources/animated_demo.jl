#!/usr/bin/env julia
"""
    animated_demo.jl

Generate an animated GIF demonstrating CTPL noise on geometric graphs and SC convergence.

Left panel: Graph with accumulated noisy overlays (alpha blending creates blur effect)
Right panel: SC strong and weak convergence with variance bands

Usage:
    julia --project=reproduction reproduction/resources/animated_demo.jl
"""

import Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using CalibratedTemperedPowerLaw
using GLMakie
using DynamicGeometricGraphs
using Graphs
using Random
using Statistics
using Printf

# Configuration
const N_SAMPLES = 100
const GRAPH_SIZE = 100.0
const N_VERTICES = 30
const NOISE_C = 5.0
const NOISE_LIMIT = 40.0
const ALPHA_OVERLAY = 0.2  # Transparency for overlaid graphs
const SEED = 42

println("CTPL Animated Demo")
println("="^50)

# Set random seed
Random.seed!(SEED)

# Create reference graph
println("Creating reference graph...")
g_ref = refgraph(GRAPH_SIZE, N_VERTICES)
n_vertices = nv(g_ref)
println("  Graph has $n_vertices vertices")

# Auto-tune lambda for desired tail probability
println("Auto-tuning noise parameter λ...")
target_prob = 0.05  # 5% tail probability
λ = autotune_lambda(NOISE_C, NOISE_LIMIT, target_prob)
println("  λ = $(round(λ, digits=4)) (c=$NOISE_C, limit=$NOISE_LIMIT)")

# Pre-compute reference spectrum
println("Computing reference spectrum...")
A_ref = adjacency_matrix(g_ref)
_, λ_ref, _ = supra_spectrum(A_ref; mode=:comb, return_vectors=true)
ref_lambda = maximum(λ_ref)

# Storage for SC tracking
sc_strong_samples = Float64[]  # d(G_ref, G_noisy_i)
sc_weak_distances = Float64[]  # All pairwise d(G_noisy_i, G_noisy_j)
perturbed_spectra = Vector{Vector{Float64}}()

# Statistics tracking
sc_strong_means = Float64[]
sc_strong_stds = Float64[]
sc_weak_means = Float64[]
sc_weak_stds = Float64[]

println("\nGenerating $N_SAMPLES noisy samples...")
print("Progress: ")
flush(stdout)

# Create figure with two panels
fig = Figure(size=(1600, 800))

# Left panel: Graph visualization
ax_graph = Axis3(fig[1, 1],
    title="CTPL Noisy Graph (α=$ALPHA_OVERLAY per sample)",
    xlabel="X", ylabel="Y", zlabel="Z",
    aspect=(1, 1, 1))

# Right panel: SC convergence
ax_sc = Axis(fig[1, 2],
    title="SC Convergence (Strong vs Weak)",
    xlabel="Number of Samples",
    ylabel="Spectral Distance (Wasserstein-2)",
    xgridvisible=true,
    ygridvisible=true)

# Plot ground truth graph (once, solid)
println("  Plotting ground truth graph...")
# Manually plot since extension might not load properly
vertices_ref = collect(Graphs.vertices(g_ref))
n_ref = length(vertices_ref)
crds_ref = zeros(n_ref, 3)
v2v_ref = Dict{Int, Int}()

for (idx, vi) in enumerate(vertices_ref)
    v2v_ref[vi] = idx
    coords = DynamicGeometricGraphs.get_vertex_coords(g_ref, vi)
    crds_ref[idx, 1:2] .= coords[1], coords[2]
    crds_ref[idx, 3] = 0.0
end

# Plot vertices
scatter!(ax_graph, crds_ref[:, 1], crds_ref[:, 2], crds_ref[:, 3],
    color=:red, markersize=15, label="Ground Truth")

# Plot edges
for e in Graphs.edges(g_ref)
    ui = Graphs.src(e)
    vi = Graphs.dst(e)
    if haskey(v2v_ref, ui) && haskey(v2v_ref, vi)
        x, y, z = crds_ref[v2v_ref[ui], :]
        xp, yp, zp = crds_ref[v2v_ref[vi], :]
        lines!(ax_graph, [x, xp], [y, yp], [z, zp],
            color=:darkred, linewidth=2)
    end
end

# Animation frames
frames = 1:N_SAMPLES

# Progress bar helper
function progress_bar(i, n)
    bar_width = 40
    progress = i / n
    filled = round(Int, bar_width * progress)
    bar = "["  * "=" ^ filled * " " ^ (bar_width - filled) * "]"
    percent = round(progress * 100, digits=1)
    print("\r  $bar $percent% ($i/$n)")
    flush(stdout)
end

record(fig, "ctpl_demo.gif", frames; framerate=10) do i
    progress_bar(i, N_SAMPLES)

    # Generate noisy graph
    g_noisy, _, _ = perturb_graph(g_ref; λ=λ, c=NOISE_C, limit=NOISE_LIMIT, seed=i)

    # Compute noisy spectrum and strong distance
    A_noisy = adjacency_matrix(g_noisy)
    _, λ_noisy, _ = supra_spectrum(A_noisy; mode=:comb, return_vectors=true)

    # Strong SC: distance to ground truth
    d_strong = w2_1d_empirical(λ_ref, λ_noisy; ref_lambda=ref_lambda)
    push!(sc_strong_samples, d_strong)
    push!(perturbed_spectra, λ_noisy)

    # Weak SC: pairwise distances between all noisy samples so far
    for j in 1:(i-1)
        d_weak = w2_1d_empirical(perturbed_spectra[j], λ_noisy; ref_lambda=ref_lambda)
        push!(sc_weak_distances, d_weak)
    end

    # Update statistics
    push!(sc_strong_means, mean(sc_strong_samples))
    push!(sc_strong_stds, i > 1 ? std(sc_strong_samples) : 0.0)

    if length(sc_weak_distances) > 0
        push!(sc_weak_means, mean(sc_weak_distances))
        push!(sc_weak_stds, length(sc_weak_distances) > 1 ? std(sc_weak_distances) : 0.0)
    else
        push!(sc_weak_means, 0.0)
        push!(sc_weak_stds, 0.0)
    end

    # Plot noisy graph with transparency (overlay on ground truth)
    # Get vertex coordinates
    vertices_list = collect(Graphs.vertices(g_noisy))
    n = length(vertices_list)
    crds = zeros(n, 3)
    v2v = Dict{Int, Int}()

    for (idx, vi) in enumerate(vertices_list)
        v2v[vi] = idx
        coords = DynamicGeometricGraphs.get_vertex_coords(g_noisy, vi)
        crds[idx, 1:2] .= coords[1], coords[2]
        crds[idx, 3] = 0.0  # 2D graph in 3D space
    end

    # Plot noisy vertices with transparency
    scatter!(ax_graph, crds[:, 1], crds[:, 2], crds[:, 3],
        color=(:blue, ALPHA_OVERLAY),
        markersize=8,
        label=(i == 1 ? "Noisy Samples" : ""))

    # Plot noisy edges with transparency
    for e in Graphs.edges(g_noisy)
        ui = Graphs.src(e)
        vi = Graphs.dst(e)

        if haskey(v2v, ui) && haskey(v2v, vi)
            x, y, z = crds[v2v[ui], :]
            xp, yp, zp = crds[v2v[vi], :]

            lines!(ax_graph, [x, xp], [y, yp], [z, zp],
                color=(:cyan, ALPHA_OVERLAY),
                linewidth=1)
        end
    end

    # Add legend to graph panel (first frame only)
    if i == 1
        axislegend(ax_graph, position=:lt)
    end

    # Update SC convergence plot
    empty!(ax_sc)

    x_vals = 1:i

    # Plot SC strong with error bars
    lines!(ax_sc, x_vals, sc_strong_means,
        color=:red, linewidth=2, label="SC Strong (Oracle)")
    band!(ax_sc, x_vals,
        sc_strong_means .- sc_strong_stds,
        sc_strong_means .+ sc_strong_stds,
        color=(:red, 0.2))

    # Plot SC weak with error bars (if we have data)
    if length(sc_weak_means) > 0 && sc_weak_means[end] > 0
        lines!(ax_sc, x_vals, sc_weak_means,
            color=:blue, linewidth=2, label="SC Weak (Observable)")
        band!(ax_sc, x_vals,
            sc_weak_means .- sc_weak_stds,
            sc_weak_means .+ sc_weak_stds,
            color=(:blue, 0.2))
    end

    # Add legend (show after we have both strong and weak data)
    if i >= 2
        axislegend(ax_sc, position=:rt)
    end

    # Update title with current statistics
    if i > 1
        ax_sc.title = @sprintf("SC Convergence | Strong: %.4f±%.4f | Weak: %.4f±%.4f",
            sc_strong_means[end], sc_strong_stds[end],
            length(sc_weak_means) > 0 ? sc_weak_means[end] : 0.0,
            length(sc_weak_stds) > 0 ? sc_weak_stds[end] : 0.0)
    end
end

println("\n\n✓ Animation complete!")
println("  Saved to: ctpl_demo.gif")
println("\nFinal Statistics:")
println("  SC Strong: $(round(sc_strong_means[end], digits=4)) ± $(round(sc_strong_stds[end], digits=4))")
if length(sc_weak_means) > 0
    println("  SC Weak:   $(round(sc_weak_means[end], digits=4)) ± $(round(sc_weak_stds[end], digits=4))")
    println("  Ratio (Weak/Strong): $(round(sc_weak_means[end] / sc_strong_means[end], digits=2))")
end
println("\nNote: SC Weak ≤ 2×SC Strong by triangle inequality")

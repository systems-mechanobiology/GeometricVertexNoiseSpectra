#!/usr/bin/env julia
"""
Motif Spectral Tightness Table Generator

Combines:
- CORRECT geometric construction from patterns_minfix2.jl
  (hub moves along bisector, spokes move away from displaced hub h')
- Spectral analysis (Laplacian eigenvalues, Weyl bound)

Outputs the table for paper: degree, δWD, actual_shift, weyl_limit, ratio

Author: Generated for GTVN paper
Date: 2025-01-16
"""

import Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using LinearAlgebra
using Printf

# =============================================================================
# CORRECT GEOMETRIC CONSTRUCTION (from patterns_minfix2.jl)
# =============================================================================

"""
    canonical_slots(d, k)

Generate d spoke positions at 60° intervals (clockwise), all at distance k from origin.
Spokes are at angles: 0°, -60°, -120°, ..., -(d-1)*60°
This matches the DynamicGeometricGraphs.generate_hub_spoke_graph convention.
"""
function canonical_slots(d::Int, k::Float64)
    # Clockwise placement: 0, -60, -120, ... degrees
    angles = [-(i-1) * 60.0 * (π/180) for i in 1:d]
    [k .* [cos(a), sin(a)] for a in angles]
end

"""
    bisector_angle(d)

Compute the bisector angle of the unoccupied arc for degree d.
Spokes are placed at 0°, -60°, -120°, ..., -(d-1)*60° (clockwise from 0°).
The unoccupied arc is the remaining angular span, and the hub moves
along its bisector.

NOTE: For d=6, the hub still moves. The unoccupied arc is 60° (from -300° to 0°),
and the bisector is at -330° (equivalently 30°). This matches the correct
geometric model in scripts/patterns_new.jl.
"""
function bisector_angle(d::Int)
    # Last spoke angle: -(d-1) * 60° in radians
    last_spoke_angle = -(d - 1) * (π / 3)

    # Occupied arc spans from 0 to last_spoke_angle (absolute value)
    occupied_arc = abs(last_spoke_angle)

    # Unoccupied arc is the remainder
    unoccupied_arc = 2π - occupied_arc

    # Bisector: from last spoke, go halfway into unoccupied arc (clockwise = negative)
    bisector = last_spoke_angle - (unoccupied_arc / 2)

    return bisector
end

"""
    extremal_motif_geometry(d; k, r)

Compute the CORRECT extremal configuration for a degree-d hub-spoke motif.

Returns: (weights_before, weights_after, hub, hub_displaced, spokes, spokes_displaced)

Geometric rules:
1. Hub at origin, d spokes at distance k in consecutive 60° slots
2. Hub moves distance r along bisector of unoccupied arc → h'
3. Each spoke moves distance r away from h' (along direction k - h')
"""
function extremal_motif_geometry(d::Int; k::Float64, r::Float64)
    # Initial configuration
    hub = [0.0, 0.0]
    spokes = canonical_slots(d, k)

    # Edge weights before = distances from hub to spokes = k for all
    weights_before = [norm(s - hub) for s in spokes]

    # Compute hub displacement (hub always moves for all degrees 1-6)
    θ = bisector_angle(d)
    hub_displaced = hub .+ r .* [cos(θ), sin(θ)]

    # Compute spoke displacements: each moves away from h' by distance r
    spokes_displaced = Vector{Vector{Float64}}(undef, d)
    for i in 1:d
        s = spokes[i]
        direction = s - hub_displaced
        direction_normalized = direction / norm(direction)
        spokes_displaced[i] = s .+ r .* direction_normalized
    end

    # Edge weights after = distances from h' to displaced spokes
    weights_after = [norm(spokes_displaced[i] - hub_displaced) for i in 1:d]

    return (
        weights_before = weights_before,
        weights_after = weights_after,
        hub = hub,
        hub_displaced = hub_displaced,
        spokes = spokes,
        spokes_displaced = spokes_displaced
    )
end

# =============================================================================
# SPECTRAL ANALYSIS (from paper, these functions are correct)
# =============================================================================

"""
    star_adjacency(weights)

Build symmetric adjacency matrix for a star graph (hub + leaves).
Hub is vertex 1, leaves are vertices 2, 3, ..., d+1.
"""
function star_adjacency(weights::AbstractVector{<:Real})
    d = length(weights)
    n = d + 1
    A = zeros(Float64, n, n)
    for (j, w) in enumerate(weights)
        leaf = j + 1
        A[1, leaf] = w
        A[leaf, 1] = w
    end
    return Symmetric(A)
end

"""
    laplacian_from_weighted_adj(A)

Compute combinatorial Laplacian L = D - A from weighted adjacency matrix.
"""
function laplacian_from_weighted_adj(A::AbstractMatrix{<:Real})
    Ad = Array(A)
    D = Diagonal(sum(Ad, dims=2)[:])
    L = Matrix(D) - Ad
    return Symmetric((L + L') / 2)
end

"""
    spectral_shift(L0, L1)

Compute maximum eigenvalue shift between two Laplacians.
"""
function spectral_shift(L0::Symmetric, L1::Symmetric)
    λ0 = eigvals(L0)
    λ1 = eigvals(L1)
    @assert length(λ0) == length(λ1)
    return maximum(abs.(λ1 .- λ0))
end

"""
    weyl_bound(L0, L1)

Compute Weyl's bound ||E||_2 where E = L1 - L0.
"""
function weyl_bound(L0::Symmetric, L1::Symmetric)
    E = Matrix(L1) - Matrix(L0)
    return opnorm(E, 2)
end

# =============================================================================
# TABLE GENERATION
# =============================================================================

"""
    motif_tightness_row(d; k, r)

Compute one row of the spectral tightness table for degree d.
"""
function motif_tightness_row(d::Int; k::Float64, r::Float64)
    # Get correct geometry
    geom = extremal_motif_geometry(d; k=k, r=r)

    # Build Laplacians
    A0 = star_adjacency(geom.weights_before)
    A1 = star_adjacency(geom.weights_after)
    L0 = laplacian_from_weighted_adj(A0)
    L1 = laplacian_from_weighted_adj(A1)

    # Compute spectral quantities
    shift = spectral_shift(L0, L1)
    weyl = weyl_bound(L0, L1)
    ratio = shift / weyl

    # Also compute δWD for completeness
    δWD = sum(geom.weights_after) - sum(geom.weights_before)

    return (
        degree = d,
        δWD = δWD,
        δWD_over_r = δWD / r,
        actual_shift = shift,
        weyl_limit = weyl,
        ratio = ratio
    )
end

"""
    generate_table(; k, r)

Generate the full spectral tightness table for degrees 1-6.
"""
function generate_table(; k::Float64, r::Float64)
    println("=" ^ 80)
    println("Motif Spectral Tightness Table")
    println("Parameters: k = $k, r = $r (tight case: k = 2r)")
    println("=" ^ 80)
    println()

    # Print header
    @printf("%6s | %10s | %12s | %12s | %12s | %8s\n",
            "Degree", "δWD/r", "δWD", "Actual Shift", "Weyl ||E||₂", "Ratio")
    println("-" ^ 80)

    rows = []
    for d in 1:6
        row = motif_tightness_row(d; k=k, r=r)
        push!(rows, row)

        @printf("%6d | %10.4f | %12.6f | %12.6f | %12.6f | %8.6f\n",
                row.degree, row.δWD_over_r, row.δWD,
                row.actual_shift, row.weyl_limit, row.ratio)
    end

    println("-" ^ 80)
    println()

    # Print LaTeX table format
    println("LaTeX table rows:")
    println("-" ^ 80)
    for row in rows
        @printf("%d & \$%.2fr\$ & %.4f & %.4f & %.4f \\\\\n",
                row.degree, row.δWD_over_r, row.actual_shift, row.weyl_limit, row.ratio)
    end

    return rows
end

# =============================================================================
# MAIN
# =============================================================================

if abspath(PROGRAM_FILE) == @__FILE__
    # Tight case: k = 2r (ε → 0)
    r = 1.0
    k = 2.0 * r

    rows = generate_table(; k=k, r=r)

    println()
    println("=" ^ 80)
    println("Verification: Geometry details for each degree")
    println("=" ^ 80)

    for d in 1:6
        geom = extremal_motif_geometry(d; k=k, r=r)
        θ = bisector_angle(d)
        θ_deg = @sprintf("%.1f°", mod(θ * 180 / π, 360))

        println()
        println("Degree $d:")
        println("  Bisector angle: $θ_deg")
        println("  Hub displaced: $(round.(geom.hub_displaced, digits=4))")
        println("  Weights before: $(round.(geom.weights_before, digits=4))")
        println("  Weights after:  $(round.(geom.weights_after, digits=4))")
    end

    # Save markdown table to file
    outpath = joinpath(@__DIR__, "table1.md")
    open(outpath, "w") do io
        println(io, "| Degree | δWD ε→0 | Actual Shift | Weyl Limit | Ratio |")
        println(io, "|--------|---------|-------------|------------|-------|")
        for row in rows
            @printf(io, "| %d | \$%.2fr\$ | %.4f | %.4f | %.4f |\n",
                    row.degree, row.δWD_over_r, row.actual_shift, row.weyl_limit, row.ratio)
        end
    end
    println("\nMarkdown table saved to: $outpath")
end

#!/usr/bin/env julia

"""
Motif Pattern Generation with Geometric Perturbations

This script generates hub-spoke motif graphs with degrees 1-6 and applies
geometric perturbations by moving each vertex a specified distance in a
specified angular direction.

# Dependencies
- DynamicGeometricGraphs.jl
- LinearAlgebra, StaticArrays, Graphs, Printf

# Usage
    julia --project=reproduction reproduction/resources/patterns_extremal_motifs.jl

# Author
Generated for GTVN paper analysis
Date: 2025-12-15
"""

import Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
using DynamicGeometricGraphs
using LinearAlgebra
using StaticArrays
using Graphs
using Printf
const DDG = DynamicGeometricGraphs

"""
    perturb_motif(g::DynamicGeometricGraph{2, Float64}, angles::Vector{Float64}, r::Float64)

Apply geometric perturbations to a hub-spoke motif graph by moving each vertex
a distance `r` in the direction specified by the corresponding angle.

# Arguments
- `g::DynamicGeometricGraph{2, Float64}`: A hub-spoke motif graph
- `angles::Vector{Float64}`: Vector of angles in radians, one per vertex (length must equal nv(g))
- `r::Float64`: Movement distance for each vertex

# Returns
- A new `DynamicGeometricGraph` with vertices moved to new positions

# Algorithm
For each vertex i with coordinates (x, y):
1. Get the angle θ from angles[i]
2. Compute displacement: dx = r * cos(θ), dy = r * sin(θ)
3. New position: (x + dx, y + dy)

# Example
```julia
g = DDG.generate_hub_spoke_graph(100.0, 3)  # degree-3 motif
n = DDG.nv(g)
angles = [0.0, π/2, π, 3π/2]  # 4 vertices (hub + 3 spokes)
r = 50.0
g_perturbed = perturb_motif(g, angles, r)
```
"""
function perturb_motif(g::DynamicGeometricGraph{2, Float64}, angles::Vector{Float64}, r::Float64)
    n = DDG.nv(g)

    # Validate input
    if length(angles) != n
        error("Number of angles ($(length(angles))) must match number of vertices ($n)")
    end

    if r < 0
        error("Movement distance r must be non-negative, got $r")
    end

    # Create a copy of the graph
    g_new = DDG.freeze(g)

    # Get all vertices (sorted to ensure consistent ordering)
    vids = sort(collect(DDG.vertices(g)))

    # Apply perturbations to each vertex
    for (i, vid) in enumerate(vids)
        # Get current coordinates
        old_coords = DDG.get_vertex_coords(g, vid)

        # Calculate displacement vector
        θ = angles[i]
        dx = r * cos(θ)
        dy = r * sin(θ)
        displacement = SVector{2, Float64}(dx, dy)

        # Calculate new coordinates
        new_coords = old_coords + displacement

        # Update the vertex position
        DDG.update_coord!(g_new, old_coords, new_coords)

        @debug "Vertex $vid: moved from $old_coords to $new_coords (angle=$(rad2deg(θ))°, r=$r)"
    end

    return g_new
end


# ============================================================================
# Test the perturb_motif function
# ============================================================================

println("Testing perturb_motif function")
println("=" ^ 60)

# Test 1: Simple degree-1 motif
println("\nTest 1: Degree-1 motif (hub + 1 spoke)")
g1 = DDG.generate_hub_spoke_graph(100.0, 1)
println("  Original graph: $(DDG.nv(g1)) vertices, $(DDG.ne(g1)) edges")

# Get vertex info
vids = sort(collect(DDG.vertices(g1)))
println("  Vertex IDs: $vids")
for vid in vids
    coords = DDG.get_vertex_coords(g1, vid)
    println("    Vertex $vid: $coords")
end

# Define angles: move hub to the right (0°), spoke to the right (0°)
angles_test1 = [0.0, 0.0]  # Both move in same direction
r = 50.0

g1_perturbed = perturb_motif(g1, angles_test1, r)
println("\n  After perturbation (both at 0°, r=$r):")
for vid in vids
    coords = DDG.get_vertex_coords(g1_perturbed, vid)
    println("    Vertex $vid: $coords")
end

# Check weighted degree
wd_before = DDG.weighted_degree(g1, 1)
wd_after = DDG.weighted_degree(g1_perturbed, 1)
println("\n  Hub weighted degree before: $wd_before")
println("  Hub weighted degree after: $wd_after")
println("  Change: $(wd_after - wd_before)")


# Test 2: Degree-3 motif with varied angles
println("\n" * "=" ^ 60)
println("\nTest 2: Degree-3 motif with different angles")
g3 = DDG.generate_hub_spoke_graph(100.0, 3)
vids3 = sort(collect(DDG.vertices(g3)))

println("  Original vertices:")
for vid in vids3
    coords = DDG.get_vertex_coords(g3, vid)
    println("    Vertex $vid: $coords")
end

# Different angles for each vertex
angles_test2 = [0.0, π/6, π/3, π/2]  # Hub, spoke1, spoke2, spoke3
g3_perturbed = perturb_motif(g3, angles_test2, r)

println("\n  After perturbation (varied angles, r=$r):")
for (i, vid) in enumerate(vids3)
    coords = DDG.get_vertex_coords(g3_perturbed, vid)
    println("    Vertex $vid: $coords (angle=$(rad2deg(angles_test2[i]))°)")
end


# Test 3: All degrees 1-6
println("\n" * "=" ^ 60)
println("\nTest 3: Generate perturbed versions for all degrees 1-6")
println("  Using zero angles (no movement) as baseline")

for degree in 1:6
    g = DDG.generate_hub_spoke_graph(100.0, degree)
    n = DDG.nv(g)

    # All vertices move at angle 0 (to the right)
    angles_zero = zeros(Float64, n)
    g_perturbed = perturb_motif(g, angles_zero, r)

    println("\n  Degree $degree: $n vertices")
    println("    Original weighted degree (hub): $(DDG.weighted_degree(g, 1))")
    println("    Perturbed weighted degree (hub): $(DDG.weighted_degree(g_perturbed, 1))")
end

println("\n" * "=" ^ 60)
println("✓ All perturb_motif tests completed successfully!")


# ============================================================================
# Calculate bisector angles for extremal configurations
# ============================================================================

println("\n" * "=" ^ 60)
println("Computing bisector angles for extremal configurations")
println("=" ^ 60)

for d in 1:6
    if d == 6
        println("\nDegree $d: Hub does not move (special case)")
        continue
    end

    # Spokes are at angles: 0°, -60°, -120°, ..., -(d-1)×60°
    # First spoke: 0°
    # Last spoke: -(d-1) × 60°
    last_spoke_angle_deg = -(d-1) * 60

    # Occupied arc: from first spoke (0°) to last spoke
    occupied_arc_deg = abs(last_spoke_angle_deg)

    # Unoccupied arc: the remainder
    unoccupied_arc_deg = 360 - occupied_arc_deg

    # Bisector: middle of unoccupied arc
    # Starting from last spoke, go half the unoccupied arc (clockwise, so negative)
    bisector_angle_deg = last_spoke_angle_deg - (unoccupied_arc_deg / 2)

    # Convert to radians
    bisector_angle_rad = deg2rad(bisector_angle_deg)

    # Normalize to [0, 2π) for display
    bisector_normalized_deg = mod(bisector_angle_deg, 360)

    println("\nDegree $d:")
    println("  Spoke angles: 0°, -60°, ..., $(last_spoke_angle_deg)°")
    println("  Occupied arc: $(occupied_arc_deg)°")
    println("  Unoccupied arc: $(unoccupied_arc_deg)°")
    println("  Bisector angle: $(bisector_angle_deg)° (or $(bisector_normalized_deg)°)")
    println("  Bisector in radians: $(bisector_angle_rad)")
end

println("\n" * "=" ^ 60)

"""
    extremal_open_arc_push(g::DynamicGeometricGraph{2,Float64}, d::Int, r::Float64;
                           hub_vid::Int=1)

1) Compute the open (unoccupied) arc of the clockwise-filled spokes (as in patterns_new.jl).
2) Move hub h along the open-arc bisector by distance r to h'.
3) For each spoke endpoint k, move it along the direction (k - h') by distance r.

Assumes hub vertex id is `hub_vid` (default 1), consistent with weighted_degree(g, 1) usage.
"""
function extremal_open_arc_push(g::DynamicGeometricGraph{2,Float64}, d::Int, r::Float64; hub_vid::Int=1)
    r < 0 && error("r must be non-negative, got $r")

    # --- bisector angle (same construction as in your script) ---
    last_spoke_angle = -(d - 1) * (π/3)               # -(d-1)*60° in radians
    occupied_arc = abs(last_spoke_angle)
    unoccupied_arc = 2π - occupied_arc
    bisector = last_spoke_angle - (unoccupied_arc / 2)

    g_new = DDG.freeze(g)

    # --- move hub: h -> h' ---
    h = DDG.get_vertex_coords(g, hub_vid)
    hprime = h + SVector{2,Float64}(r*cos(bisector), r*sin(bisector))
    DDG.update_coord!(g_new, h, hprime)

    # --- move each spoke endpoint k along direction h' -> k by distance r ---
    for vid in sort(collect(DDG.vertices(g)))
        vid == hub_vid && continue
        k = DDG.get_vertex_coords(g, vid)  # original k
        v = k - hprime
        nv = norm(v)
        nv == 0 && continue  # degenerate; should not occur unless r is extreme
        kprime = k + (r / nv) * v
        DDG.update_coord!(g_new, k, kprime)
    end

    return g_new, bisector
end


k=100
r=k/2
d=3

graphs = Dict()
for d in 1:6
    g = DDG.generate_hub_spoke_graph(k, d)
    gp, bs = extremal_open_arc_push(g, d, r)
    A = adjacency_matrix(g)
    Ap = adjacency_matrix(gp)
    @info sum(Ap[1,:] .-  A[1,:])/r, bs
    graphs[d] = (g, gp)
end

g3, g3p = graphs[3]

f=plot_before_after_glmakie(g3, g3p)
f
using GLMakie
using Graphs


using Statistics
GLMakie.activate!()
Gs, Ms, f = run_suite(; lim=200)
f
f
Ms
CairoMakie.activate!()
save("displacements.pdf", f)


k = 100.0
r = k/2
d = 3
g = DDG.generate_hub_spoke_graph(k, d)

fig, info = plot_optimality_regret(g, d; r=r, δmax_deg=180, ngrid=721)
display(fig)


function run_suite(; k=100.0, r=k/2, ntrials=2000, σ_deg=10.0, seed=1,
                   outprefix="motifs", lim=150.0)

    rng = MersenneTwister(seed)
    σ = deg2rad(σ_deg)

    graphs  = Dict{Int,Tuple{Any,Any}}()
    summary = Dict{Int,NamedTuple}()

    fig = Figure(resolution=(1200, 700))

    for d in 1:6
        g = DDG.generate_hub_spoke_graph(k, d)
        gp, bis = extremal_open_arc_push(g, d, r)
        graphs[d] = (g, gp)

        wd0 = DDG.weighted_degree(g, 1)
        wdp = DDG.weighted_degree(gp, 1)

        deltas = Vector{Float64}(undef, ntrials)
        for t in 1:ntrials
            gη, _, _ = extremal_open_arc_push_noisy(g, d, r; σ=σ, rng=rng)
            deltas[t] = DDG.weighted_degree(gη, 1) - wd0
        end

        summary[d] = (k=k, r=r, d=d, bisector=bis, wd_before=wd0, wd_after=wdp,
                      delta_opt=wdp-wd0, delta_noise_mean=mean(deltas),
                      delta_noise_max=maximum(deltas), delta_noise_p99=quantile(deltas, 0.99))

        i, j = fldmod1(d, 3)
        ax = Axis(fig[i, j];
            aspect=DataAspect(),
            title="d=$d, bis=$(round(rad2deg(bis), digits=1))°",
            xlabel="x", ylabel="y"
        )

        plot_before_after!(ax, g, gp)

        # key change: force identical limits for all subplots
        xlims!(ax, -lim, lim)
        ylims!(ax, -lim, lim)
    end

    save("$(outprefix)_before_after.png", fig)
    return graphs, summary, fig
end


function plot_before_after_glmakie(g, gp; hub_vid::Int=1, title::AbstractString="", show_vids::Bool=false)
    vids = sort(collect(DDG.vertices(g)))

    coords(gx) = Dict(v => DDG.get_vertex_coords(gx, v) for v in vids)
    c0 = coords(g)
    c1 = coords(gp)

    toP(c::SVector{2,<:Real}) = Point2f(c[1], c[2])

    # edges -> line segments
    function segments(gx_coords)
        segs = Point2f[]
        for e in edges(g)
            u, v = src(e), dst(e)
            push!(segs, toP(gx_coords[u]))
            push!(segs, toP(gx_coords[v]))
        end
        return segs
    end

    fig = Figure()
    ax = Axis(fig[1, 1]; aspect=DataAspect(), title=title, xlabel="x", ylabel="y")

    # edges
    linesegments!(ax, segments(c0); linewidth=2)                  # before (solid)
    linesegments!(ax, segments(c1); linewidth=2, linestyle=:dash) # after  (dashed)

    # vertices
    pts0 = [toP(c0[v]) for v in vids if v != hub_vid]
    pts1 = [toP(c1[v]) for v in vids if v != hub_vid]
    scatter!(ax, pts0; markersize=10)
    scatter!(ax, pts1; marker=:x, markersize=10)

    # hub highlight
    scatter!(ax, [toP(c0[hub_vid])]; marker=:star5, markersize=18)
    scatter!(ax, [toP(c1[hub_vid])]; marker=:star5, markersize=18)

    if show_vids
        for v in vids
            p = toP(c0[v])
            text!(ax, string(v); position=p, align=(:left, :bottom))
        end
    end

    return fig
end



# --------- noise model: perturb hub-direction only ----------
# (so we are testing the "bisector is optimal hub direction" claim directly)
function extremal_open_arc_push_noisy(g::DynamicGeometricGraph{2,Float64}, d::Int, r::Float64;
                                      hub_vid::Int=1, σ::Float64=0.0, rng=Random.default_rng())
    # bisector from your construction
    last = -(d - 1) * (π/3)
    occ  = abs(last)
    unocc = 2π - occ
    bis  = last - unocc/2

    θ = bis + (σ == 0 ? 0.0 : randn(rng) * σ)   # noisy hub direction

    g_new = DDG.freeze(g)

    h = DDG.get_vertex_coords(g, hub_vid)
    h′ = h + SVector{2,Float64}(r*cos(θ), r*sin(θ))
    DDG.update_coord!(g_new, h, h′)

    # spokes: deterministically away from h′ (your rule)
    for vid in sort(collect(DDG.vertices(g)))
        vid == hub_vid && continue
        k = DDG.get_vertex_coords(g, vid)
        v = k - h′
        nv = norm(v)
        nv == 0 && continue
        k′ = k + (r/nv) * v
        DDG.update_coord!(g_new, k, k′)
    end

    return g_new, bis, θ
end
using Random
# --------- run: generate, plot, save, noise-test ----------
function run_suite(; k=100.0, r=k/2, ntrials=2000, σ_deg=10.0, seed=1, outprefix="motifs")
    rng = MersenneTwister(seed)
    σ = deg2rad(σ_deg)

    graphs  = Dict{Int,Tuple{Any,Any}}()
    summary = Dict{Int,NamedTuple}()

    ps = Plots.Plot[]
    for d in 1:6
        g = DDG.generate_hub_spoke_graph(k, d)

        gp, bis = extremal_open_arc_push(g, d, r)  # your deterministic “optimal”
        graphs[d] = (g, gp)

        wd0 = DDG.weighted_degree(g, 1)
        wdp = DDG.weighted_degree(gp, 1)

        # noise test: vary hub direction around bisector, keep everything else same
        deltas = Vector{Float64}(undef, ntrials)
        for t in 1:ntrials
            gη, _, _ = extremal_open_arc_push_noisy(g, d, r; σ=σ, rng=rng)
            deltas[t] = DDG.weighted_degree(gη, 1) - wd0
        end

        summary[d] = (k=k, r=r, d=d, bisector=bis, wd_before=wd0, wd_after=wdp,
                      delta_opt=wdp-wd0, delta_noise_mean=mean(deltas),
                      delta_noise_max=maximum(deltas), delta_noise_p99=quantile(deltas, 0.99))

        push!(ps, overlay_plot(g, gp; title_str="d=$d, bis=$(round(rad2deg(bis), digits=1))°"))
    end

    bigplot = plot(ps..., layout=(2,3), size=(1200,700))
    savefig(bigplot, "$(outprefix)_before_after.png")

    # @save "$(outprefix)_graphs_and_summary.jld2" graphs summary

    return graphs, summary
end

using Random, Statistics, LinearAlgebra, StaticArrays
using GLMakie
using Graphs
const DDG = DynamicGeometricGraphs

# in-place Makie plotter onto an existing Axis
function plot_before_after!(ax::Axis, g, gp; hub_vid::Int=1, show_vids::Bool=false)
    vids = sort(collect(DDG.vertices(g)))
    coords(gx) = Dict(v => DDG.get_vertex_coords(gx, v) for v in vids)
    c0, c1 = coords(g), coords(gp)

    toP(c::SVector{2,<:Real}) = Point2f(c[1], c[2])

    function segs(cdict)
        s = Point2f[]
        for e in edges(g)  # assume same adjacency in g and gp
            u, v = src(e), dst(e)
            push!(s, toP(cdict[u])); push!(s, toP(cdict[v]))
        end
        return s
    end

    # edges
    linesegments!(ax, segs(c0); linewidth=2)                           # before
    linesegments!(ax, segs(c1); linewidth=2, linestyle=:dash)          # after

    # vertices
    scatter!(ax, [toP(c0[v]) for v in vids if v != hub_vid]; markersize=10)
    scatter!(ax, [toP(c1[v]) for v in vids if v != hub_vid]; marker=:x, markersize=10)

    # hub highlight
    scatter!(ax, [toP(c0[hub_vid])]; marker=:star5, markersize=18)
    scatter!(ax, [toP(c1[hub_vid])]; marker=:star5, markersize=18)

    if show_vids
        for v in vids
            text!(ax, string(v); position=toP(c0[v]), align=(:left, :bottom))
        end
    end
    return ax
end

function run_suite(; k=100.0, r=k/2, ntrials=2000, σ_deg=10.0, seed=1, outprefix="motifs")
    rng = MersenneTwister(seed)
    σ = deg2rad(σ_deg)

    graphs  = Dict{Int,Tuple{Any,Any}}()
    summary = Dict{Int,NamedTuple}()

    fig = Figure(resolution=(1200, 700))

    for d in 1:6
        g = DDG.generate_hub_spoke_graph(k, d)
        gp, bis = extremal_open_arc_push(g, d, r)
        graphs[d] = (g, gp)

        wd0 = DDG.weighted_degree(g, 1)
        wdp = DDG.weighted_degree(gp, 1)

        deltas = Vector{Float64}(undef, ntrials)
        for t in 1:ntrials
            gη, _, _ = extremal_open_arc_push_noisy(g, d, r; σ=σ, rng=rng)
            deltas[t] = DDG.weighted_degree(gη, 1) - wd0
        end

        summary[d] = (
            k=k, r=r, d=d, bisector=bis,
            wd_before=wd0, wd_after=wdp,
            delta_opt=wdp-wd0,
            delta_noise_mean=mean(deltas),
            delta_noise_max=maximum(deltas),
            delta_noise_p99=quantile(deltas, 0.99),
        )

        i, j = fldmod1(d, 3)  # (row, col) with 3 columns
        ax = Axis(fig[i, j];
            aspect=DataAspect(),
            title="d=$d, bis=$(round(rad2deg(bis), digits=1))°",
            xlabel="x", ylabel="y"
        )
        plot_before_after!(ax, g, gp)
    end

    save("$(outprefix)_before_after.png", fig)
    # using JLD2; @save "$(outprefix)_graphs_and_summary.jld2" graphs summary

    return graphs, summary, fig
end

using GLMakie
using LinearAlgebra, StaticArrays, Statistics
const DDG = DynamicGeometricGraphs

# helper: apply your construction with an explicit hub direction θ (radians)
function config_with_hub_angle(g::DynamicGeometricGraph{2,Float64}, d::Int, r::Float64, θ::Float64; hub_vid::Int=1)
    g_new = DDG.freeze(g)

    # move hub
    h  = DDG.get_vertex_coords(g, hub_vid)
    h′ = h + SVector{2,Float64}(r*cos(θ), r*sin(θ))
    DDG.update_coord!(g_new, h, h′)

    # push spokes away from h′
    for vid in sort(collect(DDG.vertices(g)))
        vid == hub_vid && continue
        k  = DDG.get_vertex_coords(g, vid)
        v  = k - h′
        nv = norm(v)
        nv == 0 && continue
        k′ = k + (r/nv) * v
        DDG.update_coord!(g_new, k, k′)
    end
    return g_new
end

"""
High-res GLMakie figure with two panels:
(1) ΔWD(δ) = WD(after) - WD(before) vs δ (degrees) around bisector.
(2) Regret(δ) = ΔWD(0) - ΔWD(δ) (>=0 if bisector is optimal on the scanned range).

Set ngrid odd to ensure δ=0 is included exactly.
"""
function plot_optimality_regret(g, d; r::Float64, hub_vid::Int=1, δmax_deg::Real=180, ngrid::Int=721)
    isodd(ngrid) || error("ngrid must be odd so δ=0 lies on the grid (got $ngrid).")

    # bisector from your construction
    last = -(d - 1) * (π/3)
    occ  = abs(last)
    unocc = 2π - occ
    bis  = last - unocc/2

    wd0 = DDG.weighted_degree(g, hub_vid)

    δgrid = range(-deg2rad(δmax_deg), deg2rad(δmax_deg), length=ngrid)
    Δ = Vector{Float64}(undef, ngrid)

    for (i, δ) in enumerate(δgrid)
        gp = config_with_hub_angle(g, d, r, bis + δ; hub_vid=hub_vid)
        Δ[i] = DDG.weighted_degree(gp, hub_vid) - wd0
    end

    # ΔWD at δ=0 is the middle point (since ngrid odd and symmetric grid)
    Δ0 = Δ[(ngrid + 1) ÷ 2]
    regret = Δ0 .- Δ

    fig = Figure(resolution=(1600, 700))  # high-res; you will save via CairoMakie

    ax1 = Axis(fig[1, 1];
        xlabel="δ from bisector (degrees)",
        ylabel="Δ weighted degree",
        title="ΔWD(δ) (d=$d)"
    )
    lines!(ax1, rad2deg.(δgrid), Δ; linewidth=3)
    vlines!(ax1, [0.0]; linestyle=:dash)

    ax2 = Axis(fig[1, 2];
        xlabel="δ from bisector (degrees)",
        ylabel="regret = ΔWD(0) − ΔWD(δ)",
        title="Regret curve (should be ≥ 0 if bisector is optimal on scan)"
    )
    lines!(ax2, rad2deg.(δgrid), regret; linewidth=3)
    vlines!(ax2, [0.0]; linestyle=:dash)
    hlines!(ax2, [0.0]; linestyle=:dot)

    return fig, (bisector=bis, Δ0=Δ0)
end



# Trial code_lowered
using Random, LinearAlgebra, StaticArrays
using GLMakie
const DDG = DynamicGeometricGraphs

# Rotate a 2D vector by angle ε
@inline function rot2(v::SVector{2,Float64}, ε::Float64)
    c, s = cos(ε), sin(ε)
    return SVector{2,Float64}(c*v[1] - s*v[2], s*v[1] + c*v[2])
end

# Build config given hub direction θ and per-spoke angular noises δspokes (length d)
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
        dir = rot2(dir, δspokes[spoke_idx])       # allow spoke to deviate by δ_i
        k′  = k + r * dir
        DDG.update_coord!(g_new, k, k′)
    end

    return g_new
end

# Bisector from your construction (radians)
function bisector_angle(d::Int)
    last = -(d - 1) * (π/3)
    occ  = abs(last)
    unocc = 2π - occ
    return last - unocc/2
end

"""
Compute envelope max ΔWD over spoke perturbations for each hub deviation in a grid.
Returns (δ0_grid_deg, bestΔWD_per_grid).
"""
function envelope_optimality(g, d; r::Float64, hub_vid::Int=1, δmax_deg::Real=30,
                             ngrid::Int=121, N::Int=50, seed::Int=1)

    rng = MersenneTwister(seed)
    bis = bisector_angle(d)
    wd0 = DDG.weighted_degree(g, hub_vid)

    δ0_grid = range(-deg2rad(δmax_deg), deg2rad(δmax_deg), length=ngrid)
    best = fill(-Inf, ngrid)

    for (j, δ0) in enumerate(δ0_grid)
        θ = bis + δ0
        best_j = -Inf

        for t in 1:N
            δsp = (rand(rng, d) .* 2 .- 1) .* deg2rad(δmax_deg)  # Unif[-δmax, +δmax] per spoke
            gp  = config_hub_and_spoke_noise(g, d, r, θ, collect(δsp); hub_vid=hub_vid)
            Δ   = DDG.weighted_degree(gp, hub_vid) - wd0
            best_j = max(best_j, Δ)
        end
        best[j] = best_j
    end

    return rad2deg.(collect(δ0_grid)), best
end

# Plot the envelope (peak at 0 is your message)
function plot_envelope_GL(g, d; r::Float64, δmax_deg::Real=30, ngrid::Int=121, N::Int=50, seed::Int=1)
    xdeg, env = envelope_optimality(g, d; r=r, δmax_deg=δmax_deg, ngrid=ngrid, N=N, seed=seed)

    fig = Figure(resolution=(1400, 500))
    ax = Axis(fig[1,1];
        xlabel = "hub deviation δ₀ from bisector (degrees)",
        ylabel = "max ΔWD found over spoke perturbations",
        title  = "Local-box envelope test (d=$d, spokes also perturbed)"
    )
    lines!(ax, xdeg, env; linewidth=3)
    vlines!(ax, [0.0]; linestyle=:dash)
    return fig
end


k = 100.0
r = k/2
d = 6
g = DDG.generate_hub_spoke_graph(k, d)

fig = plot_envelope_GL(g, d; r=r, δmax_deg=30, ngrid=121, N=50, seed=1)
display(fig)



### Final 

function envelope_optimality(g, d; r::Float64, hub_vid::Int=1, δmax_deg::Real=30,
                             ngrid::Int=121, N::Int=50, seed::Int=1)

    rng = MersenneTwister(seed)
    bis = bisector_angle(d)
    wd0 = DDG.weighted_degree(g, hub_vid)

    # Pre-sample spoke perturbations ONCE (common random numbers)
    δmax = deg2rad(δmax_deg)
    spoke_samples = [(rand(rng, d) .* 2 .- 1) .* δmax for _ in 1:N]

    δ0_grid = range(-deg2rad(δmax_deg), deg2rad(δmax_deg), length=ngrid)
    best = fill(-Inf, ngrid)
    det  = fill(-Inf, ngrid)  # deterministic baseline (all spoke noises = 0)

    δzero = zeros(Float64, d)

    for (j, δ0) in enumerate(δ0_grid)
        θ = bis + δ0

        # deterministic candidate (true spoke-optimal for that θ)
        gp0 = config_hub_and_spoke_noise(g, d, r, θ, δzero; hub_vid=hub_vid)
        Δ0  = DDG.weighted_degree(gp0, hub_vid) - wd0
        det[j] = Δ0

        # envelope over random spoke perturbations, BUT at least as good as deterministic
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


# Plot the envelope (peak at 0 is your message)
function plot_envelope_GL(g, d; r::Float64, δmax_deg::Real=30, ngrid::Int=121, N::Int=50, seed::Int=1)
    xdeg, env = envelope_optimality(g, d; r=r, δmax_deg=δmax_deg, ngrid=ngrid, N=N, seed=seed)

    fig = Figure(resolution=(1400, 500))
    ax = Axis(fig[1,1];
        xlabel = "hub deviation δ₀ from bisector (degrees)",
        ylabel = "max ΔWD found over spoke perturbations",
        title  = "Local-box envelope test (d=$d, spokes also perturbed)"
    )
    lines!(ax, xdeg, env; linewidth=3)
    vlines!(ax, [0.0]; linestyle=:dash)
    return fig
end


k = 100.0
r = k/2
d = 6
g = DDG.generate_hub_spoke_graph(k, d)

fig = plot_envelope_GL(g, d; r=r, δmax_deg=60, ngrid=121, N=50, seed=1)
display(fig)



### Final 3x2 figure

using GLMakie

# 1) compute-only (no figure creation)
# returns x (degrees), env (best), and optionally det if you have it
function envelope_data(g, d; r, δmax_deg=60, ngrid=121, N=50, seed=1)
    xdeg, det, env = envelope_optimality(g, d; r=r, δmax_deg=δmax_deg, ngrid=ngrid, N=N, seed=seed)
    return xdeg, det, env
end

# 2) plot into an existing axis
function plot_envelope!(ax, xdeg, det, env; title="")
    lines!(ax, xdeg, env; linewidth=3)
    lines!(ax, xdeg, det; linewidth=2, linestyle=:dash)   # optional but recommended
    vlines!(ax, [0.0]; linestyle=:dash)
    ax.title = title
    ax.xlabel = "hub deviation δ₀ (deg)"
    ax.ylabel = "max ΔWD"
    return ax
end

function plot_envelopes_grid(; k=100.0, r=k/2, δmax_deg=60, ngrid=121, N=50, seed=1)
    fig = Figure(resolution=(1400, 900), fontsize=20)

    for d in 1:6
        g = DDG.generate_hub_spoke_graph(k, d)
        xdeg, det, env = envelope_optimality(g, d; r=r, δmax_deg=δmax_deg, ngrid=ngrid, N=N, seed=seed)

        row = (d - 1) % 3 + 1   # 1..3
        col = (d - 1) ÷ 3 + 1   # 1..2

        ax = Axis(fig[row, col]; title="d=$d")
        # ax = Axis(fig[row, col]; title="d=$d", aspect=DataAspect())  # optional alternative

        lines!(ax, xdeg, env; linewidth=3)
        lines!(ax, xdeg, det; linewidth=2, linestyle=:dash)  # optional but recommended
        vlines!(ax, [0.0]; linestyle=:dash)
        ax.xlabel = "hub deviation δ₀ (deg)"
        ax.ylabel = "max ΔWD"
    end

    return fig
end


fig = plot_envelopes_grid(; k=100.0, δmax_deg=60, ngrid=121, N=50, seed=1)
display(fig)
using CairoMakie

save("sweep.pdf", fig)

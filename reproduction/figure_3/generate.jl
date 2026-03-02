#!/usr/bin/env julia

"""
Figure 3: Canonical Geometric Motifs Under Vertex Noise

Generates a 2x3 panel figure showing hub-spoke motifs of degree d=1..6,
each in their maximally sensitive (extremal) geometric configuration.

For each degree d:
  - Spokes are placed at angles 0°, -60°, ..., -(d-1)×60° (clockwise fill)
  - The hub is moved along the bisector of the *unoccupied* arc by distance r
  - Each spoke endpoint is pushed radially away from the perturbed hub by r

The bisector direction maximises the weighted-degree change at the hub.
These are the witness motifs used in Table 1 and the degree-bound theorems.

Author: Ben Cardoen
Source: scripts/patterns_new.jl
"""

import Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using DynamicGeometricGraphs
using CairoMakie
using LinearAlgebra
using StaticArrays
using Graphs

const DDG = DynamicGeometricGraphs

# Okabe-Ito colorblind-safe palette
const OI_BLUE    = colorant"#0072B2"
const OI_ORANGE  = colorant"#E69F00"
const OI_GREEN   = colorant"#009E73"
const OI_PURPLE  = colorant"#CC79A7"

# ============================================================================
# Core: extremal open-arc push
# ============================================================================

"""
    extremal_open_arc_push(g, d, r; hub_vid=1)

Move hub along the open-arc bisector by distance `r`, then push each spoke
endpoint radially away from the perturbed hub by `r`.

Returns (perturbed_graph, bisector_angle_radians).
"""
function extremal_open_arc_push(g::DynamicGeometricGraph{2,Float64}, d::Int, r::Float64;
                                 hub_vid::Int=1)
    r < 0 && error("r must be non-negative, got $r")

    # Bisector of the unoccupied arc
    last_spoke_angle = -(d - 1) * (π/3)       # -(d-1)×60° in radians
    occupied_arc     = abs(last_spoke_angle)
    unoccupied_arc   = 2π - occupied_arc
    bisector         = last_spoke_angle - (unoccupied_arc / 2)

    g_new = DDG.freeze(g)

    # Move hub h → h'
    h      = DDG.get_vertex_coords(g, hub_vid)
    hprime = h + SVector{2,Float64}(r * cos(bisector), r * sin(bisector))
    DDG.update_coord!(g_new, h, hprime)

    # Push each spoke endpoint k radially away from h'
    for vid in sort(collect(DDG.vertices(g)))
        vid == hub_vid && continue
        k  = DDG.get_vertex_coords(g, vid)
        v  = k - hprime
        nv = norm(v)
        nv == 0 && continue
        kprime = k + (r / nv) * v
        DDG.update_coord!(g_new, k, kprime)
    end

    return g_new, bisector
end

# ============================================================================
# Plotting: before/after overlay on a single Axis
# ============================================================================

function plot_before_after!(ax::Axis, g, gp; hub_vid::Int=1)
    vids = sort(collect(DDG.vertices(g)))
    coords(gx) = Dict(v => DDG.get_vertex_coords(gx, v) for v in vids)
    c0, c1 = coords(g), coords(gp)

    toP(c::SVector{2,<:Real}) = Point2f(c[1], c[2])

    function segs(cdict)
        s = Point2f[]
        for e in edges(g)
            u, v = src(e), dst(e)
            push!(s, toP(cdict[u]))
            push!(s, toP(cdict[v]))
        end
        return s
    end

    # Edges: solid blue (before), dashed orange (after)
    linesegments!(ax, segs(c0); linewidth=2.5, color=OI_BLUE)
    linesegments!(ax, segs(c1); linewidth=2.5, linestyle=:dash, color=OI_ORANGE)

    # Spoke vertices: filled circles (before), crosses (after)
    scatter!(ax, [toP(c0[v]) for v in vids if v != hub_vid];
             markersize=12, color=OI_BLUE)
    scatter!(ax, [toP(c1[v]) for v in vids if v != hub_vid];
             marker=:xcross, markersize=12, color=OI_ORANGE)

    # Hub: star markers
    scatter!(ax, [toP(c0[hub_vid])];
             marker=:star5, markersize=20, color=OI_GREEN)
    scatter!(ax, [toP(c1[hub_vid])];
             marker=:star5, markersize=20, color=OI_PURPLE)

    return ax
end

# ============================================================================
# Main: generate 2×3 panel figure
# ============================================================================

function generate_figure(; k=100.0, r=nothing, lim=200.0)
    r === nothing && (r = k / 2)

    fig = Figure(size=(1400, 900), fontsize=16)

    for d in 1:6
        g  = DDG.generate_hub_spoke_graph(k, d)
        gp, bis = extremal_open_arc_push(g, d, r)

        row, col = fldmod1(d, 3)
        ax = Axis(fig[row, col];
            aspect    = DataAspect(),
            title     = "d=$d, bis=$(round(rad2deg(bis), digits=1))°",
            titlesize = 18,
            xlabel    = "x",
            ylabel    = "y",
        )
        plot_before_after!(ax, g, gp)
        xlims!(ax, -lim, lim)
        ylims!(ax, -lim, lim)

        # Print weighted-degree change for verification (used in Table 1)
        wd0 = DDG.weighted_degree(g, 1)
        wdp = DDG.weighted_degree(gp, 1)
        Δ   = wdp - wd0
        rel = Δ / wd0
        println("  d=$d: WD_before=$(round(wd0,digits=2)), WD_after=$(round(wdp,digits=2)), " *
                "ΔWD=$(round(Δ,digits=2)), rel=$(round(rel,digits=4)), " *
                "bisector=$(round(rad2deg(bis),digits=1))°")
    end

    return fig
end

# ============================================================================
# Run
# ============================================================================

println("=" ^ 70)
println("Generating Figure 3: Canonical Motifs Under Vertex Noise")
println("=" ^ 70)

fig = generate_figure(; k=100.0, lim=200.0)

output_pdf = joinpath(@__DIR__, "Figure3.pdf")
output_png = joinpath(@__DIR__, "Figure3.png")
save(output_pdf, fig)
save(output_png, fig, px_per_unit=2)

println("\nSaved: $output_pdf")
println("Saved: $output_png")
println("=" ^ 70)

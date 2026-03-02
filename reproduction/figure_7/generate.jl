#!/usr/bin/env julia

"""
Figure 7: Graph Constructions for Motif-Frequency and Spanner Experiments

Generates a 1×3 panel figure showing three geometric graph constructions
that share the same spatial scaffold of 6 hub-spoke motifs in a 2×3 grid:

  (A) Heterogeneous motifs (degrees 1–6) — the reference graph `refgraph()`
  (B) ε-thickness violation (hub removed) — hub vertex removed + 3 crossing edges
  (C) Uniform motifs (6 × degree 4) — all motifs have the same degree

Author: Ben Cardoen
Source: archive/fig7_motif_panel/generate.jl
"""

import Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using DynamicGeometricGraphs
using CairoMakie
using Graphs
using StaticArrays

const DDG = DynamicGeometricGraphs

# ============================================================================
# Violated graph: remove hub, add crossing edges
# ============================================================================

function make_violated(g)
    gp = freeze(g)
    Graphs.rem_vertex!(gp, 1)
    Graphs.add_edge!(gp, 24, 9)
    Graphs.add_edge!(gp, 13, 3)
    Graphs.add_edge!(gp, 6, 21)
    return gp
end

# ============================================================================
# 2D graph plotter
# ============================================================================

function plot_graph_2d!(ax, g; nodesize=8, linewidth=1.5,
                        nodecolor=:black, edgecolor=:black)
    vids = collect(Graphs.vertices(g))
    v2idx = Dict{Int, Int}()
    xs = Float64[]
    ys = Float64[]

    for (i, vi) in enumerate(vids)
        v2idx[vi] = i
        coords = DDG.get_vertex_coords(g, vi)
        push!(xs, coords[1])
        push!(ys, coords[2])
    end

    for e in Graphs.edges(g)
        u, v = src(e), dst(e)
        i1, i2 = v2idx[u], v2idx[v]
        lines!(ax, [xs[i1], xs[i2]], [ys[i1], ys[i2]];
               color=edgecolor, linewidth=linewidth)
    end

    scatter!(ax, xs, ys; color=nodecolor, markersize=nodesize)
end

# ============================================================================
# Main: generate 1×3 panel figure
# ============================================================================

function generate_figure()
    G_normal  = refgraph()
    G_uniform = refgraph(; pattern=[4,4,4,4,4,4])
    G_violated = make_violated(G_normal)

    println("  Normal:   nv=$(nv(G_normal)), ne=$(ne(G_normal))")
    println("  Violated: nv=$(nv(G_violated)), ne=$(ne(G_violated))")
    println("  Uniform:  nv=$(nv(G_uniform)), ne=$(ne(G_uniform))")

    fig = Figure(size=(1200, 400), fontsize=12)

    titles = ["(A) Heterogeneous motifs (degrees 1–6)",
              "(B) ε-thickness violation (hub removed)",
              "(C) Uniform motifs (6 × degree 4)"]
    graphs = [G_normal, G_violated, G_uniform]

    for (i, (g, title)) in enumerate(zip(graphs, titles))
        ax = Axis(fig[1, i];
            title=title,
            xlabel="x",
            ylabel= i == 1 ? "y" : "",
            aspect=DataAspect(),
        )
        plot_graph_2d!(ax, g; nodesize=6, linewidth=1.2)
        if i > 1
            hideydecorations!(ax; grid=false)
        end
    end

    colgap!(fig.layout, 10)
    return fig
end

# ============================================================================
# Run
# ============================================================================

println("=" ^ 70)
println("Generating Figure 7: Graph Constructions (Motif Panel)")
println("=" ^ 70)

fig = generate_figure()

output_pdf = joinpath(@__DIR__, "Figure7.pdf")
output_png = joinpath(@__DIR__, "Figure7.png")
save(output_pdf, fig)
save(output_png, fig, px_per_unit=3)

println("\nSaved: $output_pdf")
println("Saved: $output_png")
println("=" ^ 70)

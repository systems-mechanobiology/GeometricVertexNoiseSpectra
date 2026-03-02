# Figure 6: Motif Geometry Sensitivity (Optimality Envelopes)

## Paper reference
- Label: `fig-motifs-optimality`
- File: `figures/Figure6.pdf`
- Paper figure number: 6

## Description
3×2 panel figure showing, for each motif degree d=1..6, the maximum weighted-degree change ΔWD as a function of hub deviation angle δ₀ from the bisector of the unoccupied arc.

Two curves per panel:
- **Deterministic baseline (dashed):** hub moved at angle bisector+δ₀, spokes pushed radially with zero angular noise.
- **Envelope (solid):** maximum ΔWD over N=50 random spoke-angle perturbations using common random numbers (fixed seed).

## Reproduction

### Dependencies
- Julia (tested with 1.10+)
- Shared environment: `reproduction/Project.toml`
  - `DynamicGeometricGraphs.jl`
  - `CairoMakie`, `LinearAlgebra`, `StaticArrays`, `Random`

### Generate
```bash
julia reproduction/figure_6/generate.jl
```

### Output
- `Figure6.pdf` — publication-quality PDF (matches `figures/Figure6.pdf`)
- `Figure6.png` — quick preview

## Provenance
- **Script origin:** `scripts/patterns_new.jl` (functions extracted verbatim)
- **Author:** Ben Cardoen

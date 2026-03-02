# Figure 7: Graph Constructions (Motif Panel)

## Paper reference
- Label: `fig-motifgraph`
- File: `figures/Figure7.pdf`
- Paper figure number: 7

## Description
1×3 panel figure showing three geometric graph constructions sharing the same spatial scaffold of six hub-spoke motifs arranged in a 2×3 grid:

- **(A) Heterogeneous motifs (degrees 1–6):** reference graph used in motif-frequency experiments
- **(B) ε-thickness violation:** hub removed and crossing edges added
- **(C) Uniform motifs (6 × degree 4):** frequency-controlled uniform construction

## Reproduction

### Dependencies
- Julia (tested with 1.10+)
- Shared environment: `reproduction/Project.toml`
  - `DynamicGeometricGraphs.jl`
  - `CairoMakie`, `Graphs`, `StaticArrays`

### Generate
```bash
julia reproduction/figure_7/generate.jl
```

### Output
- `Figure7.pdf` — publication-quality PDF (matches `figures/Figure7.pdf`)
- `Figure7.png` — quick preview

## Provenance
- **Script origin:** `archive/fig7_motif_panel/generate.jl`
- **Author:** Ben Cardoen

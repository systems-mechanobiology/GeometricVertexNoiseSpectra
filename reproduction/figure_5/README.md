# Figure 5: Motif-informed graph-level sensitivity and spanner effects (S3I)

## Paper reference
- Label: `fig-s3-results`
- File: `figures/new/Figure5.pdf`
- Paper figure number: 5

## Description
Evaluates S3I (stochastic spectral separability index) on controlled geometric graphs
assembled from witness motifs, testing four sources of spectral variation:

- **Panel A (top-left):** Noise amplitude — S3I between c₁=1 and c₂=[2,3,4,5], same graph
- **Panel B (top-right):** Noise radius — S3I between L₁=10 and L₂=[20,30,40,50], same graph
- **Panel C (bottom-left):** Motif frequency — heterogeneous motifs (1,2,3,4,5,6) vs uniform (4×6)
- **Panel D (bottom-right):** ε-thick vs unconstrained graph construction

Each panel shows jittered scatter (100 runs), bootstrap CI lines, and green ✓ for significance.

## Reproduction

### Dependencies
- Julia (tested with 1.10+)
- Shared environment: `reproduction/Project.toml`
  - `CalibratedTemperedPowerLaw.jl` (for `refgraph`, `strong_distances`, `SC`, `S3I_from_SC`, `bootstrap_s3i_test`, `autotune_lambda`)
  - `DynamicGeometricGraphs.jl` (for graph construction)
  - `CairoMakie`, `Statistics`, `Graphs`, `LinearAlgebra`, `StaticArrays`

### Generate
```bash
julia reproduction/figure_5/generate.jl
```

**Note:** This is computationally heavy (~100 runs × 200 samples × 4 panels).

### Output
- `Figure5.pdf` — publication-quality PDF
- `Figure5.png` — quick preview
- `s3i_panel_a_vary_c.csv` — raw S3I values, Panel A (400 rows)
- `s3i_panel_b_vary_L.csv` — raw S3I values, Panel B (400 rows)
- `s3i_panel_c_motif_freq.csv` — raw S3I values, Panel C (500 rows)
- `s3i_panel_d_spanner.csv` — raw S3I values, Panel D (500 rows)

### Parameters
```julia
n_samples = 200   # Spectral samples per S3I estimate
n_runs = 100      # Bootstrap runs
c1 = 1            # Base noise amplitude
limit = 40        # Base noise radius
p = 0.05          # Tail probability for λ auto-tuning
seed = 42         # Base random seed
```

### Verification
All 18 comparisons show ✓ (statistically significant separation at α=0.05).
Distribution shapes and CI ranges match the original figure.

## Provenance
- **Script origin:** `archive/fig10_s3i_tests/generate.jl`
- **Author:** Ben Cardoen

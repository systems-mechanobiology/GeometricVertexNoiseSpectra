# Figure 4: Vertex Noise Induces Correlated Edge Perturbations

## Paper reference
- Label: `fig-vertex-noise-is-not-edge-noise`
- File: `figures/new/Figure4.pdf`
- Paper figure number: 4

## Description
Demonstrates that i.i.d. vertex noise induces non-i.i.d. (correlated) edge noise,
even in a minimal 3-spoke hub motif.

- **Row 1, Panel A:** Strong oracle — true vertex positions with uncertainty circles (r=50)
- **Row 1, Panel B:** Weak oracle — 30 noisy realisations overlaid (semi-transparent)
- **Row 2:** Correlation matrices (3×3, E1/E2/E3) varying c = 1, 3, 5 with ε=0, k=2r
- **Row 3:** Correlation matrices varying ε = 0, 5, 10 with c=1 fixed

Key finding: off-diagonal correlations are non-zero and structured (positive and negative),
confirming vertex noise cannot be modeled as independent edge noise.

## Reproduction

### Dependencies
- Julia (tested with 1.10+)
- Shared environment: `reproduction/Project.toml`
  - `CalibratedTemperedPowerLaw.jl` (for `autotune_lambda`, `levy_perturb`)
  - `CairoMakie`, `Statistics`, `LinearAlgebra`, `StaticArrays`

### Generate
```bash
julia reproduction/figure_4/generate.jl
```

### Output
- `Figure4.pdf` — publication-quality PDF
- `Figure4.png` — quick preview
- `correlations.csv` — all 54 correlation values (6 matrices × 9 entries) for audit

### Verification
All 18 off-diagonal correlations match the original figure at 2 d.p. precision.
Seeds are identical to original script, producing reproducible results.

### Parameters
```julia
k = 100.0       # Edge length
r = 50.0        # Noise radius (k = 2r)
n_samples = 500 # Samples per correlation estimate
δ = 0.05        # Tail probability for λ auto-tuning
seed = 42       # Base random seed
```

## Notes
- 3-row layout (was 2×4): graph panels get more space, heatmaps more readable
- Okabe-Ito colorblind-safe palette
- CSV audit trail added (not in original)

## Provenance
- **Script origin:** `scripts/graph_noise_model_figure.jl`
- **Author:** Ben Cardoen

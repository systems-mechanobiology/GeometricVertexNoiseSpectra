# Figure 3: Canonical Geometric Motifs Under Vertex Noise

## Paper reference
- Label: `fig-motifs-worst-case`
- File: `figures/new/Figure3.pdf`
- Paper figure number: 3

## Description
Six panels (2x3 grid) showing hub-spoke motifs of degree d=1 through d=6
in their maximally sensitive (extremal) geometric configuration:
- **Blue solid** lines/dots: noise-free vertex positions and edges
- **Orange dashed** lines/crosses: perturbed positions and edges
- **Green star**: hub vertex (noise-free)
- **Purple star**: hub vertex (perturbed)

For each degree d, spokes are placed at 60-degree intervals (clockwise).
The hub is displaced along the bisector of the unoccupied arc — the direction
that maximises weighted-degree change. Spoke endpoints are then pushed
radially away from the perturbed hub.

## Critical importance
These extremal configurations define the witness motifs used in:
- **Table 1**: maximum relative weighted-degree change per degree
- **Theorem 4.1 / 4.2**: degree-dependent spectral perturbation bounds

If the bisector calculation or push direction is wrong, those results are invalid.

## Reproduction

### Dependencies
- Julia (tested with 1.10+)
- Shared environment: `reproduction/Project.toml`
  - `DynamicGeometricGraphs.jl`
  - `CairoMakie`, `LinearAlgebra`, `StaticArrays`, `Graphs`

### Generate
```bash
julia reproduction/figure_3/generate.jl
```

### Output
- `Figure3.pdf` — publication-quality PDF (matches `figures/new/Figure3.pdf`)
- `Figure3.png` — quick preview
- Console output: weighted-degree changes per degree (cross-check with Table 1)

### Verification
The script prints WD_before, WD_after, ΔWD, and relative change for each degree.
These values must match Table 1 in the paper.

## Provenance
- **Script origin:** `scripts/patterns_new.jl` (functions: `extremal_open_arc_push`, `plot_before_after!`, `run_suite`)
- **Author:** Ben Cardoen

## Notes
- Colors updated to Okabe-Ito colorblind-safe palette (consistent with Figure 2)
- The file is named Figure3.pdf because that is what the paper references (`figures/new/Figure3.pdf`)

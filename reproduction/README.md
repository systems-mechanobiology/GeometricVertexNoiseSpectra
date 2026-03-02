# Reproduction

Self-contained scripts to reproduce every figure and table in the paper.
Each subfolder contains a `generate.jl` script and a README documenting the algorithm, parameters, and verification steps.

## Figures and Tables

| Folder | Paper Element | Description |
|--------|--------------|-------------|
| `figure_1/` | Figure 1 | Graph model schematic (Inkscape, not generated) |
| `figure_2/` | Figure 2 | CTPL noise model (PDF, tail dominance, auto-tuned λ) |
| `figure_3/` | Figure 3 | Canonical motifs under vertex noise (d=1–6) |
| `figure_4/` | Figure 4 | Vertex noise induces correlated edge perturbations |
| `figure_5/` | Figure 5 | S3I graph-level sensitivity |
| `figure_6/` | Figure 6 | Motif angular sensitivity (optimality envelopes) |
| `figure_7/` | Figure 7 | Graph constructions (heterogeneous, violated, uniform) |
| `table_1/`  | Table 1  | Motif spectral tightness |

## Running

Each script activates the shared Julia environment (`Project.toml` / `Manifest.toml` in this directory).
Install dependencies once, then run any script:

```bash
julia --project=reproduction -e 'using Pkg; Pkg.instantiate()'
julia reproduction/figure_X/generate.jl
```

**Note:** The manifest uses local (path) checkouts of:
- `DynamicGeometricGraphs.jl`
- `CalibratedTemperedPowerLaw.jl`

These are expected to be cloned alongside this repository (see the paper Appendix / repo-level README for links).
If your checkout layout differs, update the corresponding `path = ...` entries in `reproduction/Manifest.toml` or `Pkg.develop` the packages from their locations.

## Resources

`resources/` contains reference implementations from which per-figure scripts were derived.
See `resources/README.md` for details.

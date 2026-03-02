# Table 1: Motif Spectral Tightness

## Paper reference
- Label: `tbl-degree-summary`
- Paper table number: 1

## Description
Maximum relative change in weighted degree for maximally affected subgraphs
under CTPL noise, together with observed spectral displacement and Weyl bounds.

Columns:
- **Degree (d):** Hub-spoke motif degree (1-6)
- **δWD ε→0:** Weighted-degree change coefficient in the tight case (k = 2r)
- **Actual Shift:** Maximum eigenvalue shift of the graph Laplacian
- **Weyl Limit:** Weyl bound ||E||₂ on eigenvalue perturbation
- **Ratio:** Actual Shift / Weyl Limit (1.0 = Weyl-tight)

## Critical importance
This table depends directly on the extremal motif configurations shown in Figure 3.
The geometric construction (bisector direction, spoke push rule) must match exactly.
Cross-checked against Figure 3's `extremal_open_arc_push`:
- Bisector angles match (mod 360° equivalence verified)
- Relative WD changes match to 4 decimal places
- All table values match the paper at printed precision

## Key results
- Degrees 1 and 2 are **Weyl-tight** (ratio = 1.0000)
- d=6 δWD ≈ d=5 δWD because the 60° unoccupied arc constrains hub movement

## Reproduction

### Dependencies
- Julia (tested with 1.10+)
- `LinearAlgebra` (standard library only)

### Generate
```bash
julia reproduction/table_1/generate.jl
```

### Output
- Console: full table with geometry verification details
- `table1.md`: markdown-formatted table for comparison with paper

### Parameters
```julia
r = 1.0      # Noise magnitude
k = 2.0 * r  # Edge length (tight case: ε → 0)
```

## Provenance
- **Script origin:** `archive/table1_motif_spectral/generate.jl`
- **Author:** Ben Cardoen

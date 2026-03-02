# Resources

Reference implementations used as source material for the per-figure reproduction scripts.

## patterns_extremal_motifs.jl

Canonical reference implementation for hub-spoke extremal motif geometry.
Contains the core functions extracted (verbatim) into the Figure 3, 6, and Table 1 reproduction scripts:

- `extremal_open_arc_push()` — move hub along unoccupied-arc bisector, push spokes radially
- `config_with_hub_angle()` — generalised hub direction for angular sweeps
- `config_hub_and_spoke_noise()` — hub + per-spoke angular perturbation
- `bisector_angle()` — compute bisector of unoccupied arc for degree d
- `envelope_optimality()` — sweep hub deviation with spoke perturbation envelope
- `plot_envelopes_grid()` — 3×2 panel production figure (Figure 6)

Parameters: k=100 (spoke length), r=k/2 (perturbation radius), degrees d=1–6.

## animated_demo.jl

Generates the CTPL noise demonstration GIF (`ctpl_demo.gif`) shown in the repository README.
Requires GLMakie, CalibratedTemperedPowerLaw.jl, and DynamicGeometricGraphs.jl.

- Left panel: reference graph (red) with accumulated noisy overlays (blue/cyan, alpha-blended)
- Right panel: SC strong (oracle) and weak (observable) convergence with variance bands

Parameters: 100 samples, k=100, 30 vertices, c=5, limit=40, seed=42.

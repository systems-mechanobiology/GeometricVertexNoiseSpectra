# Figure 2: CTPL Noise Model

## Paper reference
- Label: `fig-tl`
- File: `figures/new/Figure2.pdf`

## Description
Three-panel figure illustrating the calibrated tempered power-law (CTPL) noise model:
- **Panel A (PDF Comparison):** PDFs showing the effect of exponential tempering lambda on a Levy distribution, preserving a power-law regime while suppressing extreme displacements beyond a truncation scale.
- **Panel B (Tail Dominance):** Survival functions demonstrating that the calibrated CTPL distribution retains tail dominance relative to moment-matched Gaussian and Gamma alternatives under a fixed tail probability constraint.
- **Panel C (Auto-Tuned lambda):** Empirical auto-tuning of lambda as a function of scale parameter c, obtained via rejection sampling to enforce P(X > r_max) <= delta.

## Reproduction

### Dependencies
- Julia (tested with 1.10+)
- [CalibratedTemperedPowerLaw.jl](https://github.com/systems-mechanobiology/CalibratedTemperedPowerLaw.jl) — provides `autotune_lambda`, `sample_tempered_levy`, `levy_pdf`, `tempered_levy_pdf`
- Julia packages: `CairoMakie`, `Statistics`, `Distributions`, `Random`, `Printf`

### Generate
```bash
julia reproduction/figure_2/generate.jl
```

### Output
- `Figure2.pdf` — publication-quality PDF (matches `figures/new/Figure2.pdf`)
- `Figure2.png` — quick preview

## Notes
- 2-row layout: Panel A full width on top, Panels B+C side by side below
- Okabe-Ito colorblind-safe palette (no red/green)
- Legend inside top-right of Panel A to maximize plot area

## Provenance
- **Script origin:** `archive/fig2_ctpl_distribution/generate.jl`
- **Author:** Ben Cardoen

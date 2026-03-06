# Geometric Vertex and Topological Noise Spectra for Dynamic Geometric Graphs.

[![Render Paper](https://github.com/systems-mechanobiology/GeometricVertexNoiseSpectra/actions/workflows/render-paper.yml/badge.svg)](https://github.com/systems-mechanobiology/GeometricVertexNoiseSpectra/actions/workflows/render-paper.yml)

This repository holds all reproduction content for our paper on vertex noise in geometric graphs in presence of tempered power law noise model and the corresponding limits and approximations of the spectral impact.

![CTPL Demo Animation](ctpl_demo.gif)

*Animated demonstration of CTPL noise effects on geometric graphs. Left: accumulated noisy graph overlays showing noise distribution. Right: SC strong (oracle) and weak (observable) convergence with variance bands.*

## Building the paper

Note, this is a quarto-powered paper, you will need quarto installed, see dependencies.

### PDF (clean build)

```bash
./build_paper.sh
```

Cleans stale LaTeX auxiliary files and renders to PDF. Output: `paper.pdf`. Works on MacOS / Linux, not on Windows. For Windows, use explicit removal of stale files + quarto render.

### HTML

```bash
quarto render paper.qmd --to html
```

### Live preview

```bash
quarto preview paper.qmd
```

Opens a local server with hot-reload — the browser refreshes automatically on each save to `paper.qmd` or any referenced figure.

## Reproducing figures

Each figure and table has a self-contained generation script in `reproduction/`.
Scripts self-activate the shared environment, but dependencies must be instantiated once:

```bash
julia --project=reproduction -e 'using Pkg; Pkg.instantiate()'
julia reproduction/figure_X/generate.jl
```

See [`reproduction/README.md`](reproduction/README.md) for the full index.
`reproduction/Project.toml` and `reproduction/Manifest.toml` pin the Julia environment; it expects local checkouts of `DynamicGeometricGraphs.jl` and `CalibratedTemperedPowerLaw.jl` (see Dependencies). No external data sources are needed.

## Reviewing

This repository uses GitHub Issues as the structured review channel. Reviewers do not need to commit or open PRs — comments are captured in issues using tagged categories (FIX, CLARIFY, CHANGE, QUESTION, BLOCKER, NIT).

- [`agents/reviewing-guide.md`](agents/reviewing-guide.md) — full reviewer guide and workflow
- [`ISSUE_TEMPLATE/review.md`](ISSUE_TEMPLATE/review.md) — issue template for structured review rounds

## Dependencies

- [Quarto](https://quarto.org/) (paper rendering)
- [Julia 1.10+](https://julialang.org/) (figure reproduction)
- [DynamicGeometricGraphs.jl](https://github.com/systems-mechanobiology/DynamicGeometricGraphs.jl)
- [CalibratedTemperedPowerLaw.jl](https://github.com/systems-mechanobiology/CalibratedTemperedPowerLaw.jl) (CTPL.jl)
- [GitHub CLI](https://cli.github.com/) (`gh`) — optional; needed for the issue-based review workflow

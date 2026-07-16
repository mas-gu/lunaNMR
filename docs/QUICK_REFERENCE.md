# LunaNMR Quick Reference Card

One-page cheat sheet for common operations and parameters.

---

## Launch
```bash
cd lunaNMR_v1o0
python3 launch_lunaNMR.py
```

## Main Workflow
```
Load Spectrum → Detect Peaks → Fit All Peaks → Export Results
```

---

## Detection Parameters (GUI defaults)

| Parameter | Default | Purpose |
|-----------|---------|---------|
| Search window 1H / X (ppm) | 0.070 / 0.070 | Both dims; keep tight for low-S/N (e.g. saturated hetNOE) |
| Noise Threshold | 3.0 | S/N multiplier |
| Centroid window X / Y (ppm) | 0.01 / 0.10 | Max position shift |

## Fitting Controls

| Toggle | Default | Effect |
|--------|---------|--------|
| PS2D multi-peak | ✓ ON | Auto 2D fitting for overlaps |
| Fix Positions | ☐ OFF | Zero drift (use for series) |
| Fix Linewidths | ☐ OFF | Lock widths to reference |
| Parallel Processing | ✓ ON | ~2.7× faster |
| Per-Peak Linewidth Reuse | ☐ OFF | Learned LW as initial guess |

---

## Overlap Thresholds (auto-clustering)

Peaks within these elliptical radii cluster for 2D fitting.

| Nucleus | F1 / Y (ppm) | F2 / 1H (ppm) |
|---------|--------------|---------------|
| 15N | ±0.4 | ±0.04 |
| 13C | ±0.15 | ±0.04 |

Full clustering/fitting algorithm: [ALGORITHMS.md](ALGORITHMS.md).

---

## Quality Indicators (marker colors)

| Color | R² | Category |
|-------|-----|----------|
| 🟩 Lime | ≥ 0.9 | Excellent |
| 🟢 Green | ≥ 0.8 | Good |
| 🟠 Orange | ≥ 0.5 | Fair |
| 🔴 Red | < 0.5 | Poor / Failed / unfitted |

Marker colors follow the result `quality` field (`_assign_quality_from_r2`, `core_integrator.py`). The standalone 1D Voigt fitter uses a stricter internal scale (0.95/0.85/0.70) that does not drive marker color.

---

## Series Integration (Titration/Relaxation)

1. Fit reference spectrum
2. ✓ Fix Positions ON
3. Set as Reference
4. Process Series
5. Export Series Results → CSV

Writes `series_results_*/`: `comprehensive_peak_tracking.csv` (wide), `series_analysis_tidy.csv` (long), `peak_intensity_matrix.csv`. Fitted positions to 4 decimals (CSP needs it).

---

## Analysis Modules (Modules menu)

| Module | Purpose |
|--------|---------|
| **DynamiXs** | T₁/T₂, methyl T₂ (bi-exp), spectral density, model-free, CPMG |
| **Kd / Titration** | Binding Kd — CSP (quadratic isotherm) + intensity decay |

---

## Headless CLI

```bash
python -m lunaNMR <subcommand>          # or: lunanmr <subcommand>
```

| Command | Purpose |
|---------|---------|
| `series` | Process a series/titration → intensity/position matrices |
| `dynamixs t1t2 / methyl-t2 / hetnoe / density / modelfree` | Relaxation + model-free |
| `kd` / `export kd` | Kd titration fit / figures from a fit JSON |
| `project inventory / remove` | Inspect / prune a `.lunaNMR` bundle |
| `batch` | Folder-wide detect + fit |

Full flags: [CLI.md](CLI.md).

---

## Project Bundles

`.lunaNMR` directories save/reopen full state. Fit surfaces are stored slim and
**reconstructed on load** (spectrum must be loaded first). DynamiXs and Kd persist as
multiple named analyses; reopen from the Project Browser (File ▸ Project Contents…).

---

## File Formats

| Format | Extension | Software |
|--------|-----------|----------|
| Bruker | `.2ii`, `.2rr` | TopSpin |
| NMRPipe | `.ft`, `.pipe` | NMRPipe |
| Varian | `.fid`, `.ft` | VnmrJ |
| SPARKY | `.ucsf` | SPARKY |

---

## Performance (150 peaks, 6-core CPU)

| Mode | Fitting |
|------|---------|
| Sequential | ~120 s |
| Parallel | ~45 s (2.7×) |

---

## Code References

| Component | File |
|-----------|------|
| Peak Detection | `core/enhanced_peak_picker.py` |
| Clustering + Orchestrator | `core/core_integrator.py` (`identify_overlap_clusters`) |
| PS2D 2D Fitting | `core/ps2d_2d_fitter.py` |
| 1D Voigt / fallback | `core/enhanced_voigt_fitter.py` |
| Parallel Processing | `core/parallel_voigt_processor.py` |
| Multi-Spectrum Overlay | `gui/dialogs/multi_spectrum_viewer_dialog.py` |

---

## Documentation Index

| Topic | File |
|-------|------|
| Installation | [INSTALLATION.md](INSTALLATION.md) |
| 5-min tutorial | [README.md](README.md#getting-started-5-minutes) |
| GUI features | [lunaNMR_guide.html](lunaNMR_guide.html) |
| Series integration | [MULTI_SPECTRUM_GUIDE.md](MULTI_SPECTRUM_GUIDE.md) |
| Headless CLI | [CLI.md](CLI.md) |
| Algorithms & fitting | [ALGORITHMS.md](ALGORITHMS.md) |
| Architecture | [ARCHITECTURE.md](ARCHITECTURE.md) |

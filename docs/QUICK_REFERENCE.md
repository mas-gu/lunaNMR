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
Load Data → Detect → Fit Spectrum → File ▸ Export Results
```

---

## Detection Parameters (GUI defaults)

| Parameter | Default | Purpose |
|-----------|---------|---------|
| Search window 1H / X (ppm) | 0.070 / 0.070 | Both dims, 15N and 13C (`ps2d_config.SEARCH_WINDOW_PPM`); keep tight for low-S/N (e.g. saturated hetNOE) |
| Centroid window X / Y (ppm) | 0.01 / 0.10 | Max position shift |

## Fitting Controls

| Toggle | Default | Effect |
|--------|---------|--------|
| PS2D multi-peak | always on (no toggle) | Auto 2D fitting for overlaps |
| Fix Positions | ☐ OFF | Zero drift (use for series) |
| Fix Linewidths | ☐ OFF | Lock widths to reference |
| Use Parallel Processing | ✓ ON | ~2.7× faster on large peak lists; same results as sequential |
| Adaptive Parameter Optimization | ☐ OFF | radF1/radF2 grid search; Fit Spectrum only — a series never runs it |
| Use Spectrum 1 Linewidths as Initial Guess | ☐ OFF | Learned LW as initial guess |
| Auto-add nearby unmatched peaks | ✓ ON | Adds detected peaks missing from the peak list (dummies), in a series too |

---

## Overlap Thresholds (auto-clustering)

Each peak carries an ellipse with these semi-axes; two peaks cluster for 2D fitting when their ellipses touch (normalised elliptical distance ≤ 2, i.e. < 0.08 ppm ¹H / < 0.8 ppm ¹⁵N on axis).

| Nucleus | F1 / Y (ppm) | F2 / 1H (ppm) |
|---------|--------------|---------------|
| 15N | ±0.4 | ±0.04 |
| 13C | ±0.15 | ±0.04 |

Full clustering/fitting algorithm: [ALGORITHMS.md](ALGORITHMS.md).

---

## Quality Indicators (marker colors)

| Color | R² | Category |
|-------|-----|----------|
| 🟩 Lime | > 0.9 | Excellent |
| 🟢 Green | > 0.8 | Good |
| 🟠 Orange | > 0.6 (0.5–0.6 at S/N ≤ 5) | Fair |
| 🔴 Red | below Fair | Poor / Failed / unfitted |

Failed is R² < 0.2 at S/N ≤ 5, or a fit that did not succeed; a near-noise peak is Poor, never Failed. Marker colors follow the label the fit carries, given by the one rule every fitter and the series use: `fit_quality_label` (`utils/quality_categories.py`); full table in [ALGORITHMS.md §4](ALGORITHMS.md).

---

## Series Integration (Titration/Relaxation)

1. Fit the reference spectrum (Fit Spectrum)
2. **Fit Series** → choose the Peak Source (Cascade / Detected / Independent; Independent is refused for a titration) and the X-axis mode (Off / Time series / Titration point)
3. ▶ Start Integration (✕ Cancel stops it within a cluster; a cancelled series writes nothing)

Fix Positions holds every peak at its seed: leave it off for a titration (CSP).

Writes `series_results_*/`: `comprehensive_peak_tracking.csv` (wide), `series_analysis_tidy.csv` (long), `peak_intensity_matrix.csv` / `peak_volume_matrix.csv` / `peak_detection_matrix.csv` / `peak_intensity_detected_matrix.csv`, `series_metadata.json`. Fitted positions to 4 decimals (CSP needs it).

---

## Analysis Modules (Modules menu)

| Module | Purpose |
|--------|---------|
| **DynamiXs** | T₁/T₂, methyl T₂ (bi-exp), T₁ρ, spectral density, model-free |
| **Kd / Titration** | Binding Kd — CSP (quadratic isotherm) + intensity decay |
| **1D Integration** | 1D peak picking, then height / region sum / Voigt area across a series |

---

## Headless CLI

```bash
python -m lunaNMR <subcommand>          # `lunanmr` needs `pip install -e .`
```

| Command | Purpose |
|---------|---------|
| `diagnose` | Read-only pre-flight over a dataset before running anything |
| `series` | Process a series/titration → intensity/position matrices |
| `dynamixs t1t2 / methyl-t2 / hetnoe / t1rho / density / modelfree` | Relaxation + model-free |
| `kd` / `export kd` | Kd titration fit / figures from a fit JSON |
| `integrate-1d` | 1D peak integration across a series |
| `peaks shift` | Apply the registration offset `diagnose` measured |
| `project build / inventory / export / remove` | Build a `.lunaNMR` bundle from a run; inspect / copy out / prune one |
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
| 5-min tutorial | [README.md](README.md#getting-started) |
| GUI features | [lunaNMR_guide.html](lunaNMR_guide.html) |
| Series integration | [MULTI_SPECTRUM_GUIDE.md](MULTI_SPECTRUM_GUIDE.md) |
| Headless CLI | [CLI.md](CLI.md) |
| Algorithms & fitting | [ALGORITHMS.md](ALGORITHMS.md) |
| Architecture | [ARCHITECTURE.md](ARCHITECTURE.md) |

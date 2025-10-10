# Algorithm Reference – lunaNMR v0.9

This document summarises the major numerical routines used by lunaNMR: peak detection, PS2D simultaneous Voigt fitting, overlap resolution, consensus fitting, and quality assessment.

## 1. Peak Detection Pipeline

1. **Pre-processing** (`core/enhanced_peak_picker.py`)
   - Noise estimation via robust statistics (corner sampling or MAD).
   - Baseline correction options (ArPLS, polynomial, rolling median).
2. **Candidate search**
   - Prominence-based detection combined with derivative checks.
   - Graph clustering (`networkx`) merges nearby candidates.
3. **Context metrics** recorded for QA and ML:
   - S/N, local baseline slope, neighbouring peak count, isolation scores.
4. **Output schema**
   - Each candidate carries positions (`ppm_x`, `ppm_y`), initial linewidths, detection confidence, and is stored for downstream fits.

## 2. PS2D Simultaneous Voigt Fitting

Implemented in `core/ps2d_2d_fitter.py` and powered by the Faddeeva-based primitives from `ps2d_style_fitter.py`.

### 2.1 Data Selection

Before fitting, `ps2d_data_selector.py` builds a union of elliptical masks around an overlap group:

```
radius² = (ΔF1 / radF1)² + (ΔF2 / radF2)² ≤ 1
```

Typical defaults: `radF1 = 0.4 ppm (15N/13C)`, `radF2 = 0.04 ppm (1H)`.

### 2.2 Parameter Layout

For `n` peaks the solver receives an array of length `8 × n`:

```
[pos_f1, lw_lor_f1, lw_gau_f1,
 pos_f2, lw_lor_f2, lw_gau_f2,
 intensity, spare]
```

Widths use FWHM units; Gaussian widths are converted to `σ = FWHM / √(8 ln 2)` internally, Lorentzian to `γ = FWHM / 2`.

### 2.3 Five-Stage Levenberg–Marquardt

| Stage | Parameters floating | Purpose |
|-------|---------------------|---------|
| 0 | Intensities only | Stabilise amplitudes with positions/widths fixed |
| 1 | Widths + intensities | Adjust envelope while keeping centres fixed |
| 2 | Positions (± bounds) + widths + intensities | Refine centres if `fix_positions` is false |
| 3 | *Skipped* in 2D path | Reserved for 3D legacy variants |
| 4 | Full refinement | Joint optimisation with adaptive damping |

Adaptive damping constants, step multipliers, and stopping conditions are ported from the reference C++ implementation. Analytical intensity derivatives are used; other derivatives rely on relative perturbations (`p × 1.000001`).

### 2.4 Bounds Strategy

- Position bounds: `pos_f1 ± max(1.5 × FWHM_f1, 0.1 ppm)`, `pos_f2 ± max(1.5 × FWHM_f2, 0.01 ppm)`.
- Linewidth bounds: `[0.5, 5.0] × median` for F1 (Lor/Gau), `[0.5, 2.0] × median` for F2.
- Intensities clamped to `[0.01, 5.0] ×` initial estimate.

`fix_positions` / `fix_linewidths` flags (controlled by the GUI and parameter manager) shrink the corresponding bounds to a tiny epsilon window.

### 2.5 Outputs

Each fitted peak reports:

- `pos_f1`, `pos_f2`, Lorentzian/Gaussian FWHM in both dimensions.
- `intensity` (volume), `height` (evaluated at centre), `volume`.
- Quality metrics: total iterations, χ², R² computed on masked points.

The fitted 2D surface and per-peak surfaces can be reconstructed for plotting and ML logging.

## 3. Overlap Resolution & Staged 1D Fitting

`core/overlap_resolver_engine.py` provides an alternative when PS2D is disabled or the spectrum is effectively 1D.

1. **Model selection** (`ModelSelectionEngine`): tries multiple peak counts and scores them via AIC/BIC.
2. **Staged fitting** (`staged_fitting_strategy.py`): three warm-start stages – positions locked, constrained (± window fraction), then full freedom – each using SciPy’s `curve_fit` with Voigt models from `multi_peak_models.py`.
3. **QA**: jackknife resampling for uncertainties and correlation analysis to highlight parameter coupling.

Stage-2 position windows derive from the fitting windows supplied by the GUI/processor (≈±0.03 ppm in 1H, ≈±0.4 ppm in 15N regions).

## 4. Consensus & Simplified Parameter Modes

- **Simplified parameter manager** (`utils/simplified_parameter_manager.py`): exposes a handful of sliders (`sensitivity`, `window_scale`, `quality_target`) and maps them to the legacy 25+ parameters. Supports nucleus-specific defaults and adaptive thresholds.
- **Consensus fitting** (`core/consensus_fitting_engine.py`): runs multiple detection/fitting strategies (sequential, simultaneous, iterative) and selects the best model via information criteria and quality scores. Integrates seamlessly with the simplified manager.

## 5. Machine Learning Hooks

- `ml/ps2d_training_collector.py` registers callbacks with the PS2D fitter to capture 90+ features per fit (spectral widths, asymmetry, residual stats, detection context).
- Batch runs (`batch_processing/batch_processor.py`) can persist training samples to `ml_training_data/` for later modelling.

## 6. Quality Assessment

Across engines, lunaNMR tracks:

- **R²** (global or masked, depending on fit region).
- **χ²** and reduced χ² (with degrees-of-freedom adjustments in staged fits).
- **Residual statistics** (RMS, relative residual heatmaps in 2D plots).
- **Correlation warnings** when parameter correlations exceed thresholds (default 0.95).
- **Detection/fitting status flags** for the peak navigator and exports.

## 7. Export Schema

Processed peaks share a common dictionary structure:

```
{
  'assignment': str,
  'pos_f1': float, 'pos_f2': float,
  'lw_gau_f1': float, 'lw_lor_f1': float,
  'lw_gau_f2': float, 'lw_lor_f2': float,
  'intensity': float, 'height': float, 'volume': float,
  'r_squared': float,
  'amplitude': float,
  'method': '2d_simultaneous_multi_peak' | 'consensus' | 'staged_1d',
  'region_2d': {...} (optional),
  'x_fit': {...}, 'y_fit': {...}
}
```

Maintaining this schema ensures GUI tables, exports, and ML pipelines remain compatible regardless of which engine produced the results.

## 8. Further Reading

- `docs/AUTOMATED_FITTING_IMPLEMENTATION.md` for simplified/consensus design decisions.
- `docs/batch_processing_ml_guide.md` for batch-ML integration details.
- `core/ps2d_2d_fitter.py` docstrings for in-code references.

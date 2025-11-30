# Fine-Tuning PS2D 2D Voigt Fitting – Developer Playbook

This guide consolidates the tunable parameters, code entry points, and experimentation ideas for improving the PS2D multi-peak fitting pipeline. File and line references point to the current codebase (`lunaNMR_v1o0`). Quantitative suggestions below note expected qualitative changes in convergence, fit stability, and runtime so you can predict trade-offs before running sweep experiments.

---

## 1. Key Files & Control Flow

| Purpose | File & Section |
|---------|----------------|
| Entry point for PS2D fitting | `lunaNMR/core/ps2d_2d_fitter.py` (`Ps2dMultiPeakFitter2D.fit_multi_peak_2d`) |
| Parameter mapping, Faddeeva calls, derivatives | `lunaNMR/core/ps2d_style_fitter.py` (Voigt primitives) |
| Mask construction / overlap grouping | `lunaNMR/core/ps2d_data_selector.py`, `ps2d_exact_overlap_detector.py` |
| GUI toggles & user-facing knobs | `lunaNMR/gui/main_gui.py:205-215`, `gui_components.py` |
| Simplified parameter integration | `lunaNMR/utils/simplified_parameter_manager.py` |
| Overlap resolver fallback (staged 1D) | `lunaNMR/core/overlap_resolver_engine.py` |
| ML collector | `lunaNMR/ml/ps2d_training_collector.py` |

Understanding how these blocks interact will help you plug in experiments without breaking the result schema expected by downstream components.

---

## 2. Parameter Bounds & Movement Constraints

### 2.1 Position Bounds (F1 & F2)
- **Code**: `ps2d_2d_fitter.py:257-276`
- **Default**: `pos_f1 ± max(1.5 × FWHM_f1, 0.1 ppm)` and `pos_f2 ± max(1.5 × FWHM_f2, 0.01 ppm)`
- **Expectation**: With reasonably accurate initial guesses, fewer than 20% of peaks reach the boundary. Enlarging the bound increases the risk of drifting into neighbouring clusters but may recover mis-centred detections.

**Ideas to test**:
1. **Looser windows** – raise the multiplier from `1.5` to `2.0` or `2.5`. Expect up to +30% more positional drift, which can help unresolved doublets but may elongate convergence (iterations +10-15%). Monitor residual maps for new artefacts.
2. **Tighter windows** – reduce to `1.2` or even `1.0`. This should stabilise fits for clean datasets (R² spread narrows) but might reject peaks when detection seeds are off by >0.05 ppm; you may see pragmatic-success flags drop.
3. **Data-driven bounds** – scale bounds by detection confidence `c` (0–1): e.g. `margin = base × (0.8 + 0.4 × (1 - c))`. Low confidence (`c=0.2`) allows ~1.1× larger movement; high confidence (`c=0.9`) tightens to ~0.9×. Expect fewer outliers after Stage 2.
4. **Asymmetric freedom** – for HSQC, allow `pos_f2` to move more (e.g. `ΔF2_high = 2.5 × FWHM` and `ΔF2_low = 1.5 × FWHM`). This addresses drift along 1H where shimming errors dominate while keeping F1 anchored.

### 2.2 Linewidth Bounds
- **Code**: same block (`median × {0.5, 5.0}` for F1; `median × {0.5, 2.0}` for F2).
- **Expectation**: If true linewidth ratios vary widely, Stage 1 will hit bounds, flattening peaks.

**Adjustments**:
- Raise F2 upper multiplier from `2.0` to `3.0` for datasets with field inhomogeneity (expect R² improvements of 0.02–0.05 at cost of occasional over-broad peaks).
- Set lower multiplier to `0.3` if narrow peaks collapse (R² dips due to underestimation). Be cautious: values <0.3 may cause numerical instability.
- Link bounds to overlap severity (from `ps2d_exact_overlap_detector.py`): `upper = median × (2 + overlap_score)` where `overlap_score` ∈ [0,1]. Overlaps get extra width budget; isolated peaks stay tight.

---

## 3. LM Optimiser Controls

### 3.1 Damping Parameters
- **Code**: `ps2d_style_fitter.py:63-71`
  - `LAMBDA_INIT = 0.001`
  - `LAMBDA_UP = 10.0`
  - `LAMBDA_DOWN = 0.1`
  - `MAX_ITER = 250`
  - `CONV_TOL = 1e-6`, `SLOW_PROGRESS_TOL = 1e-7`, `SLOW_PROGRESS_LIMIT = 15`

**Levers & expected behaviour**:
- Lower `LAMBDA_INIT` to `1e-4` → larger initial steps; expect faster convergence when starting far from optimum but higher overshoot risk. Monitor for oscillatory χ² traces.
- Increase `LAMBDA_INIT` to `5e-3` when derivatives are noisy; the solver will take smaller steps, reducing divergence at the cost of ~15–20% more iterations.
- Reduce `LAMBDA_UP` from `10` to `4` to prevent the damping from exploding after a failed step; this keeps progress smoother but might slow recovery if the model needs a large damping increase.
- Increase `LAMBDA_DOWN` from `0.1` to `0.3` to allow faster step-length growth after improvements.
- Raise `MAX_ITER` to `400` for large clusters (>6 peaks). Convergence probability increases ~5–10% while runtime scales roughly linearly.
- Bump `SLOW_PROGRESS_LIMIT` to `30` if high-precision convergence is needed (e.g., sub-ppm accuracy). Combine with logging to ensure the solver isn't spinning.

### 3.2 Derivative Step Size
- **Code**: `ps2d_style_fitter.py:74` (`DERIV_STEP_MULTIPLIER = 1.000001`)
- **Interpretation**: Finite-difference step ≈ `1e-6 × parameter`.

**Experiments**:
- Increase to `1.00001` (step ≈ `1e-5`). Expect smoother Jacobians and fewer numerical artefacts for very broad peaks, but small features may get under-resolved.
- Decrease to `1.0000001` (step ≈ `1e-7`). Sharper gradients leading to aggressive LM updates; combine with higher `LAMBDA_INIT` to prevent divergence.

### 3.3 Stage Sequencing
- **Stages**: `ps2d_2d_fitter.py:272-425`

**Variants & effects**:
- Run Stage 1 twice with updated linewidth bounds (first pass `median × 2`, second pass `median × 1.2`). Helps when initial width guesses are poor; expect Stage 2 to start closer to optimum.
- Insert a custom Stage 3 that fixes intensities and widths but frees positions with smaller margins (`±0.3 × FWHM`); may improve centre accuracy by focusing the optimiser.
- Skip Stage 0 for high-SNR isolated peaks to save time (integration runtime -10% with negligible accuracy change).

### 3.4 Alternative Solvers
- Replace LM with `scipy.optimize.least_squares(..., method='trf', loss='soft_l1')`. Soft-L1 mitigates outliers; expect 1.1–1.3× runtime but better robustness to baseline artefacts.
- Explore trust-region reflective with `max_nfev=500`. Good for large bound changes; test on multi-peak clusters where LM diverges.
- GPU: If you implement CuPy-based Faddeeva, expect >5× speedup on >100k-point masks but budget time to convert arrays and handle precision.

---

## 4. Elliptical Mask & Data Selection

- **Code**: `ps2d_data_selector.py`, invoked in `core_integrator.py:2375-2389`
- **Default radii**: `radF1 = 0.4`, `radF2 = 0.04` ppm (~33% tighter than original PS2D C++ defaults of 0.6/0.06)

**Tuning strategies**:
1. Multiply radii by 1.25 for clusters with >4 peaks. Expect more context captured (R² up to +0.03) but slightly higher runtime; check that masks don’t overlap neighbouring clusters.
2. Shrink radii to 0.3/0.03 for high-SNR, well-separated peaks to reduce noise influence.
3. Add anisotropic scaling: `radF1 = base × (1 + 0.5 × overlap_score)`, `radF2 = base × (1 + 1.0 × overlap_score)`. Heavily overlapped clusters get larger masks especially in F2.
4. Apply intensity thresholding inside masks, e.g. keep top 85% contour levels. Lowers χ² sensitivity to noise but may clip broad shoulders.
5. Stage-wise masks: Use 0.5/0.05 for Stage 0-1 then tighten to 0.35/0.035 for Stage 4. Should stabilise final refinement.

---

## 5. Initial Guesses & Normalisation

- **Initial parameters**: computed in `core_integrator.py:2321-2367`
- **Normalisation**: `region['intensity'] / max_intensity`

**Enhancements**:
- Try mixed normalisation: divide by `95th percentile` instead of absolute max to reduce influence of outliers.
- For Stage 0 intensities, start at `0.5 × detected peak height` instead of `0.8` to avoid overshoot; helps clusters where detection overestimates amplitude.
- Introduce jitter: `params += Normal(0, 0.02 × bounds_width)` before Stage 2 for Monte Carlo restarts. Run 3–5 restarts and choose best R²; expect improved accuracy on pathological clusters.
- Use 1D pre-fit: Fit F2 cross-section with `ps2d_multi_peak_fitter.py` results to seed positions/widths (especially when `fix_positions` is False).

---

## 6. PS2D Style Fitter (Numerics)

- **Voigt evaluation**: `ps2d_style_fitter.py:81-137`
- **Jacobian**: `ps2d_style_fitter.py:211-319`

**Tweaks**:
- Cache per-peak Voigt surfaces (`y_peak`) between derivative evaluations to cut redundant computation (already partly done, but can be extended across stages).
- Switch Numba to `parallel=True` for loops in `compute_multi_voigt_jacobian_2d` when cluster size ≥4; check for thread-safety of `wofz` (requires object mode).
- Introduce weighting: pass a weight matrix `w(f1,f2)` to the model function (`model_function` in `ps2d_2d_fitter.py`) to down-weight low-intensity edges, e.g. `w = exp(-(distance/radius)^2)`.
- Evaluate gradient scaling: multiply Jacobian columns by `1 / (parameter_scale)` before solving LS to improve conditioning.

---

## 7. GUI & User Parameters

- **Toggle mapping**: `main_gui.py:205-215`
- **Simplified manager**: `simplified_parameter_manager.py`

**Experiments**:
- Add advanced sliders to the GUI (developer mode) to expose `position multiplier`, `λ init`, `mask radii`. Present default, min, max values with tooltips describing expected effects.
- Store per-peak overrides in the integrator (e.g. if user manually tightens a bound on a difficult peak, persist to JSON for later runs).
- Provide a “retry with aggressive settings” button that automatically sets `λ_init=1e-4`, `Λ_up=6`, `mask×1.3`.

---

## 8. Overlap Resolver Fallback

- **Usage**: triggered when PS2D fails or is disabled.
- **Bounds**: Stage 2 uses ±10–20% of fitting window (`staged_fitting_strategy.py:120-220`).

**Fine-tuning**:
- Align PS2D and staged bounds by passing PS2D-derived linewidths and centres to the staged fitter; helps comparative analysis.
- Replace default tolerances with data-driven ones: e.g. `position_delta = max(0.01, 0.5 × FWHM)` for 1H. Improves stability on narrow features.
- Combine staged results with PS2D as verification: if both agree within 0.01 ppm, flag peak as reliable; if not, prompt for manual review.

---

## 9. Batch Processor Hooks

- **File**: `lunaNMR/batch_processing/batch_processor.py`
- **Presets**: S/N thresholds, expected peaks, quality thresholds.

**Advanced ideas**:
- Add CLI flags to override PS2D parameters per run (`--ps2d-pos-multiplier`, `--ps2d-lambda-init`). Useful for implementing grid sweeps.
- Implement result logging (CSV) capturing runtime, iterations, damping trajectories, residual stats for each cluster to analyse offline.
- Add a “stress-test mode” that runs each spectrum under multiple parameter combinations (e.g. 3×3 grid) and records best-performing configuration.

---

## 10. Exotic Experiments

1. **Adaptive Bounds** – After Stage 1, recompute median linewidths and shrink bounds by 20% if parameters stayed interior. Expected effect: faster Stage 4 convergence, potential risk of constraining legitimate shifts.
2. **Temperature Scaling** – Multiply bound widths by `0.9` each stage to gradually tighten (simulated annealing). Anticipated outcome: more precise final centres but higher chance of early rejection if initial stage is poor.
3. **Penalty Terms** – Augment χ² with `(linewidth_i - median)^2 × α` (`α` ~ 1e3). Encourages consistent linewidths across a cluster; may suppress true heterogeneity.
4. **Multi-Resolution Fitting** – Downsample the mask to 50% resolution for Stage 0-1, then revert to full grid for Stage 4. Expected runtime drop ~30% with similar accuracy.
5. **Peak Ordering** – Sort peaks by intensity before flattening params to improve conditioning. Should reduce degeneracy when merging results.
6. **Hybrid 1D/2D** – Run 1D fits on F2 cross-sections between stages; use derived positions to nudge Stage 2 initial parameters (weighted average). Helps when peaks overlap heavily in one dimension but are separable in the other.
7. **Randomised Restarts** – Execute LM with 3–5 random seeds (perturb positions ±0.03 ppm). Keep solution with best masked R². Expect improved robustness at cost of `×N` runtime.
8. **GPU Faddeeva** – Investigate GPU libraries (`cupyx.scipy.special.wofz`). Would require reworking Numba wrappers; payoff is significant for large masks (>200×200 grid).
9. **Alternative Loss Functions** – Try Huber or Tukey losses via TRF solver to handle outliers; can reduce sensitivity to spikes or baseline ripples.

---

## 11. ML-Assisted Improvements

### Current Capabilities
- `ps2d_training_collector.py` stores per-fit features: linewidths, residual stats, overlap metrics, detection confidence, LM iteration counts.
- Batch processor can export datasets (CSV/Parquet) for offline modelling.

### How ML Can Help
1. **Parameter Suggestions** – Train regression models to predict optimal `(λ_init, position_multiplier, mask radii)` given features such as cluster size, S/N, linewidth variance. Expect to cut manual tuning cycles.
2. **Failure Prediction** – Classification model (e.g., gradient boosting) to flag fits likely to fail; take preventive actions (increase bounds, re-run 1D pre-fit).
3. **Mask Customisation** – Use ML to predict ideal radii per cluster (target variables derived from human-validated fits). Should reduce manual tweaking.
4. **Post-Fit Corrections** – ML could adjust centres slightly based on residual asymmetry. Must ensure corrections remain physically plausible; treat as advisory.
5. **Cluster Typing** – Multi-class model to tag overlap patterns (doublet, triplet, broad cluster). Use tags to load tailored parameter templates or stage sequences.

### Limitations / Cautions
- ML requires labelled ground truth; manual curation is expensive. Start with small, high-quality datasets.
- Spectrometer drift may change optimal parameters; incorporate metadata (field strength, temperature) and consider continual learning.
- Real-time application should stay deterministic: ML suggestions should feed into deterministic PS2D pipelines rather than replace them.
- Over-reliance on ML can mask underlying numerical issues; treat it as a co-pilot, not an autopilot.

### Next Steps for ML Integration
- Extend collector to log LM damping trajectory (`λ` per iteration) for better training signals.
- Build notebooks (`report/`) comparing feature importance vs. convergence quality to prioritise which parameters to learn.
- Prototype a CLI assistant that sweeps small parameter grids guided by ML priors and reports back the highest-R² configuration per dataset.

---

Fine-tuning PS2D is an iterative process. Use version-controlled configuration files or notebooks to track experiments, and document findings directly in this playbook as you discover settings that consistently improve fit quality for your datasets.

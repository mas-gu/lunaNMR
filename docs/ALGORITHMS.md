# Algorithm Reference

Prominence-based peak detection with top-contour centroid refinement, PS2D simultaneous 2D Voigt fitting (5-stage Levenberg-Marquardt), and hierarchical overlap clustering.

Quality (R²): Excellent ≥0.9, Good 0.8–0.9, Fair 0.5–0.8, Poor 0.2–0.5, Failed <0.2.

---

## 1. Peak Detection

`core/enhanced_peak_picker.py`

**Pipeline**: noise estimation → candidate detection → overlap clustering → centroid refinement.

- **Noise**: corner sampling / MAD.
- **Candidates**: local maxima via `scipy.ndimage.maximum_filter` + prominence/S-N thresholding.
- **Search window** default **0.070 ppm in both 1H and 15N** (`main_window.search_window_x/y`). Wide 15N windows are harmful on low-S/N spectra (peaks fit to noise → Height 0 → residues dropped).

### 1.1 Top-Contour Centroid

`calculate_top_contour_center` (picker) / `_calculate_top_contour_center` (`core_integrator.py`). Sub-pixel positioning for flat-top peaks.

1. Find max intensity I_max.
2. Select pixels within ±5% band `[0.95·I_max, 1.05·I_max]` (`intensity_band=0.05`).
3. Intensity-weighted center `x_c = Σ(x_i·I_i)/Σ I_i` (same for y).
4. Clip shift to GUI `max_shift_x/y`; fall back to pixel max if band empty.

GUI defaults (`main_window`): `centroid_window_x_ppm=0.01` (1H), `centroid_window_y_ppm=0.1` (15N).

---

## 2. PS2D Simultaneous 2D Voigt Fitting

`core/ps2d_2d_fitter.py` — `Ps2dMultiPeakFitter2D.fit_multi_peak_2d`. Optimizer `PS2DStyleFitter` in `ps2d_style_fitter.py`.

### 2.1 Data Selection
Union of elliptical windows around each peak: `(ΔF1/radF1)² + (ΔF2/radF2)² ≤ 1`.
Defaults (`ps2d_config.py`, 15N): radF1=0.4, radF2=0.04. 13C: radF1=0.15, radF2=0.04.

### 2.2 Parameter Vector (8 per peak)
`[pos_f1, lw_lor_f1, lw_gau_f1, pos_f2, lw_lor_f2, lw_gau_f2, intensity, spare]` (spare fixed 0).

Internal Voigt units: σ = FWHM/2.3548 (Gaussian), γ = FWHM/2 (Lorentzian); profile via Faddeeva `wofz`.

### 2.3 Five-Stage Levenberg-Marquardt

| Stage | Floating | Notes |
|-------|----------|-------|
| 0 | intensities only | DISABLED (caused intensity collapse) — fitting starts at Stage 1 |
| 1 | linewidths + intensities (positions fixed) | skipped if `fix_linewidths=True` |
| 2 | positions + linewidths + intensities | skipped if `fix_positions=True` |
| 3 | reserved / unused | |
| 4 | all parameters (global refine) | fixed params restored exactly after `np.clip` bounds |

GUI **Fix Positions**→skip Stage 2; **Fix Linewidths**→skip Stage 1. Stages 1/4 support optional L/G-ratio and intensity-ratio soft penalties (applied only for too-close clusters).

**Optimizer** (`PS2DStyleFitter`): additive damping `α + λ·I`; λ init 1e-4, ×5 on rejection, ×0.2 on acceptance. Analytical Jacobian (`compute_multi_voigt_jacobian_2d`). Parameter normalization `self.scales = max(|params|, 1e-8)` (Hessian condition ~10¹⁰→10³–10⁴); Jacobian scaled by `scales`, covariance denormalized via `outer(scales, scales)`.

### 2.4 Bounds
- **Positions**: `pos ± max(0.04·FWHM, pos_margin)`; `pos_margin` 15N F1=0.02, F2=0.01 (titration mode passes wider margins). Cascade mode additionally clips to `±max_drift` of the reference position.
- **Linewidths**: min 0.03 (F1/15N), 0.005 (F2/1H); max = `min(config_max, learned_median + α·MAD)` from two-pass stats, config_max 15N F1=0.8, F2=0.08.
- **Intensities**: `[0.001, 5.0] × max cluster intensity`.

### 2.5 Output (per peak)
`pos_f1/pos_f2`, `lw_lor_f1/lw_gau_f1/lw_lor_f2/lw_gau_f2` (FWHM), `intensity` (volume), `height`, `r_squared`, `method='2d_simultaneous_multi_peak'`.

---

## 3. Overlap Clustering & Routing

`core_integrator.py:identify_overlap_clusters` — agglomerative closest-pair merging (NOT networkx / connected-components). Called once before fitting; results identical in sequential and parallel mode.

**Overlap test** (`_check_peaks_overlap`): elliptical distance `√((Δx/tx)² + (Δy/ty)²) ≤ 2.0`, with thresholds `overlap_threshold_x/y` (15N: x=0.04 (1H), y=0.4 (15N); 13C: x=0.04, y=0.15).

**Merging**: start singletons → repeatedly merge the closest overlapping pair whose merged **diameter** stays within `max_cluster_diameter_x=0.20` (1H), `max_cluster_diameter_y=1.5` (15N); oversized clusters (>`max_cluster_size`, default 20) are split. Diameter/size constraints deliberately bound transitive closure so it can't form mega-clusters. Result = disjoint partition, each peak in exactly one cluster.

**Routing**: `needs_2d = len(cluster) > 1`.
- size 1 → 1D consensus / staged Voigt (`enhanced_peak_fitting`).
- size >1 → `fit_overlap_group_2d` → PS2D; falls back to 1D cross-sections on failure.

**Parallel** (`parallel_voigt_processor.py::fit_all_peaks_parallel`): same clustering runs once, then whole clusters distribute across workers; results consolidated by integer `peak_number`. ~2.7× speedup (6 cores). **Two-pass linewidth learning**: pass 1 fits isolated peaks → median linewidths (+MAD); pass 2 seeds multi-peak clusters with them (falls back to `ps2d_config.py` defaults).

---

## 4. Quality Assessment

`_assign_quality_from_r2` (`core_integrator.py`, applied to PS2D and 1D results alike): Excellent ≥0.9, Good ≥0.8, Fair ≥0.5, Poor ≥0.2, else Failed. Colors: green (Excellent/Good), yellow (Fair), red (Poor/Failed/unfitted).

**Acceptance**: PS2D pragmatic `r_squared > 0.2` (or formal convergence / chi² reduction); 1D `min_r_squared` default 0.3.

**Failed peaks** stored with placeholders (`fitted=False, r_squared=0, intensity=0, quality='Failed'`) to keep 1:1 peak_list↔results index mapping for the Navigator.

---

## 5. Parameter Management

`utils/parameter_manager.py` (+ `simplified_parameter_manager.py`).

**Simplified mode** (`use_simplified_mode`) exposes 3 knobs → maps to legacy params via nucleus-adaptive thresholds:
- `sensitivity` (legacy, no longer affects detection)
- `window_scale` (default 1.0) — fitting-window sizing
- `quality_target` (default 0.85) — feeds `min_r_squared`

GUI spinbox values for `centroid_window_x/y_ppm` **override** simplified defaults via `update_from_gui_variables` → `get_effective_parameters` → `integrator.gui_params`.

---

## 6. Complexity

- Detection O(N) data points; centroid O(W·H) per peak.
- PS2D O(I·n·M) (iterations × cluster peaks × masked points); bottleneck for large clusters. Numba JIT gives 3–5× on the 2D model/Jacobian.

---

## 7. File Formats

Bruker `.2ii`/`.2rr`, Varian/Agilent `.fid`/`.ft`, NMRPipe `.ft`/`.pipe`, SPARKY `.ucsf` — loaded via `nmrglue` with automatic axis-ordering detection.

---

## 8. Downstream Models (Modules)

Fitted series feed DynamiXs and Kd (GUI Modules menu / CLI).

| Analysis | Model |
|----------|-------|
| T₁ / T₂ | Mono-exponential decay |
| Methyl T₂ | Shared-amp bi-exp `I(t)=0.5·A·[e^(−t/T2a)+e^(−t/T2b)]`, fit in (A, T2_avg, ΔT2); T2a (slow/TROSY) is the reportable T₂ |
| hetNOE | `I_sat/I_unsat` per residue (+ QC plot vs residue) |
| Spectral density / model-free | Reduced spectral density → Lipari-Szabo (single/dual field) |
| Kd — CSP | 1:1 quadratic isotherm `Δδ = Δδ_max·fraction_bound(L,P0,Kd)` |
| Kd — intensity | `I = I_inf + (I0−I_inf)·e^(−L/Kd)` |

All Kd fits bound Kd ≥ 0. Headless equivalents: [CLI.md](CLI.md).

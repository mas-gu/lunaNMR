# Algorithm Reference

Prominence-based peak detection with top-contour centroid refinement, PS2D simultaneous 2D Voigt fitting (5-stage Levenberg-Marquardt), and hierarchical overlap clustering.

Quality (R²): Excellent >0.9, Good >0.8, Fair >0.6; below that Poor above S/N 5, else Poor ≥0.2, Failed <0.2 (a near-noise peak Poor). See §4.

---

## 1. Peak Detection

`core/enhanced_peak_picker.py`

**Pipeline**: noise estimation → candidate detection → overlap clustering → centroid refinement.

- **Noise**: corner sampling / MAD.
- **Candidates**: local extrema via `scipy.ndimage.maximum_filter` + prominence/S-N thresholding, **per sign** (`_detect_peaks_by_threshold`): maxima above +threshold and minima below −threshold, so a negative peak beside a positive lobe is kept.
- **Search window** default **0.070 ppm in 1H; the indirect window is the nucleus's: 0.3 ppm for 15N, 0.11 for 13C** (`ps2d_config` `search_window_f1`; startup, Reset, a nucleus switch and the CLI series, whose indirect default is resolved once spectrum 1 sets the nucleus). A 15N point can be 0.074 ppm (600 MHz hetNOE planes), so 0.070 could not reach a peak one point off its reference. **On every later spectrum of a series** a reference only searches the 0.070 × 0.070 box around its LISTED position, so a peak that has decayed below threshold cannot take a stronger unlisted neighbour (the wide box gave JA1 T2 44LysH a neighbour's peak and a T2 8× its own). Measured on the saturated/unsaturated hetNOE pair at the GUI threshold: 47 residues kept against 41, three lost (a wider box reaches a different maximum).
- **Reference matching** (`core_integrator._detect_peaks_reference_based`): **competitive**. Every (listed peak, candidate) pair inside the search box is ranked by its distance **in linewidths**, `√((Δ1H/radF2)² + (ΔX/radF1)²)` with the PS2D radii of the spectrum's nucleus (15N 0.04/0.4, 13C 0.04/0.15) — not raw ppm, which counts an indirect ppm as a 1H ppm — and the nearest pairs are settled first across all listed peaks, one maximum per peak. A contested maximum goes to the peak nearest it, not to whichever comes first in the list: in list order, a peak with no maximum of its own would take the one its neighbour sits on. A listed peak left without a candidate is read at its listed position; dummies come from the maxima no listed peak took.

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
| 1 | linewidths + intensities (positions fixed) | always runs; with `fix_linewidths=True` the linewidths are held and only intensities float |
| 2 | positions + linewidths + intensities | skipped if `fix_positions=True` |
| 3 | reserved / unused | |
| 4 | all parameters (global refine) | positions held if `fix_positions`, linewidths if `fix_linewidths`; fixed params restored exactly after `np.clip` bounds |

GUI **Fix Positions**→skip Stage 2 and hold positions in Stage 4; **Fix Linewidths**→hold linewidths in every stage (Stage 1 still fits the intensities). A near-noise peak (`weak_peak_indices`, §2.6) has its position and linewidths held in every stage. Stages 1/4 support optional L/G-ratio and intensity-ratio soft penalties (applied only for too-close clusters).

**Optimizer** (`PS2DStyleFitter`): additive damping `α + λ·I`; λ init 1e-4, ×5 on rejection, ×0.2 on acceptance. Analytical Jacobian (`compute_multi_voigt_jacobian_2d`). Parameter normalization `self.scales = max(|params|, 1e-8)` (Hessian condition ~10¹⁰→10³–10⁴); Jacobian scaled by `scales`, covariance denormalized via `outer(scales, scales)`.

### 2.4 Bounds
- **Positions**: `pos ± max(0.04·FWHM, pos_margin)`; `pos_margin` 15N F1=0.02, F2=0.01 (13C F1=0.01). Titration mode passes a per-step window of 4 × the reference spectrum's median linewidth (capped at 1.0 / 0.2 ppm; `titration_pos_margin` 0.3 / 0.05 before linewidths are learned), and each peak's own **territory** (`ps2d_territories.territory_margins`: 0.49 × the distance to its nearest neighbour, in units of that window) through `fit_multi_peak_2d(pos_margins=...)`, taken exactly with no floor — all-or-nothing per cluster. Cascade mode additionally clips to `±max_drift` of the reference position (15N F1=0.025, F2=0.005) when **Enforce max drift from reference** is on (default); titration mode turns that limit off.
- **Linewidths**: min 0.03 (F1/15N), 0.005 (F2/1H); max = `min(config_max, learned_median + α·MAD)` from two-pass stats (α sample-size dependent, default 5; MAD floored at 5% of the typical linewidth), config_max 15N F1=0.8, F2=0.08. `fit_multi_peak_2d(lw_margins=...)` can also cap each peak's **total** `lw_lor + lw_gau` per dimension (`cap_linewidth_totals`, a projection that scales the pair together and keeps the lor/gau ratio); no lunaNMR series sets it — only the overlay refit does, from a point's `linewidth_caps` in `series_metadata.json`.
- **Intensities**: cluster-relative, signed by the seed (§2.6): `[0.001, 5.0] × max cluster intensity` for a positive peak, the mirror image for a negative one, `[−5, 5] × max` for a near-noise peak.

### 2.5 Output (per peak)
`pos_f1/pos_f2`, `lw_lor_f1/lw_gau_f1/lw_lor_f2/lw_gau_f2` (FWHM), `intensity` (volume), `height` (the signed extremum of the fitted surface), `r_squared`, `method='2d_simultaneous_multi_peak'`, plus `detection_sigma`, `position_pinned`, `intensity_pinned` (§4).

### 2.6 Sign and Near-Noise Peaks
`core_integrator.peak_sign` classifies each peak from the largest-magnitude sample the fit can reach (its territory or position margin, never less than the seed radius):
- **+1 / −1** at 3σ or more: the seed and intensity box carry that sign; the box never contains zero, where the derivatives with respect to position and linewidth vanish.
- **0** (near-noise): position and linewidths are held, the intensity spans ±5 × max, so the fit is linear in intensity and crossing zero costs nothing. A cluster of only near-noise peaks is reported as solved (`near_noise_only`) and labelled Poor, never Failed.

Without a noise level every peak is treated as positive. Mixed-sign clusters are not handled.

---

## 3. Overlap Clustering & Routing

`core_integrator.py:identify_overlap_clusters` — agglomerative closest-pair merging (NOT networkx / connected-components). Called once before fitting; results identical in sequential and parallel mode.

**Overlap test** (`_check_peaks_overlap`): elliptical distance `√((Δx/tx)² + (Δy/ty)²) ≤ 2.0`, with thresholds `overlap_threshold_x/y` (15N: x=0.04 (1H), y=0.4 (15N); 13C: x=0.04, y=0.15).

**Merging**: start singletons → repeatedly merge the closest overlapping pair whose merged **diameter** stays within `max_cluster_diameter_x=0.20` (1H), `max_cluster_diameter_y=1.5` (15N); oversized clusters (>`max_cluster_size`, default 20) are split. Diameter/size constraints deliberately bound transitive closure so it can't form mega-clusters. Result = disjoint partition, each peak in exactly one cluster.

**Routing**: `needs_2d = len(cluster) > 1`.
- size 1 → `enhanced_peak_fitting` → single-peak PS2D fit.
- size >1 → `fit_overlap_group_2d` → PS2D.

There is no 1D cross-section fallback: a 2D fit that fails is reported Failed. Without a peak list to place the peak among its neighbours (`all_peaks_context`), `fit_peak_voigt_2d` returns `{'success': False, 'error': 'no 2D fit context: …'}`, and the reason is logged.

**One pipeline, sequential or parallel** (`parallel_voigt_processor.py::fit_all_peaks_parallel`): the clustering runs once, then whole clusters are fitted as tasks of `_parallel_cluster_worker` — over a `spawn` pool in parallel mode, as direct calls in this process in sequential mode (`ParallelVoigtProcessor(in_process=True)`, also used below 3 peaks). Results are consolidated by integer `peak_number`, and the two modes give identical results; parallel pays off only on large lists (~2.7× on 6 cores). **Two-pass linewidth learning**: pass 1 fits isolated peaks → median linewidths (+MAD); pass 2 seeds multi-peak clusters with them (R² > 0.85, at least 3 good fits; otherwise `ps2d_config.py` defaults).

**Adaptive radius search** (`adaptive_optimizer.py`, between pass 1 and pass 2): a grid over radF1/radF2 at 1.0 / 1.5 / 2.0 × the median linewidth. **Off by default** (`use_adaptive_optimization`); the Expert-Mode box applies to Fit Spectrum only, a series from the dialog never runs it, and the CLI opts in through `series --params`. It writes its radii into the global PS2D config, so a fresh run (no locked clusters or learned statistics, and every series) first calls `restore_base_ps2d_config()` — the base being the nucleus defaults plus the Expert panel's Apply (`mark_ps2d_config_as_base`).

**Cancel**: a `threading.Event` passed down as `stop=` (`utils/cancellation.py`), checked before each cluster, each detection peak and each adaptive-search fit, and polled while waiting on a pool task. A stopped fit or series raises `Cancelled` and returns nothing: no fallback, no partial results, no `series_results_*` folder.

---

## 4. Quality Assessment

One rule labels every fit: `fit_quality_label(r_squared, success, near_noise, snr)` (`utils/quality_categories.py`). The fitters call it (the 2D path for an isolated peak or a cluster member, the parallel worker, the overlay refit); the series output, the results browser and the markers show the label the fit carries and never recompute it. For a cluster member, R² is the cluster's and S/N its own `detection_sigma`.

| Condition (in order) | Label |
|---|---|
| fit did not succeed, or no R² | Failed |
| R² > 0.9 | Excellent |
| R² > 0.8 | Good |
| R² > 0.6 | Fair |
| S/N > 5 | Poor |
| R² ≥ 0.2, or a near-noise peak | Poor |
| otherwise | Failed |

The Kd loaders and the intensity matrices drop a Failed point, so the last two rows decide what is measured. A near-noise peak's R² is ~0 by construction (a held shape against noise) and its near-zero intensity is the measurement, so it is Poor, never Failed. Marker colors: lime (Excellent), green (Good), orange (Fair), red (Poor/Failed/unfitted).

**Acceptance**: PS2D pragmatic `r_squared > 0.2` (or formal convergence / chi² reduction / a near-noise-only cluster); 1D `min_r_squared` default 0.3.

**Failed peaks** stored with placeholders (`fitted=False, r_squared=0, intensity=0, quality='Failed'`) to keep 1:1 peak_list↔results index mapping for the Navigator.

### What acceptance cannot tell you

`success` is `formal_convergence OR r² > 0.2 OR chi² reduced 100×`, and the first arm is the optimiser's own `done_count >= ndone` — "χ² stopped changing" — which is **true immediately on flat noise**. Three positions chosen in genuinely empty regions of a real spectrum all returned `success=True`; at the weakest point of a real titration, 34 of 68 successful fits sat below three sigma.

R² does not cover the gap either. It is a squared **Pearson correlation** (shape only, blind to a fit whose every intensity is twice too large) and it is computed **per cluster**, so a member with nothing at its position rides the R² its healthy neighbours earn — measured, 0.58σ under a cluster R² of 0.9884.

Each fitted peak therefore carries three further fields, stamped once in `fit_overlap_group_2d`:

| field | meaning | `None` when |
|---|---|---|
| `detection_sigma` | fitted height ÷ noise σ | the spectrum has no noise estimate |
| `position_pinned` | the fitted position rests on its bound — the seed was further out than `pos_margin` allows, and the **linewidth stretched** to cover the gap, so the volume is wrong too | the fit did not report its box |
| `intensity_pinned` | the fitted intensity rests on its interval's edge. That floor is **cluster-relative** (`max_initial_intensity * 0.001`), so such a value is a product of its **neighbours'** seeds and its own seeded linewidths | as above |

`None` is not `False`: silence about a bound is not a denial of one.

### Zero is information, and only sometimes

Two states look identical in a series — both report a near-zero height — and they need opposite handling:

| | CSP / positions | I/I₀ (T1, T2, titration) |
|---|---|---|
| peak genuinely **decayed** below the noise | **exclude** — there is no position to measure a shift from | **KEEP — zero is the measurement** |
| fit stopped on its **intensity floor** | exclude | **exclude** — the value came from a bound, not from data |
| **position** stopped on its margin | exclude — the position is the margin | **keep** — the height does not depend on it |

Dropping a decayed residue from an intensity series fabricates a plateau that is not there; keeping one on its intensity floor publishes its neighbours' seeds as a measurement. A pinned position is kept for intensity because, measured on a 600 MHz hetNOE pair, most sit on the 0.020 ppm ¹⁵N margin (a quarter of one ¹⁵N point): widened to 0.05 ppm, the freed peaks' heights moved by a median 0.0%. `quality_categories.usable_points(frame, observable)` is this table as code, and `dynamiXs_Kd/kd_input.load_titration_tidy` applies it — blanking every column for the intensity floor, and only the ppm columns for a decay or a pinned position.

**These flags cross three boundaries and each drops silently**: the four sites that copy fields out of a 2D result onto a peak record (only one is on the series path), `multi_spectrum_processor.TIDY_COLUMNS` (the single list every tidy row is built from), and `_convert_to_integration_format` (which rebuilds `peak_data` from scratch and drops anything it does not name). A field must be present at all three or it reaches `series_analysis_tidy.csv` as an empty column.

---

## 5. Parameter Management

`utils/parameter_manager.py`: one parameter set, with a default, range and type per parameter.

GUI spinbox values (search and centroid windows, fixed positions/linewidths, custom linewidths) **override** the defaults: `main_window.sync_parameters_to_integrator` writes them into `current_params`, and `get_integrator_parameters` groups them into the detection, fitting and `integrator.gui_params` sets.

---

## 6. Complexity

- Detection O(N) data points; centroid O(W·H) per peak.
- PS2D O(I·n·M) (iterations × cluster peaks × masked points); bottleneck for large clusters. Numba JIT gives 3–5× on the 2D model/Jacobian.

---

## 7. File Formats

Bruker `.2ii`/`.2rr`, Varian/Agilent `.fid`/`.ft`, NMRPipe `.ft`/`.pipe`, SPARKY `.ucsf` — loaded via `nmrglue` with automatic axis-ordering detection.

---

## 8. Downstream Models (Modules)

Fitted series feed DynamiXs and Kd (GUI Modules menu / CLI). DynamiXs fits the intensity or volume matrix, where an intensity-pinned or Failed point is left empty.

| Analysis | Model |
|----------|-------|
| T₁ / T₂ | Mono-exponential decay |
| Methyl T₂ | Shared-amp bi-exp `I(t)=0.5·A·[e^(−t/T2a)+e^(−t/T2b)]`, fit in (A, T2_avg, ΔT2); T2a (slow/TROSY) is the reportable T₂ |
| hetNOE | `I_sat/I_unsat` per residue (+ QC plot vs residue) |
| Spectral density / model-free | Reduced spectral density → Lipari-Szabo (single/dual field), fitted by bounded least squares from a grid of starts including Rex > 0; NOE predicted from the calculated T₁; errors NaN without Monte Carlo refits |
| Kd — CSP | 1:1 quadratic isotherm `Δδ = Δδ_max·fraction_bound(L,P0,Kd)` |
| Kd — intensity | `I = I_inf + (I0−I_inf)·e^(−L/Kd)` |

All Kd fits bound Kd ≥ 0. Headless equivalents: [CLI.md](CLI.md).

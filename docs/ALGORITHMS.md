# Algorithm Reference

**TL;DR**: LunaNMR uses prominence-based peak detection with top-contour centroid refinement, PS2D simultaneous 2D Voigt fitting (5-stage Levenberg-Marquardt), and hierarchical overlap clustering. Quality thresholds: Excellent (R²≥0.9), Good (0.8-0.9), Fair (0.5-0.8), Poor (0.2-0.5), Failed (<0.2).

---

## 1. Peak Detection

### 1.1 Pipeline
```
Noise Estimation → Candidate Detection → Overlap Clustering → Position Refinement
```

**Noise estimation**: Corner sampling or MAD (Median Absolute Deviation).

**Candidate detection**: `scipy.ndimage.maximum_filter` with prominence thresholding.

**Overlap clustering**: Hierarchical graph-based clustering using overlap thresholds:
- **15N**: ±0.04 ppm (F2), ±0.4 ppm (F1)
- **13C**: ±0.04 ppm (F2), ±0.1 ppm (F1)

Clusters are disjoint partitions—each peak belongs to exactly one cluster.

### 1.2 Top Contour Centroid

**Purpose**: Sub-pixel positioning for flat-top peaks.

**Algorithm** (`_calculate_top_contour_center`):
1. Find peak maximum intensity I_max
2. Select pixels with I ∈ [0.95×I_max, 1.05×I_max]  (±5% intensity band)
3. Calculate intensity-weighted center:
   ```
   x_c = Σ(x_i × I_i) / Σ(I_i)
   y_c = Σ(y_i × I_i) / Σ(I_i)
   ```
4. Apply safety constraint: clip shift to GUI-specified max_shift_x/y

**GUI Parameters** (`Peak Centroid Detection`):
- `centroid_window_x_ppm`: Max shift in F2 (default: 0.01 ppm)
- `centroid_window_y_ppm`: Max shift in F1 (default: 0.1 ppm)

**Output**: Console prints "Centroid shift: Δ=X.XXXX ppm from pixel max" when shift >0.001 ppm.

**Complexity**: O(W×H) where W,H are search window dimensions.

---

## 2. PS2D Simultaneous 2D Voigt Fitting

### 2.1 Data Selection

Elliptical windows around each peak in overlap group:
```
(ΔF1/radF1)² + (ΔF2/radF2)² ≤ 1
```
**Defaults**: radF1=0.4 ppm (15N), radF2=0.04 ppm (1H).

### 2.2 Parameter Vector (8 parameters × n peaks)
```
[pos_f1, lw_lor_f1, lw_gau_f1, pos_f2, lw_lor_f2, lw_gau_f2, intensity, spare]
```

**FWHM to internal units**:
- Gaussian: σ = FWHM / √(8 ln 2) ≈ FWHM / 2.3548
- Lorentzian: γ = FWHM / 2

### 2.3 Five-Stage Levenberg-Marquardt

| Stage | Parameters Floating | fix_positions | fix_linewidths |
|-------|---------------------|---------------|----------------|
| 0 | Intensities only | N/A | N/A |
| 1 | Widths + intensities | N/A | Respected (Stage skipped if True) |
| 2 | Positions + widths + intensities | **Skipped if True** | Respected |
| 3 | (Reserved, unused) | N/A | N/A |
| 4 | All parameters | **Positions fixed** | **Widths fixed** |

**Absolute Constraint** (when `fix_positions=True`):
- Stage 2 is entirely skipped
- Stage 4 adds positions to `fixed_params` dictionary
- After `np.clip()` bounds enforcement, fixed parameters are restored to exact values
- **Result**: Zero position drift (within float precision ~10⁻¹⁵)

**Damping**: Adaptive λ starting at 0.001, multiplied by 10 on rejection, divided by 10 on acceptance.

**Derivatives**: Analytical for intensity (∂y/∂A = y/A), finite difference for others.

### 2.4 Bounds

**Positions**:
```
pos_f1 ∈ [initial ± max(1.5×FWHM_f1, 0.1 ppm)]
pos_f2 ∈ [initial ± max(1.5×FWHM_f2, 0.01 ppm)]
```

**Linewidths**:
```
F1: [0.5, 5.0] × median_FWHM
F2: [0.5, 2.0] × median_FWHM
```

**Intensities**: [0.01, 5.0] × initial estimate.

### 2.5 Output
```python
{
  'pos_f1': float, 'pos_f2': float,
  'lw_lor_f1': float, 'lw_gau_f1': float,  # FWHM units
  'lw_lor_f2': float, 'lw_gau_f2': float,
  'intensity': float,  # Volume
  'height': float,     # Peak amplitude at center
  'r_squared': float,
  'method': '2d_simultaneous_multi_peak'
}
```

---

## 3. Overlap Resolution Strategy

**Overlap detection**: Two-circle touching test
```
(|Δx| ≤ 2×threshold_x) AND (|Δy| ≤ 2×threshold_y)
```

**Cluster formation**:
1. Build overlap graph (peaks = nodes, overlaps = edges)
2. Find connected components using `networkx`
3. Each component is a disjoint cluster

**Fitting**:
- Single-peak clusters → 1D cross-section fitting
- Multi-peak clusters → PS2D 2D simultaneous fitting (max 6 peaks/cluster)

**Parallel mode**: Distributes clusters (not individual peaks) across workers. Identical clustering to sequential mode.

---

## 4. Quality Assessment

### 4.1 Quality Hierarchy (Centralized at `core_integrator.py:4297-4317`)

| Quality | R² Range | Color |
|---------|----------|-------|
| Excellent | ≥ 0.9 | Green |
| Good | [0.8, 0.9) | Green |
| Fair | [0.5, 0.8) | Yellow |
| Poor | [0.2, 0.5) | Red |
| Failed | < 0.2 | Red |

### 4.2 Acceptance Thresholds

**PS2D 2D fitting**: R² > 0.2 (pragmatic threshold for difficult overlaps)

**1D fitting**: R² > 0.5

### 4.3 Failed Peak Storage

Peaks that fail fitting are stored with placeholder values:
```python
{
  'quality': 'Failed',
  'fitted': False,
  'r_squared': 0.0,
  'intensity': 0.0,
  'method': 'none',
  'failure_reason': 'Fitting did not converge or failed acceptance criteria'
}
```

This maintains 1:1 index mapping between peak_list and fitted_results for Navigator compatibility.

---

## 5. Parameter Management

### 5.1 Simplified Mode (Default)

3-5 core parameters:
- `sensitivity`: Detection threshold multiplier (0.1-0.9)
- `window_scale`: Fitting window sizing (0.1-10.0)
- `quality_target`: Target R² (0.5-0.95)

Maps to 25+ legacy parameters via nucleus-specific adaptive thresholds.

### 5.2 GUI Parameter Override

**Critical**: When simplified mode is enabled, GUI values for `centroid_window_x_ppm` and `centroid_window_y_ppm` **override** calculated defaults.

**Flow**:
```
GUI spinboxes → ParameterManager.update_from_gui_variables()
             → current_params['centroid_window_x_ppm']
             → ParameterManager.get_effective_parameters()
             → Override simplified_manager values
             → integrator.gui_params
```

---

## 6. Parallel Processing

**Cluster-based workflow**:
1. Call `identify_overlap_clusters()` **once** (deterministic, sequential)
2. Distribute clusters across workers
3. Each worker processes entire clusters (not individual peaks)
4. Consolidate results using integer `peak_number` matching (not float coordinates)

**Performance**: 2.7× speedup with 6 cores.

**Output**: Identical to sequential mode (same clustering, same R² values).

---

## 7. File Format Support

**Bruker**: `.2ii`, `.2rr` (TopSpin)
**Varian/Agilent**: `.fid`, `.ft`
**NMRPipe**: `.ft`, `.pipe`
**SPARKY**: `.ucsf`

Loaded via `nmrglue` with automatic axis ordering detection.

---

## 8. Complexity & Performance

**Peak detection**: O(N) where N = number of data points
**Centroid calculation**: O(W×H) per peak (W,H = window size)
**PS2D fitting**: O(I × n × M) where I = iterations, n = peaks in cluster, M = masked data points
**Parallel consolidation**: O(P) where P = total peaks

**Bottleneck**: PS2D 2D fitting for large overlap clusters (>6 peaks). Numba JIT provides 3-5× speedup.

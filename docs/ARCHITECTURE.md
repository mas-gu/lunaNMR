# Architecture Guide

**TL;DR**: MVC architecture with GUI (tkinter), processors (workflow orchestration), and core engines (detection, PS2D fitting, overlap resolution). Parallel mode distributes overlap clusters across workers. All paths share identical result schema for compatibility.

---

## 1. Directory Structure

```
lunaNMR_v1o0/
├── launch_lunaNMR.py           # Entry point
├── docs/                       # Documentation
├── lunaNMR/                    # Main package
│   ├── gui/                    # Tkinter GUI + visualization
│   ├── processors/             # Workflow orchestration
│   ├── core/                   # Detection, fitting, clustering
│   ├── utils/                  # Config, parameters, file I/O
│   ├── ml/                     # Training data collectors
│   └── validation/             # Installation tests
└── modules/                    # Optional add-ons (DynamiXs)
```

---

## 2. Layer Responsibilities

| Layer | Modules | Role |
|-------|---------|------|
| **GUI** | `gui/` | Tkinter UI, plots, Peak Navigator, Voigt Analysis tab |
| **Workflow** | `processors/` | Orchestrate detection + fitting engines for single/batch/series |
| **Core** | `core/` | Peak detection, PS2D 2D fitting, overlap clustering |
| **Utils** | `utils/` | Parameter manager, config, file formats |
| **ML** | `ml/` | Capture PS2D diagnostics for training datasets |

---

## 3. Core Engines

### 3.1 Peak Detection (`core/enhanced_peak_picker.py`)
- Prominence-based detection via `scipy.ndimage.maximum_filter`
- Graph-based clustering using `networkx`
- Top contour centroid refinement (`_calculate_top_contour_center`)

### 3.2 PS2D Stack
- **Data selector** (`ps2d_data_selector.py`): Elliptical masks around clusters
- **Voigt primitives** (`ps2d_style_fitter.py`): Faddeeva-based Voigt profiles (Numba-accelerated)
- **2D fitter** (`ps2d_2d_fitter.py`): 5-stage Levenberg-Marquardt optimizer
- **Overlap detector** (`ps2d_exact_overlap_detector.py`): Two-circle touching test for clustering

### 3.3 Parameter Management (`utils/parameter_manager.py`)
- **Simplified mode**: 3-5 sliders map to 25+ legacy parameters
- **GUI override**: Centroid windows always use GUI spinbox values (not calculated defaults)

---

## 4. Data Flow

### 4.1 Sequential Mode (Default)

```
GUI action
  → processors.single_spectrum_processor
  → core.detect_peaks_sn_native()
      → _detect_peaks_by_threshold()
          → _calculate_top_contour_center() [centroid refinement]
  → core.identify_overlap_clusters() [hierarchical graph clustering]
  → For each cluster:
      Single peak  → core.enhanced_peak_fitting() [1D cross-sections]
      Multi-peak   → core.ps2d_2d_fitter.fit_multi_peak_2d() [2D simultaneous]
  → processors.format_results()
  → GUI tables/plots
```

### 4.2 Parallel Mode

```
GUI action
  → processors.single_spectrum_processor
  → core.identify_overlap_clusters() [ONCE - creates disjoint clusters]
  → parallel_voigt_processor.distribute_clusters()
      → Worker Pool (each processes entire clusters)
          → For each cluster:
              Single peak  → 1D fitting
              Multi-peak   → PS2D 2D fitting
      → Consolidate results (integer peak_number matching)
  → processors.format_results()
  → GUI tables/plots
```

**Key differences**:
- Clustering algorithm **identical** (deterministic)
- Workers receive **entire clusters**, not individual peaks
- Result consolidation uses `peak_number` (integer), not coordinates (float)
- Performance: 2.7× speedup with 6 cores

---

## 5. GUI Architecture

**MVC pattern**:
- **Model**: Integrator state (`core_integrator.py`), peak lists, fitting results
- **View**: `main_gui.py`, `gui_components.py`, `visualization.py`, `spectrum_browser.py`
- **Controller**: GUI callbacks → processors → update model → refresh view

**Key components**:
- **Peak Navigator** (`gui_components.py`): Browse peaks 1-N with quality color coding
- **Voigt Analysis Tab**: Renders 1D cross-sections or 2D PS2D contour dashboard
- **Spectrum Browser**: Tabbed viewer (Peak Table, Voigt Analysis, Metadata)

**Layout**:
```
┌────────────────────────────────────────────────────┐
│ Menu Bar                                           │
├───────────────┬────────────────────┬───────────────┤
│ Left Panel    │ Center Panel       │ Right Panel   │
│ (Controls)    │ (Spectrum + Plot)  │ (Navigator)   │
│               │                    │               │
│ Detection     │ 2D Contours        │ Peak List     │
│ Fitting       │ + Peak Markers     │ Quality       │
│ Parameters    │ + Ellipses         │ Navigation    │
├───────────────┴────────────────────┴───────────────┤
│ Status Bar: Project | Nucleus | Quality | Progress │
└────────────────────────────────────────────────────┘
```

---

## 6. Workflow Processors

| Processor | Purpose | File |
|-----------|---------|------|
| `SingleSpectrumProcessor` | Interactive analysis of one spectrum | `processors/single_spectrum_processor.py` |
| `MultiSpectrumProcessor` | Titration/condition series | `processors/multi_spectrum_processor.py` |
| `BatchProcessor` | CLI batch processing of folders | `batch_processing/batch_processor.py` |

All processors:
- Share identical core engines (PS2D, detection, overlap clustering)
- Use same parameter manager
- Produce same result schema
- **Guarantee parity** between GUI and CLI workflows

---

## 7. Parameter Flow

```
GUI Spinboxes (centroid_window_x_ppm = 0.01)
  ↓
ParameterManager.update_from_gui_variables()
  ↓
current_params['centroid_window_x_ppm'] = 0.01
  ↓
ParameterManager.get_effective_parameters()
  ↓
if simplified_mode:
    simplified_manager.get_simplified_parameters()  # calculates defaults
    Override with GUI values:                        # CRITICAL
        params['gui_params']['centroid_window_x_ppm'] = current_params[...]
  ↓
integrator.gui_params = params['gui_params']
  ↓
_detect_peaks_by_threshold() uses gui_params.get('centroid_window_x_ppm')
```

**Result**: GUI values **always override** calculated defaults in simplified mode.

---

## 8. Result Schema

All fitting engines return consistent dictionaries:

```python
{
  'assignment': str,
  'pos_f1': float,             # F1 position (15N/13C)
  'pos_f2': float,             # F2 position (1H)
  'lw_lor_f1': float,          # Lorentzian FWHM F1
  'lw_gau_f1': float,          # Gaussian FWHM F1
  'lw_lor_f2': float,          # Lorentzian FWHM F2
  'lw_gau_f2': float,          # Gaussian FWHM F2
  'intensity': float,          # Volume
  'height': float,             # Peak amplitude
  'r_squared': float,          # Quality metric
  'quality': str,              # 'Excellent'|'Good'|'Fair'|'Poor'|'Failed'
  'method': str,               # '2d_simultaneous_multi_peak'|'consensus'|'1d_staged'
  'fitted': bool,              # True if fitting succeeded
  'failure_reason': str        # (if fitted=False)
}
```

**Compatibility**: GUI tables, exports, ML collectors all expect this schema.

---

## 9. Extension Points

### Add New Fitting Engine
1. Create module in `core/` implementing `fit_peak()` method
2. Return result dictionary matching schema (section 8)
3. Register in `core_integrator.py` workflow

### Add New Processor
1. Subclass from `processors/base_processor.py` (or similar pattern)
2. Reuse `core_integrator` API
3. Add GUI entry point in `main_gui.py`

### Add ML Collector
1. Create collector in `ml/`
2. Register callback with PS2D fitter or batch processor
3. Export training samples to `ml_training_data/`

---

## 10. File Format Handlers

**Location**: `utils/file_manager.py`

**Supported**:
- Bruker TopSpin (`.2ii`, `.2rr`)
- NMRPipe (`.ft`, `.pipe`)
- Varian/Agilent (`.fid`, `.ft`)
- SPARKY (`.ucsf`)

**Loader**: `nmrglue` with automatic axis ordering detection.

---

## 11. Quality Assurance

### Failed Peak Handling
Peaks that fail fitting are stored with placeholder values (`quality='Failed'`, `fitted=False`) to maintain 1:1 index mapping between `peak_list` and `fitted_results`.

**Why**: Peak Navigator relies on index-based lookups. Dropping failed peaks breaks clicking navigation.

### Result Consolidation (Parallel Mode)
Uses integer `peak_number` matching instead of float coordinates to avoid precision errors.

**Old** (broken): `results_cache[(7.9158, 116.988)]` → lookup fails due to precision
**New** (fixed): `results_by_number[26]` → exact match

---

## 12. Performance Characteristics

| Component | Complexity | Bottleneck |
|-----------|------------|------------|
| Peak detection | O(N) | scipy.ndimage.maximum_filter |
| Centroid | O(W×H) per peak | Window size |
| Overlap clustering | O(P²) worst-case | Number of peaks P |
| PS2D fitting | O(I×n×M) | Iterations I, peaks n, masked points M |
| Parallel distribution | O(C) | Number of clusters C |

**Numba acceleration**: 3-5× speedup for PS2D Voigt primitive calculations.

---

## 13. Dependency Graph

```
launch_lunaNMR.py
  ├── gui/main_gui.py
  │     ├── gui/gui_components.py (Peak Navigator)
  │     ├── gui/visualization.py (plotting)
  │     └── gui/spectrum_browser.py (Voigt Analysis)
  │
  ├── processors/single_spectrum_processor.py
  │     ├── core/core_integrator.py
  │     │     ├── core/enhanced_peak_picker.py
  │     │     ├── core/ps2d_2d_fitter.py
  │     │     └── core/parallel_voigt_processor.py
  │     └── utils/parameter_manager.py
  │
  └── utils/config_manager.py
```

---

## 14. Thread Safety

**GUI**: Single-threaded (tkinter main loop)
**Parallel processing**: Multiprocessing (not threading) via `multiprocessing.Pool`
**State management**: Each worker gets copy of integrator state (no shared memory)

---

## 15. Testing Strategy

**Unit tests**: Core engines (detection, fitting, clustering)
**Integration tests**: Full workflows (load → detect → fit → export)
**Regression tests**: Compare sequential vs parallel results (must be identical)

**Test data**: `test_data/` folder (if present) or synthetic Voigt peaks.

---

## Maintenance Guidelines

1. **Parity**: Keep GUI and batch processors synchronized (same engines, parameters, results)
2. **Schema**: Extend result dictionary (section 8) when adding new metrics
3. **Documentation**: Update this file when adding/removing directories or major components
4. **Compatibility**: Maintain backward compatibility for result exports (CSV, JSON formats)

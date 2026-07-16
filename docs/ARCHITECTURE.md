# Architecture Guide

**TL;DR**: GUI (PySide6/Qt6) → processors (workflow orchestration) → core engines (detection, PS2D 2D fitting, overlap clustering). A headless CLI (`python -m lunaNMR`) shares the same engines. Parallel mode distributes whole overlap clusters across worker processes. GUI and CLI produce an identical result schema.

---

## 1. Directory Structure

```
lunaNMR_v1o0/
├── launch_lunaNMR.py       # Launcher (tkinter selector → Qt app)
├── docs/
├── lunaNMR/                # Main package
│   ├── cli.py              # Headless CLI dispatcher (python -m lunaNMR)
│   ├── __main__.py         # → cli.main()
│   ├── gui/                # PySide6/Qt6: main_window.py, base/, components/, dialogs/, styles/
│   ├── processors/         # Single/multi-spectrum workflow orchestration
│   ├── core/               # Detection, Voigt/PS2D fitting, overlap clustering, parallel processor
│   ├── batch_processing/   # CLI folder batch automation
│   ├── ml/                 # ML/statistics parameter prediction + collectors + gui/
│   ├── integrators/        # In-place advanced/series integrators
│   ├── utils/              # Config, params, file I/O, project bundles, output routing
│   └── validation/         # Installation + memory checks
└── modules/
    └── dynamiXs_v2o0/      # DynamiXs (relaxation/dynamics) + dynamiXs_Kd/ (Kd/titration)
```

tkinter is used **only** for the launcher's application-selector dialog; the app is entirely Qt.

---

## 2. Layer Responsibilities

| Layer | Modules | Role |
|-------|---------|------|
| **GUI** | `gui/` | Qt6 UI, plots, Peak Navigator, Voigt Analysis tab, dialogs |
| **Workflow** | `processors/` | Drive detection + fitting for single/series; batch in `batch_processing/` |
| **Core** | `core/` | Peak detection, PS2D 2D fitting, overlap clustering, parallel dispatch |
| **Utils** | `utils/` | Parameter manager, config, file formats, project bundles |
| **ML** | `ml/` | ML/stats parameter prediction; PS2D/peak training collectors |

---

## 3. Core Engines

### 3.1 Peak Detection (`core/enhanced_peak_picker.py`, `core_integrator.py`)
- Prominence/S-N detection via `scipy.ndimage.maximum_filter`
- Top-contour centroid refinement: `core_integrator._calculate_top_contour_center()`

### 3.2 PS2D Stack
- **Data selector** (`ps2d_data_selector.py`): union of elliptical masks around a cluster
- **Voigt primitives** (`ps2d_style_fitter.py`): Faddeeva (`scipy.special.wofz`), Numba-accelerated (`forceobj`)
- **2D fitter** (`ps2d_2d_fitter.py::fit_multi_peak_2d`): 5-stage Levenberg-Marquardt
- **Overlap clustering** (`core_integrator.py::identify_overlap_clusters`): agglomerative closest-pair merging under diameter/size constraints → disjoint clusters

### 3.3 Parameter Management (`utils/parameter_manager.py`)
- Simplified mode maps a few sliders onto the full legacy parameter set (`simplified_parameter_manager.py`)
- GUI spinbox values always override calculated defaults (see §7)

---

## 4. Data Flow

### 4.1 Sequential (default)
```
GUI → processors.SingleSpectrumProcessor
  → core_integrator.detect_peaks_sn_native()  → _detect_peaks_by_threshold() → _calculate_top_contour_center()
  → core_integrator.identify_overlap_clusters()          # agglomerative closest-pair clustering
  → per cluster:  single peak → enhanced_peak_fitting() (1D)
                  multi-peak  → ps2d_2d_fitter.fit_multi_peak_2d() (2D)
  → GUI tables/plots
```

### 4.2 Parallel (`core/parallel_voigt_processor.py::ParallelVoigtProcessor`)
```
identify_overlap_clusters()  (ONCE, disjoint clusters)
  → fit_all_peaks_parallel()
      → _create_cluster_tasks()  → _execute_parallel_cluster_fitting()  (multiprocessing.Pool, whole clusters per worker)
      → _consolidate_cluster_results()  (integer peak_number matching)
```

**Key differences**: identical clustering; workers receive whole clusters (not peaks); consolidation keys on integer `peak_number` (not float coords); ~2.7× speedup on 6 cores.

---

## 5. GUI Architecture

- **Model**: integrator state (`core_integrator.py`), peak lists, fit results
- **View**: `gui/main_window.py` (`LunaNMRMainWindow(BaseWindow)`), `gui/components/peak_navigator.py`, `spectrum_plotter.py`, `voigt_analysis_plotter.py`, `gui/dialogs/`
- **Controller**: GUI callbacks → processors → update model → refresh view

Dialogs subclass `BaseDialog(QDialog)`; stylesheets from `gui/styles/main.qss`. Layout: left controls | center spectrum+contours+markers | right Peak Navigator (quality-colored); status bar shows project/nucleus/quality/progress.

---

## 6. Workflow Processors

| Processor | Purpose | File |
|-----------|---------|------|
| `SingleSpectrumProcessor` | One spectrum, interactive | `processors/single_spectrum_processor.py` |
| `MultiSpectrumProcessor` | Titration/condition series | `processors/multi_spectrum_processor.py` |
| `BatchProcessor` | CLI folder batch | `batch_processing/batch_processor.py` |

All share the same core engines, parameter manager, and result schema — GUI/CLI parity.

## 6.1 CLI, Modules & Project Bundles

- **CLI** (`lunaNMR/cli.py`, `python -m lunaNMR`): wraps the same processors/backends, no Qt (lazy `__init__`, matplotlib forced to `Agg`). Subcommands: `series`, `dynamixs {t1t2,methyl-t2,hetnoe,density,modelfree}`, `kd`, `export kd`, `project {inventory,remove}`, `batch`. See [CLI.md](CLI.md).
- **Modules** (`modules/dynamiXs_v2o0/`): DynamiXs relaxation/dynamics + `dynamiXs_Kd/` Kd/titration; launched from the GUI Modules menu, mirrored by the CLI.
- **Project bundles** (`utils/project_manager.py`, `.lunaNMR`): state-based save/load. Fit surfaces stored slim, reconstructed on load (after the spectrum is loaded); DynamiXs and Kd persist as multiple named analyses.

---

## 7. Parameter Flow

```
GUI spinboxes → ParameterManager.update_from_gui_variables() → current_params
  → get_effective_parameters()
      simplified mode: compute defaults, then override with GUI values (CRITICAL)
  → integrator.gui_params → _detect_peaks_by_threshold() reads gui_params
```
GUI values **always override** calculated defaults in simplified mode.

---

## 8. Result Schema

All fitting engines return dicts with a consistent shape (GUI tables, exports, ML collectors depend on it):

```python
{
  'assignment': str, 'pos_f1': float, 'pos_f2': float,     # F1=15N/13C, F2=1H
  'lw_lor_f1': float, 'lw_gau_f1': float,                  # Lorentzian/Gaussian FWHM F1
  'lw_lor_f2': float, 'lw_gau_f2': float,                  # ... F2
  'intensity': float, 'height': float, 'r_squared': float,
  'quality': str, 'method': str, 'fitted': bool, 'failure_reason': str
}
```
`method` example: `'2d_simultaneous_multi_peak'`. Failed peaks are kept in place (`fitted=False`) to preserve 1:1 index mapping with `peak_list` — the Peak Navigator uses index-based lookups.

---

## 9. File Format Handlers

`utils/file_manager.py`, loaded via `nmrglue` with automatic axis-ordering detection:
Bruker TopSpin (`.2ii`/`.2rr`), NMRPipe (`.ft`/`.pipe`), Varian/Agilent (`.fid`/`.ft`), SPARKY (`.ucsf`).

---

## 10. Threads & Performance

- **GUI**: Qt event loop; `QThread` workers for non-blocking detection/fitting.
- **Parallel fitting**: `multiprocessing.Pool` (process-level; each worker gets a copy of integrator state, no shared mutable memory).
- **Numba**: 3-5× on PS2D Voigt primitives.
- Hot spots: detection O(N) (`maximum_filter`), centroid O(W×H)/peak, clustering O(P²) worst-case, PS2D O(I·n·M).

---

## Maintenance

1. Keep GUI and batch/CLI on the same engines, params, and result schema (parity).
2. Extend the §8 result dict when adding metrics; keep CSV/JSON exports compatible.
3. Update this file when adding/removing directories or major components.

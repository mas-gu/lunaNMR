# Architecture Guide

**TL;DR**: GUI (PySide6/Qt6) → processors (workflow orchestration) → core engines (detection, PS2D 2D fitting, overlap clustering). A headless CLI (`python -m lunaNMR`) shares the same engines. Sequential and parallel are one pipeline: parallel distributes whole overlap clusters across worker processes, sequential runs the same cluster tasks in process. GUI and CLI produce an identical result schema.

---

## 1. Directory Structure

```
lunaNMR_v1o0/
├── launch_lunaNMR.py       # Launcher: opens the Qt main window
├── docs/
├── lunaNMR/                # Main package
│   ├── cli.py              # Headless CLI dispatcher (python -m lunaNMR)
│   ├── __main__.py         # → cli.main()
│   ├── gui/                # PySide6/Qt6: main_window.py, base/, components/, dialogs/, styles/
│   ├── processors/         # Single/multi-spectrum workflow orchestration
│   ├── core/               # Detection, Voigt/PS2D fitting, overlap clustering, parallel processor
│   ├── batch_processing/   # CLI folder batch automation
│   ├── integrators/        # In-place advanced/series integrators
│   ├── utils/              # Config, params, file I/O, project bundles, output routing
│   └── validation/         # Installation + memory checks, `diagnose` pre-flight (spectra_check.py)
└── modules/
    ├── dynamiXs_v2o0/      # DynamiXs (relaxation/dynamics) + dynamiXs_Kd/ (Kd/titration)
    └── integration_1d_v1o0/ # 1D peak integration across a series
```

tkinter is used **only** for the launcher's startup error dialogs; the app is entirely Qt.

---

## 2. Layer Responsibilities

| Layer | Modules | Role |
|-------|---------|------|
| **GUI** | `gui/` | Qt6 UI, plots, Peak Navigator, Voigt Analysis tab, dialogs |
| **Workflow** | `processors/` | Drive detection + fitting for single/series; batch in `batch_processing/` |
| **Core** | `core/` | Peak detection, PS2D 2D fitting, overlap clustering, parallel dispatch |
| **Utils** | `utils/` | Parameter manager, config, file formats, project bundles |

---

## 3. Core Engines

### 3.1 Peak Detection (`core/enhanced_peak_picker.py`, `core_integrator.py`)
- Prominence/S-N detection via `scipy.ndimage.maximum_filter`, per sign: maxima above +threshold and minima below −threshold (negative peaks are kept)
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

### 4.1 Detection
```
GUI → core_integrator.detect_peaks_sn_native()  → _detect_peaks_by_threshold() → _calculate_top_contour_center()
```

### 4.2 Fitting: one pipeline, sequential or parallel
```
GUI → processors.SingleSpectrumProcessor.process_peak_list(stop=)
  → restore_base_ps2d_config()        # fresh run: no locked clusters / learned statistics
  → enhanced_voigt_fitter.enhanced_peak_fitting_parallel(in_process=, stop=)
  → parallel_voigt_processor.ParallelVoigtProcessor.fit_all_peaks_parallel()
      → core_integrator.identify_overlap_clusters()   (ONCE, disjoint clusters)
      → PASS 1: isolated peaks → linewidth statistics
      → adaptive radius search (only when use_adaptive_optimization; off by default)
      → PASS 2: multi-peak clusters with the learned statistics
      each pass: _create_cluster_tasks() → _execute_parallel_cluster_fitting()
                 → _parallel_cluster_worker(task) per cluster:
                     single peak → enhanced_peak_fitting() (1D)
                     multi-peak  → fit_overlap_group_2d() → ps2d_2d_fitter.fit_multi_peak_2d() (2D)
      → _consolidate_cluster_results()  (integer peak_number matching)
  → GUI tables/plots
```

**Sequential vs parallel**: with `in_process=True` (sequential, or fewer than 3 peaks) each cluster task is a direct call to `_parallel_cluster_worker`; otherwise the tasks go to a `multiprocessing` pool from a `spawn` context. Same clustering, passes and locking, so the two give identical results; parallel pays off only on large peak lists (~2.7× on 6 cores). Workers rebuild the integrator from a serialised copy (`_serialize_integrator_params`, `_prepare_shared_context`).

**Cancel** is a `threading.Event` passed down as `stop=` (`utils/cancellation.py`: `Cancelled`, `raise_if_stopped`). It is checked before each cluster, each detection peak and each adaptive-search fit, and polled while waiting on pool tasks; a stopped run raises `Cancelled` and returns nothing — no fallback, no partial results.

---

## 5. GUI Architecture

- **Model**: integrator state (`core_integrator.py`), peak lists, fit results
- **View**: `gui/main_window.py` (`LunaNMRMainWindow(BaseWindow)`), `gui/components/peak_navigator.py`, `spectrum_plotter.py`, `voigt_analysis_plotter.py`, `gui/dialogs/`
- **Controller**: GUI callbacks → processors → update model → refresh view

Dialogs subclass `BaseDialog(QDialog)`; stylesheets from `gui/styles/main.qss`. Layout: left controls | center spectrum+contours+markers | right Peak Navigator (rows colored by R²); the status bar shows status messages (workflow mode, zoom, errors).

---

## 6. Workflow Processors

| Processor | Purpose | File |
|-----------|---------|------|
| `SingleSpectrumProcessor` | One spectrum, interactive | `processors/single_spectrum_processor.py` |
| `MultiSpectrumProcessor` | Titration/condition series | `processors/multi_spectrum_processor.py` |
| `BatchProcessor` | CLI folder batch | `batch_processing/batch_processor.py` |

All share the same core engines, parameter manager, and result schema — GUI/CLI parity.

## 6.1 CLI, Modules & Project Bundles

- **CLI** (`lunaNMR/cli.py`, `python -m lunaNMR`): wraps the same processors/backends, no Qt (lazy `__init__`, matplotlib forced to `Agg`). Subcommands: `diagnose`, `series`, `dynamixs {t1t2,methyl-t2,hetnoe,t1rho,density,modelfree}`, `kd`, `export kd`, `integrate-1d`, `peaks shift`, `project {inventory,export,remove,build}`, `batch`. The CLI passes no stop event. See [CLI.md](CLI.md).
- **Modules** (`modules/dynamiXs_v2o0/`, `modules/integration_1d_v1o0/`): DynamiXs relaxation/dynamics, `dynamiXs_Kd/` Kd/titration and 1D integration, each its own entry of the GUI Modules menu, mirrored by the CLI.
- **Project bundles** (`utils/project_manager.py`, `.lunaNMR`): state-based save/load. Fit surfaces stored slim, reconstructed on load (after the spectrum is loaded); each series carries its `batch_results.json`, tidy CSV, intensity and volume matrices, tracking CSV and `series_metadata.json`, and a series opened from a project reads and saves only those; DynamiXs and Kd persist as multiple named analyses. Opening a project starts a clean session.

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

All fitting engines return dicts with a consistent shape (GUI tables and exports depend on it):

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
- **Parallel fitting**: `multiprocessing` pool from a `spawn` context (process-level; each worker rebuilds the integrator from a serialised copy, no shared mutable memory).
- **Numba**: 3-5× on PS2D Voigt primitives.
- Hot spots: detection O(N) (`maximum_filter`), centroid O(W×H)/peak, clustering O(P²) worst-case, PS2D O(I·n·M).

---

## Maintenance

1. Keep GUI and batch/CLI on the same engines, params, and result schema (parity).
2. Extend the §8 result dict when adding metrics; keep CSV/JSON exports compatible.
3. Update this file when adding/removing directories or major components.

# Architecture Guide – lunaNMR v0.9

This document summarizes the main subsystems in `lunaNMR_v0o9`, how they interact, and where new features such as the PS2D 2D fitting pipeline and automated ML hooks live.

## 1. High-Level Layout

```
lunaNMR_v0o9/
├── launch_lunaNMR.py           # Entry-point launcher
├── docs/                       # Contributor & user documentation
├── data_example/               # Reference spectra
├── lunaNMR/                    # Primary Python package
│   ├── gui/                    # Tkinter GUI and visualization
│   ├── processors/             # Workflow orchestration (single, batch, series)
│   ├── core/                   # Numerical engines: detection, PS2D fitting, overlap resolver
│   ├── batch_processing/       # CLI batch runner + ML data capture
│   ├── integrators/            # Specialized in-memory integrators
│   ├── utils/                  # Config, parameter, and file managers
│   ├── ml/                     # Training data collectors for PS2D
│   └── validation/             # Installation diagnostics
├── modules/                    # Optional add-ons (e.g., DynamiXs)
└── batch_ml/, ml_training_data/ # Output and training artefacts
```

## 2. Layered Responsibilities

| Layer | Modules | Responsibilities |
|-------|---------|------------------|
| **Presentation** | `lunaNMR/gui` | Tkinter UI, interactive plots, peak navigator, spectrum browser, Voigt analysis dashboard |
| **Workflow** | `lunaNMR/processors`, `lunaNMR/batch_processing` | Turn detection + fitting engines into end-user flows (single spectrum, multi-condition, batch, CLI) |
| **Core Services** | `lunaNMR/core` | Peak detection, PS2D 2D multi-peak fitting, overlap grouping, consensus fitting, jackknife validation, correlation analysis |
| **Domain Utilities** | `lunaNMR/utils`, `lunaNMR/integrators` | Parameter management, config storage, file I/O, nucleus-specific helpers, specialized integrators |
| **Machine Learning** | `lunaNMR/ml`, `batch_ml/`, `ml_training_data/` | Capture PS2D fitting diagnostics, export training samples, notebooks and artefacts |

Processors call into the core engines, then hand results to the GUI or batch exporters. The batch processor shares the same engines, guaranteeing parity between interactive and scripted workflows.

## 3. Core Engine Highlights

- **Enhanced Peak Picker** (`core/enhanced_peak_picker.py`): multi-method detection with clustering, S/N analysis, and adaptive windows.
- **PS2D Fitting Stack**:
  - `ps2d_data_selector.py` – constructs elliptical masks around overlap clusters.
  - `ps2d_style_fitter.py` – Faddeeva-based Voigt primitives with optional Numba acceleration.
  - `ps2d_2d_fitter.py` – five-stage Levenberg–Marquardt optimizer for simultaneous 2D fits.
  - `ps2d_exact_overlap_detector.py` / `ps2d_exact_overlap_integration.py` – geometric overlap grouping and integration wrappers.
  - `overlap_resolver_engine.py` – staged 1D backup pipeline with model-selection, jackknife, and correlation analysis.
- **Consensus & Simplified Modes** (`core/consensus_fitting_engine.py`, `utils/simplified_parameter_manager.py`): reduce the parameter space to a handful of intuitive controls.
- **Validation & QA**: `jackknife_validator.py`, `parameter_correlation_analyzer.py`, and PS2D training collector feed both user diagnostics and ML pipelines.

## 4. Workflow Orchestration

- **Single Spectrum Processor** orchestrates detection, PS2D routing, and GUI table population.
- **Multi Spectrum Processor** iterates over conditions (temperature, titration, etc.) while sharing parameter presets.
- **Series Processor** handles kinetic/relaxation series and exports time-course plots.
- **Batch Processor** (CLI or GUI-triggered) applies the exact same PS2D and consensus engines, logging rich metadata and optionally persisting ML samples.

## 5. GUI Architecture

The GUI follows an MVC-inspired split:

- **Model**: Integrator state, fitted peak lists, batch results.
- **View**: `main_gui.py` dashboards, `gui_components.py` widgets, `visualization.py` plotting helpers, `spectrum_browser.py` interactive viewer.
- **Controller**: Callbacks that forward requests to processors and update the shared integrator.

The “Voigt Analysis” tab renders either 1D cross sections or the PS2D 2×2 contour dashboard depending on the fitting method returned.

## 6. Extension Points

- **New processors**: subclass `BaseProcessor`‑style utilities in `processors/` and reuse the integrator API.
- **Additional fitting strategies**: integrate with `core/overlap_resolver_engine.py` or provide new engines that conform to the same result schema (`pos_f1`, `pos_f2`, widths, intensities, quality metrics).
- **Machine learning**: drop collectors into `lunaNMR/ml/` and register them with the PS2D or batch processors.
- **Optional modules**: place self-contained apps inside `modules/` and expose them through `launch_lunaNMR.py`.

## 7. Runtime Data Flow (PS2D path)

```
GUI action → processors.single_spectrum_processor
          → core.integrated_detection_fitter (detection + overlap analysis)
          → core.ps2d_data_selector (mask)
          → core.ps2d_2d_fitter (stage 0–4 fit)
          → core.ps2d_training_collector (optional ML sample)
          → processors format results → GUI tables / plots
```

Fallback behaviour swaps in the staged 1D engine (`staged_fitting_strategy.py`) via the overlap resolver when PS2D is disabled or the data is purely 1D.

## 8. Keeping the Architecture Healthy

- Maintain parity between GUI and batch processors so that regression tests cover both code paths.
- When introducing new engines, extend the shared result schema (`pos_f1`, `pos_f2`, `lw_lor_*`, `lw_gau_*`, `intensity`, quality flags) to keep downstream consumers compatible.
- Update this document and `PACKAGE_STRUCTURE.md` whenever directories are added, renamed, or deprecated.

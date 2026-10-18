# Multi-Spectrum Series Guide — LunaNMR v1.0

Process a set of related spectra (titration, relaxation, temperature, time-series) with one
reference peak list: fit the reference peaks in every spectrum and export per-point
intensities/volumes/positions as CSV.

---

## Entry points

| Where | Effect |
|-------|--------|
| Main window **Fit Series** button (`on_fit_series`) | Opens `SeriesIntegrationDialog` |
| Results ▸ **Multi-Spectrum Overlay Viewer…** | Opens `MultiSpectrumViewerDialog` (visualization) |
| CLI | `python -m lunaNMR series …` (see `CLI.md`) |

The dialog runs processing in a `ProcessingWorker` (a QObject moved to a QThread); the
overlay viewer is independent and reads finished results.

---

## Configure the run (`SeriesIntegrationDialog`)

1. **Load a reference peak list** (main window) — its positions/assignments are the peaks
   tracked across the whole series.
2. **Pick the spectra folder** — all compatible NMR files are collected and natural-sorted.
3. Choose **Peak Source** and **X-axis mode** (below), then **▶ Start Integration**. Parallel
   processing follows the main window's Expert Mode **Use Parallel Processing** box; sequential
   mode runs the same pipeline in process and gives an identical series.

The series never runs the adaptive radius search (the dialog logs "Adaptive radius search:
off for a series"); the Expert Mode box applies to Fit Spectrum only. Each run starts from
the base PS2D config (nucleus defaults plus the Expert panel's Apply). Spectrum 1's
detection adds nearby unmatched peaks (dummies) unless Expert Mode **Auto-add nearby
unmatched peaks** is unticked; they are fitted at every point.

**✕ Cancel** stops the run before its next cluster. A cancelled series writes nothing — no
`series_results_*` folder — and says so in the log.

### Peak Source (radio group, default **Cascade**)

| Mode | Behavior |
|------|----------|
| **Cascade** | Start from reference/detected peaks; refine positions per spectrum, carrying spectrum N-1's positions forward (only mode that propagates drift). |
| **Detected** | Use the peaks already detected+fitted in the current spectrum ("Fit All Peaks"). |
| **Independent** | Full detect+fit per spectrum from the *original* reference (no position drift). Same pipeline as GUI "Fit Spectrum". Refused for a titration. |

CLI `--peak-source` additionally accepts `reference` (transfer reference positions verbatim).

### X-axis mode (combo: Off / Time / Titration)

Extracts a per-spectrum x value from filenames and uses it for column headers / DynamiXs.

| Mode | `series_mode` | Filename tokens |
|------|---------------|-----------------|
| Time | `time` | delays, e.g. `_50ms`, `_1s` → headers `50`, `1000` (ms); duplicates get `_2` |
| Titration | `titration` | points, e.g. `_1o0`, `_0o5` |

**Titration override** (`resolve_titration_tracking`): titration mode forces **cascade**
tracking and disables the absolute drift limit — cascade is the only mode that carries a
peak's position spectrum-to-spectrum as the peak moves during binding. Independent is
refused outright (dialog and CLI). Each peak also gets its own position window (its
territory: 0.49 × the distance to its nearest neighbour), so no peak can reach another.

---

## Fix Positions / Fix Linewidths

Read from main-window state into `gui_params` (default **OFF**). They gate PS2D L-M stages:

- **Fix Positions** → skips Stage 2 (position refinement) and holds positions in the final
  stage; positions stay at their seed (the reference, or the previous point in cascade).
- **Fix Linewidths** → holds the linewidths in every stage; Stage 1 still runs and fits
  the intensities, so only intensities vary.

Guidance: relaxation/titration usually leave positions to the chosen peak-source mode and
linewidths free; fix linewidths only when they are known to be constant.

---

## Output files

Written to `series_results_<name>_<timestamp>/` (`series_results_<timestamp>/` without a series
name) in the folder that holds the first spectrum file.

| File | Shape | Contents |
|------|-------|----------|
| `comprehensive_peak_tracking.csv` | wide, one row/reference peak | `{pt}_Position_X/Y`, `{pt}_Height`, `{pt}_Volume` per point |
| `series_analysis_tidy.csv` | long, one row/(spectrum×peak) | `spectrum_name, assignment, peak_number, ppm_x, ppm_y, height, volume, snr, quality, r_squared, detection_sigma, position_pinned, intensity_pinned` |
| `peak_intensity_matrix.csv` | matrix, rows=peaks | fitted heights per point (an intensity-pinned or Failed point is empty) |
| `peak_volume_matrix.csv` | matrix | volumes per point (same blanking) |
| `peak_intensity_detected_matrix.csv` | matrix | detected (raw) heights per point |
| `peak_detection_matrix.csv` | matrix | 1 where the peak was detected at that point (3σ or more), else 0 |
| `series_metadata.json` | JSON | each matrix column's source spectrum and parsed value, `value_units`, `repeat_scale`, the series' base `ps2d_config` and each point's `territories` |

Fitted **positions are formatted to 4 decimals** (`ppm_x/ppm_y` in tidy;
`_format_dataframe_for_csv` for tracking) — CSP/titration analysis needs the precision.

**Downstream consumers:** DynamiXs (intensity or volume matrix, per its **Fit on: Height /
Volume** choice; delay units from `series_metadata.json`), Kd module (tidy + comprehensive),
overlay viewer.

---

## Overlay viewer (`MultiSpectrumViewerDialog`)

- Overlays contours from any subset of loaded spectra (per-spectrum color, toggle visibility).
- Per-spectrum fitted-peak markers, each drawn in that spectrum's color (so overlaid
  lists stay distinguishable). A global **Show Peak Markers** checkbox gates all of them;
  a per-spectrum "Peak lists" panel toggles each. Assignment labels show only for the
  active spectrum.
- **Edit positions** checkbox: two-click peak move (select nearest, then click true
  position); middle-click nudges without edit mode. The refit uses the series' recorded PS2D
  config (from `series_metadata.json`) with its own position window, not the territories.
  **Save corrected** finds each point's columns through `series_metadata.json` (a series
  without it is refused) and writes the refit's quality, pinned flags and `detection_sigma`.
  On a series opened from a project it writes into the bundle copy.

---

## Troubleshooting

| Symptom | Likely cause / fix |
|---------|--------------------|
| Peaks drift when they shouldn't | Enable **Fix Positions**, or use Independent mode (no drift). |
| Peaks fail to follow a moving titration peak | Use Titration x-axis mode (forces cascade; per-step window from the reference linewidths, bounded by each peak's territory). |
| Peak missing in later spectra | Intensity fell below noise, or real (exchange/binding); check raw contour in overlay. |
| Inconsistent R² across series | Lower S/N or changed peak shape in that spectrum; inspect in overlay/navigator. |
| Parallel not faster | Overhead on small sets; use for many peaks/spectra. |

---

## Code references

| Component | File | Key symbols |
|-----------|------|-------------|
| Config dialog + worker | `gui/dialogs/series_integration_dialog.py` | `SeriesIntegrationDialog`, `ProcessingWorker` |
| Series engine | `processors/multi_spectrum_processor.py` | `MultiSpectrumProcessor.process_nmr_series`, `resolve_titration_tracking`, `_create_comprehensive_output_files`, `_create_tidy_results_file` |
| Overlay viewer | `gui/dialogs/multi_spectrum_viewer_dialog.py` | `MultiSpectrumViewerDialog` |
| CLI | `lunaNMR/cli.py` | `series` subcommand |

**See also:** `CLI.md`, `ALGORITHMS.md`, `ARCHITECTURE.md`, `lunaNMR_guide.html`.

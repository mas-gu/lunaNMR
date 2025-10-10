# GUI User Guide – lunaNMR v0.9

This guide walks through the graphical application that ships with `lunaNMR_v0o9`, highlighting PS2D fitting controls, overlap diagnostics, batch entry points, and the Voigt analysis dashboards.

## 1. Launching the Application

```bash
cd lunaNMR_v0o9
python3 launch_lunaNMR.py
```

Choose **LunaNMR** in the launcher dialog. The optional **DynamiXs** module appears if present in `modules/`.

## 2. Main Window Layout

```
┌─────────────────────────────────────────────────────────────┐
│ Menu bar (File, Analysis, View, Tools, Help)                │
├─────────────────────────────────────────────────────────────┤
│ Shortcut tiles: Single Spectrum | Series Analysis | Browser │
│                       | Configuration | Batch ML            │
├─────────────────────────────────────────────────────────────┤
│ Status cards, recent files, and quick-start hints           │
├─────────────────────────────────────────────────────────────┤
│ Status bar: current project, PS2D state, fitting quality    │
└─────────────────────────────────────────────────────────────┘
```

Key navigation tiles open dedicated workflows:
- **Single Spectrum Analysis** – interactive detection + PS2D fitting.
- **Series Analysis / Batch Processing** – multi-sample automation.
- **Spectrum Browser** – tabbed Voigt analysis, spectroscopy view.

## 3. Core Controls & Toggles

### Detection & Fitting Panel

Located on the right-hand side of the single spectrum workspace.

| Control | Purpose |
|---------|---------|
| `Detect Peaks` | Runs the enhanced peak picker. |
| `PS2D multi-peak` (toggle) | Enables the simultaneous 2D fitter for overlap clusters. Leave enabled for all 2D HSQC/HMQC work. |
| `PS2D linewidth reuse` | Reuses reference linewidths across a spectrum/series for increased stability. |
| `Simplified mode` | Activates the 3-parameter control surface (sensitivity, window scale, quality). |
| `Baseline Method` | Switch between automatic ArPLS, polynomial, or manual settings. |
| `Window X/Y` | Adjust fitting windows (affects both PS2D mask radii and staged 1D fallback bounds). |

PS2D status indicators in the status bar confirm when the multi-peak engine is active.

### Overlap Diagnostics

When PS2D is enabled, detected peaks are clustered by ellipsoidal overlap. The UI displays:
- Number of peaks in each cluster.
- Estimated separation in F1/F2.
- Whether the cluster was routed through PS2D or the staged 1D fallback.

## 4. Running a Single Spectrum Fit

1. **Load data** (File → Open Spectrum) and optionally a peak list (File → Open Peak List).
2. **Configure parameters** in the detection panel. For most HSQC datasets, leave PS2D toggled on and use simplified mode defaults.
3. **Detect peaks**. Review the table for assignments, detected coordinates, S/N.
4. **Fit selected peaks** using the “Fit Peaks” or context-menu actions. The GUI routes overlapping groups to PS2D automatically.
5. **Inspect results** in the results table and Voigt analysis tab (see below).

## 5. Voigt Analysis Dashboard

Open the Voigt Analysis tab (via toolbar or the Spectrum Browser). Two display modes exist:

- **1D mode** (when a staged fit is returned): stacked plots of ¹H and ¹⁵N cross sections, fitted curves, residuals, and window annotations.
- **2D PS2D mode**: 2×2 layout showing experimental contours (top left), individual fitted peaks (top right), residual heatmap (bottom left), and a parameter summary (bottom right). The peak markers in both top panels use the fitted coordinates (`pos_f2`, `pos_f1`) so you can compare centre shifts directly; the original detection coordinates remain in the results table for reference.

Tooltips and legends indicate R² values, window sizes, and colour coding per peak. Use the navigator controls to step through peaks or clusters.

## 6. Spectrum Browser Highlights

- Presents a tabbed interface (Peak Table, Voigt Analysis, Metadata) for quick inspections.
- Accepts the same result dictionaries as the main GUI; PS2D contour plots are available here as well.
- Handy for reviewing batch outputs: load a `.json` export and browse without re-running fits.

## 7. Batch & Series Workflows

- **Batch Processing** tile opens a wizard for folder selection, nucleus presets, S/N thresholds, and output destinations. Batch runs reuse the interactive engines (PS2D, consensus, simplified manager) ensuring parity.
- **Series Analysis** supports relaxation or titration series. Enable PS2D linewidth reuse if you need consistent linewidths across time points.

Progress indicators, log panes, and export buttons are shared across workflows for a consistent experience.

## 8. Tips & Troubleshooting

- If PS2D fails to converge, the system falls back to staged 1D fitting. Warnings appear in the log panel; you can retry with adjusted windows or temporarily disable PS2D for that peak.
- Use the simplified parameter mode when onboarding new users—only three sliders control the full parameter set.
- The status bar shows the active nucleus, current parameter preset, and the last fit quality. Clicking the status widgets opens the underlying configuration dialogs.
- Need raw numbers? Use `File → Export Results…` to write CSV/JSON files containing both detection coordinates and final fitted centres/widths.

Keeping the GUI in sync with the underlying processors ensures reproducible results between interactive and automated analyses.

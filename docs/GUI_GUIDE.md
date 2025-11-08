# GUI User Guide

**TL;DR**: Launch with `python3 launch_lunaNMR.py`. Main workflow: Load Spectrum → Detect Peaks → Fit All Peaks. PS2D 2D fitting is automatic for overlaps. Use simplified mode (3 sliders) for easy parameter control. Quality color coding: Green (Excellent/Good), Yellow (Fair), Red (Poor/Failed).

---

## 1. Launch
```bash
cd lunaNMR_v0o9
python3 launch_lunaNMR.py
```

Choose **LunaNMR**. (DynamiXs module appears if installed in `modules/`.)

---

## 2. Main Window Layout

```
┌─────────────────────────────────────────────────────┐
│ Menu Bar: File, Analysis, View, Tools, Help         │
├─────────────────────────────────────────────────────┤
│ Left:  Controls (detection, fitting, parameters)    │
│ Center: Spectrum plot (2D contours + peak markers)  │
│ Right:  Peak Navigator (browse peaks 1-N)           │
├─────────────────────────────────────────────────────┤
│ Status Bar: Project | Nucleus | Quality | Progress  │
└─────────────────────────────────────────────────────┘
```

---

## 3. Core Controls

### 3.1 Detection Parameters

| Control | Default | Purpose |
|---------|---------|---------|
| `1H/15N (ppm)` | 0.03 / 0.2 | Search window size |
| `Noise Threshold` | 0.5 | S/N multiplier for detection |
| `X×Y (pixels)` | 3 × 1 | Minimum peak separation |

### 3.2 Peak Centroid Detection

| Control | Default | Purpose |
|---------|---------|---------|
| `Centroid Window X (ppm)` | 0.01 | Max shift in F2 (1H) for top contour centroid |
| `Centroid Window Y (ppm)` | 0.10 | Max shift in F1 (15N/13C) for top contour centroid |

**Effect**: Peak positions refined to geometric center of ±5% intensity band. Console shows "Centroid shift: Δ=X.XXXX ppm from pixel max" for shifts >0.001 ppm.

**Recommended values**:
- Tight (default): 0.01 / 0.10 ppm — prevents jumping to wrong peaks
- Moderate: 0.02 / 0.20 ppm — allows more correction for flat-top peaks
- Loose: 0.04 / 0.40 ppm — maximum correction (use with caution)

### 3.3 Fitting Toggles

| Toggle | Purpose |
|--------|---------|
| `PS2D multi-peak` | Enables 2D simultaneous fitting for overlaps (default: ON) |
| `Fix Positions` | **Absolute constraint**: zero position drift during fitting |
| `Fix Linewidths` | Locks linewidths to reference values |
| `Use Parallel Processing` | Distributes clusters across 75% of CPU cores |
| `Simplified Mode` | 3-slider control (sensitivity, window scale, quality target) |

**Fix Positions behavior**:
- Stage 2 (position refinement) is **skipped entirely**
- Stage 4 locks positions to fixed_params dictionary
- After bounds clipping, positions are restored to exact values
- **Result**: Zero drift (within float precision ~10⁻¹⁵)
- **Use case**: Series integration with consistent peak positions across spectra

---

## 4. Workflow

### Single Spectrum
1. **File → Open Spectrum** (Bruker `.2ii`, NMRPipe `.ft`, SPARKY `.ucsf`)
2. **File → Open Peak List** (optional reference peaks)
3. Click **Detect Peaks** → wait for console completion
4. Click **Fit All Peaks** → PS2D routes overlaps to 2D simultaneous fitting
5. Review **Peak Navigator** → click peaks to inspect Voigt analysis

### Series Analysis
1. **Analysis → Start Series Integration**
2. Select folder containing multiple spectra
3. Configure first spectrum parameters (PS2D, fix positions, etc.)
4. Click **Process Series** → identical fitting applied to all spectra
5. Export results: `File → Export Series Results`

---

## 5. Peak Navigator

**Location**: Right panel in main window, right panel in Voigt Analysis tab (Spectrum Browser)

**Features**:
- List of all detected peaks with assignment, coordinates, height
- Quality color coding:
  - **Green**: Excellent (R²≥0.9) or Good (R²≥0.8)
  - **Yellow**: Fair (R²≥0.5)
  - **Red**: Poor (R²≥0.2) or Failed (R²<0.2)
- Navigation buttons: **◀ Previous** | **🔬 Analyze** | **Next ▶**
- Click peak → auto-centers spectrum and switches to Voigt Analysis tab

---

## 6. Voigt Analysis Tab

**Access**: Toolbar button or double-click peak in Navigator

### 1D Mode (staged fitting)
- Top: 1H cross-section (experimental + fitted)
- Middle: 15N cross-section
- Bottom: Residuals
- Annotations: R², window sizes, linewidths

### 2D Mode (PS2D multi-peak)
```
┌───────────────┬───────────────┐
│ Experimental  │ Fitted Peaks  │ ← Contour plots
│ Contours      │ (individual)  │
├───────────────┼───────────────┤
│ Residual      │ Peak          │ ← Residual heatmap
│ Heatmap       │ Navigator     │    + navigator
└───────────────┴───────────────┘
```

**Peak markers**: Use fitted coordinates (`pos_f2`, `pos_f1`) for direct comparison with experimental contours.

---

## 7. Quality Indicators

### Color Coding
| Quality | R² Range | Color | Meaning |
|---------|----------|-------|---------|
| Excellent | ≥ 0.9 | Green | High confidence |
| Good | [0.8, 0.9) | Green | Acceptable |
| Fair | [0.5, 0.8) | Yellow | Review recommended |
| Poor | [0.2, 0.5) | Red | Low confidence |
| Failed | < 0.2 | Red | Fitting failed |

### Console Messages
- `Centroid shift: Δ=0.0123 ppm from pixel max` — Position refined by centroid
- `✅ 125 peaks detected` — Detection complete
- `116/125 matched, 9 references retained` — Reference peak matching summary
- `Cluster 5: 3 peaks → PS2D 2D fitting` — Overlap resolved via PS2D

---

## 8. Overlap Visualization

**Enable**: Check **"Show Ellipses"** in Detection panel

**Ellipse Colors**:
- **Cyan dashed**: Data selector region (PS2D extracts data from this area)
- **Orange solid**: Overlap threshold (two circles touch = peaks will be clustered)
- **Magenta solid**: Fitting region (PS2D optimizes parameters within this window)

**Overlap threshold defaults** (two-circle touching test):
- **15N**: ±0.04 ppm (F2), ±0.4 ppm (F1)
- **13C**: ±0.04 ppm (F2), ±0.1 ppm (F1)

Adjust in **PS2D Configuration** menu if needed.

---

## 9. Simplified Mode

**Toggle**: `Simplified Mode` checkbox

**Controls** (3 sliders replace 25+ parameters):
- `Sensitivity` (0.1-0.9): Lower = more sensitive detection
- `Window Scale` (0.1-10.0): Fitting window size multiplier
- `Quality Target` (0.5-0.95): Target R² for acceptance

Maps to nucleus-specific adaptive thresholds internally. **GUI centroid parameters always override** calculated defaults.

---

## 10. Parallel Processing

**Enable**: Check **"Use Parallel Processing (75% cores)"**

**Behavior**:
- Calls `identify_overlap_clusters()` **once** (deterministic, same as sequential)
- Distributes clusters (not individual peaks) across workers
- Results **identical** to sequential mode
- Performance: 2.7× speedup with 6 cores

**Note**: All GUI parameters (fix positions, fix linewidths, centroid windows) are synchronized to workers.

---

## 11. Batch Processing

**Access**: **Tools → Batch Processing** menu

**Features**:
- Process folder of spectra with same parameters
- Supports all file formats (Bruker, NMRPipe, SPARKY, Varian)
- Parallel processing supported
- ML training data collection (optional)
- Export: CSV, JSON, NMRPipe peak lists

---

## 12. Keyboard Shortcuts

| Key | Action |
|-----|--------|
| `Ctrl+O` | Open Spectrum |
| `Ctrl+L` | Open Peak List |
| `Ctrl+D` | Detect Peaks |
| `Ctrl+F` | Fit All Peaks |
| `Ctrl+S` | Save Results |
| `Ctrl+E` | Export Results |
| `←/→` | Navigate Previous/Next Peak |
| `Enter` | Analyze Selected Peak |

---

## 13. Troubleshooting

**Problem**: PS2D fitting fails (R² < 0.2)
- **Solution**: Check overlap ellipses (may be too tight), increase window scale, or temporarily disable PS2D for that peak.

**Problem**: Peaks detected at wrong positions
- **Solution**: Tighten centroid windows (decrease X/Y ppm values) to prevent jumping.

**Problem**: Parallel mode slower than sequential
- **Solution**: Small datasets (<50 peaks) have overhead. Use parallel only for >100 peaks.

**Problem**: "Fix Positions" not working
- **Solution**: Ensure checkbox is checked before clicking "Fit All Peaks". Console should show Stage 2 skipped.

---

## 14. Export Formats

**CSV**: Peak table with assignment, coordinates, linewidths, intensities, quality
**JSON**: Full result dictionaries (including fitted surfaces, residuals)
**NMRPipe**: `.tab` peak list for use in NMRDraw/NMRView
**SPARKY**: `.list` assignment file

---

## 15. File Formats Supported

| Format | Extension | Software |
|--------|-----------|----------|
| Bruker TopSpin | `.2ii`, `.2rr` | TopSpin |
| NMRPipe | `.ft`, `.pipe` | NMRPipe |
| Varian/Agilent | `.fid`, `.ft` | VnmrJ |
| SPARKY | `.ucsf` | SPARKY |

Loaded via `nmrglue` with automatic axis ordering detection.

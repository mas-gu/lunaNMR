# Multi-Spectrum Analysis Guide - LunaNMR v0.9

**TL;DR**: Use **Analysis → Start Series Integration** for titrations, relaxation experiments, and time-series. Configure first spectrum, then apply identical fitting to all others with fixed peak positions. Results exported as CSV with one row per peak per spectrum.

---

## Overview

Multi-spectrum analysis processes datasets where:
- **Same sample** analyzed under different conditions (temperature, pH, ligand concentration)
- **Same peaks** expected across all spectra (positions consistent)
- **Different intensities/linewidths** due to experimental conditions

**Common Applications:**
- Titration experiments (ligand binding)
- Relaxation experiments (T1, T2, NOE)
- Temperature series
- Time-resolved experiments

---

## Workflow

### 1. Access Multi-Spectrum Viewer

**From Main GUI:**
```
Fit Series
```

**What Opens:**
- **Multi-Spectrum Viewer** window
- Left panel: Spectrum list with toggle buttons
- Right panel: Overlay plot with all spectra
- Bottom panel: Peak table and fitting controls

### 2. Load Spectra

**Load Reference Spectrum**
1. Click **Load Data**
2. Select directory containing spectra
3. All compatible files loaded automatically (.2ii, .ft, .ucsf, etc.)
4. Spectra sorted by filename (natural sort)


**Spectrum Naming:**
- Default: Filename without extension
- Edit by double-clicking spectrum in list
- Use meaningful names (e.g., "0mM", "5mM", "10mM" for titration)

### 3. Configure First Spectrum

**Critical: First spectrum defines reference parameters for entire series**

**Steps:**

1. **Select First Spectrum**
   - Click spectrum name in list
   - Spectrum displays in overlay plot

2. **Open Peak List** (optional but recommended)
   - Depending on analysis mode add a reference Peak List
   - Reference peaks loaded
   - Used as starting positions for detection

3. **Detect Peaks**
   - Click **Detect**
   - Review detected peaks in table
   - **Manually edit/add/remove peaks if needed**

5. **Fit First Spectrum**
   - Click **Fit Spectrum**
   - Review quality in Peak Navigator
   - **Verify all peaks fit well**

   ## Process Series

   ### 4. Automated Processing

   **Click "Process Series" button**

   **What Happens:**

   ```
   For each spectrum (2, 3, ..., N):
       1. Load spectrum data
       2. Transfer reference peak positions
       3. Transfer reference parameters (fix_positions, etc.)
       4. Run clustering (identical thresholds → identical clusters)
       5. Fit all peaks with PS2D (positions fixed)
       6. Store results
       7. Update progress bar

   After all spectra processed:
       8. Consolidate results table
       9. Enable export options
   ```

   **Progress Indicators:**
   - Progress bar shows current spectrum
   - Console output shows detailed fitting steps
   - Estimated time remaining displayed

   **Parallel Processing:**
   - If enabled, clusters within each spectrum processed in parallel
   - Spectra processed sequentially (not parallelized across spectra)
   - Typical speedup: 2-3× vs sequential

   ### 5. Monitor Quality

   **During Processing:**
   - Watch console for warning messages
   - Check for consistent R² values across series
   - Look for outlier spectra (may indicate experimental issues)

   **After Processing:**
   - Review Peak Table (one row per peak per spectrum)
   - Sort by quality to find problematic fits
   - Use color coding:
     - Green: Excellent/Good (R² ≥ 0.8)
     - Yellow: Fair (R² ≥ 0.5)
     - Red: Poor/Failed (R² < 0.5)

   ---

   ## Visualize Results

   ### Overlay Plot

   **Features:**
   - Open Multi Spectrum Overlay viewer
   - All spectra overlaid with different colors
   - Toggle visibility per spectrum (checkboxes in list)
   - Zoom/pan with matplotlib toolbar
   - Peaks marked at reference positions

   **Color Assignment:**
   - Automatic color cycling (distinguishable palette)
   - Manual color change: Right-click spectrum → Choose Color
   - Intensity normalization: Tools → Normalize Intensities

   ### Peak Navigator

   **Browse Peaks Across Series:**
   1. Click peak in Peak Navigator
   2. Overlay plot centers on that peak
   3. See peak intensity variation across spectra
   4. Voigt Analysis tab shows fit quality per spectrum

   ---
---

## Advanced users - Fix Positions vs Fix Linewidths

### Fix Positions (ALWAYS ON for Series)

**Purpose**: Ensure peaks stay at same chemical shifts across series

**Effect:**
- Stage 2 (position refinement) skipped entirely
- Positions locked to reference values (zero drift)
- Only intensities and linewidths allowed to vary

**When to Use:**
- **Always** for series integration
- Assumes chemical environment is constant
- Assumes no pH/temperature-induced shifts

**When NOT to Use:**
- If chemical shifts expected to change (pH titration with ionizable groups)
- If temperature effects on chemical shift are significant
- Single-spectrum analysis (no need to constrain)

### Fix Linewidths (OPTIONAL for Series)

**Purpose**: Force consistent peak shapes across series

**Effect:**
- Stage 1 (linewidth optimization) skipped
- Linewidths locked to reference values
- Only intensities allowed to vary

**When to Use:**
- Relaxation experiments where linewidths should be constant
- High S/N data where linewidths are well-determined
- Want to isolate intensity changes only

**When NOT to Use:**
- Exchange experiments (linewidths expected to vary)
- Temperature series (viscosity affects linewidths)
- Binding experiments (linewidths report dynamics changes)

**Recommendation**: Leave OFF unless you have specific reason to constrain linewidths.

---


## Advanced Features

### 7. Manual Peak Editing

**Before Processing Series:**
1. Configure first spectrum as usual
2. Review detected peaks in table
3. **Edit peak table**:
   - Add missing peaks: Right-click → Add Peak
   - Remove artifacts: Select row → Delete
   - Adjust positions: Edit F1/F2 cells directly
4. **Set as Reference** with edited peaks
5. Process series with curated peak list

**Benefit:** Ensures consistent peak tracking even if detection misses some peaks in later spectra.

### 8. Per-Spectrum Parameter Override

**If one spectrum needs different parameters:**
1. Process series as usual
2. Note problematic spectrum
3. Right-click spectrum in list → **Refit with Custom Parameters**
4. Adjust parameters in dialog
5. Refit that spectrum only
6. Results updated in table

**Use Case:** One spectrum has lower S/N → increase window scale just for that spectrum.

### 9. Reference Spectrum Selection

**Change reference after initial processing:**
1. Select different spectrum in list
2. Click **Set as Reference**
3. Option: **Re-process All with New Reference**
4. Peaks transferred from new reference
5. Series re-processed

**Use Case:** First spectrum has poor quality, better to use spectrum #5 as reference.

---

## Troubleshooting

### Problem: Peaks Shift Across Series (Even with Fix Positions)

**Symptom:** Peak positions vary by > 0.01 ppm across spectra

**Causes:**
1. Fix Positions checkbox not checked before processing
2. Different overlap thresholds used (changed between spectra)
3. Bug in parameter transfer

**Solution:**
- **Verify** Fix Positions checkbox is ON in first spectrum
- **Re-process series** from scratch
- Check console for "Stage 2 skipped" message (confirms Fix Positions active)

### Problem: Some Peaks Missing in Later Spectra

**Symptom:** Peak present in spectrum 1, absent in spectrum 5

**Causes:**
1. Peak intensity below noise in spectrum 5
2. Peak genuinely absent (exchange, binding event)
3. Detection failed to find peak

**Solution:**
- **Manual peak addition**: Edit first spectrum, add peak manually, re-process
- **Lower noise threshold** for that spectrum (per-spectrum override)
- **Check experimental conditions**: May be real effect (peak broadened beyond detection)

### Problem: Inconsistent R² Across Series

**Symptom:** R² = 0.92 in spectrum 1, R² = 0.65 in spectrum 3 for same peak

**Causes:**
1. Lower S/N in spectrum 3
2. Peak shape changed (linewidth variation)
3. Overlap changed (neighboring peak appeared/disappeared)

**Solution:**
- **Review Voigt Analysis** for spectrum 3 (click peak in navigator)
- If lower S/N: Accept lower R² (quality still "Fair" if > 0.5)
- If peak shape changed: Disable Fix Linewidths
- If overlap: Check for missing peaks in detection

### Problem: Parallel Processing Slower Than Sequential

**Symptom:** Processing time longer with parallel ON

**Causes:**
1. Dataset too small (< 50 peaks)
2. Overhead of worker spawning exceeds benefit
3. Disk I/O bottleneck (all workers reading spectrum file)

**Solution:**
- Disable parallel for small datasets
- Use parallel only for > 100 peaks per spectrum
- Pre-load all spectra into memory (if RAM allows)

---

## Performance Characteristics

| Dataset | Spectra | Peaks/Spectrum | Time (Sequential) | Time (Parallel) |
|---------|---------|----------------|-------------------|-----------------|
| Small | 10 | 50 | ~5 min | ~3 min |
| Medium | 20 | 125 | ~40 min | ~15 min |
| Large | 50 | 300 | ~6 hours | ~2.5 hours |

**Assumptions:**
- 6-core CPU
- Numba installed (3-5× speedup)
- Fix Positions ON (faster than allowing drift)

**Recommendations:**
- Use parallel for > 100 peaks per spectrum
- Install Numba: `pip install numba`
- Process overnight for large datasets (>30 spectra with >200 peaks)

---

## Example Workflows

### Titration Experiment (Ligand Binding)

**Setup:**
- 15 spectra: 0 mM, 0.1 mM, 0.2 mM, ..., 1.4 mM ligand
- 120 peaks expected
- Want to track intensity changes vs concentration

**Workflow:**
1. Load all 15 spectra
2. Configure first spectrum (0 mM):
   - Detect peaks (reference list from apo state)
   - Fix Positions: ON
   - Fix Linewidths: OFF (linewidths may change with binding)
3. Process Series
4. Export CSV
5. In Excel/R:
   - Pivot table: rows = ligand concentration, columns = peak assignment
   - Plot intensity vs [ligand] for each peak
   - Fit binding curves (Kd estimation)

### Relaxation Experiment (T1)

**Setup:**
- 12 spectra: delay times 10 ms, 50 ms, 100 ms, ..., 1000 ms
- 80 peaks expected
- Want to extract T1 relaxation times

**Workflow:**
1. Load all 12 spectra
2. Configure first spectrum (10 ms):
   - Detect peaks
   - Fix Positions: ON
   - Fix Linewidths: ON (linewidths should be constant)
3. Process Series
4. Export CSV
5. In Python/MATLAB:
   - For each peak: intensity vs delay time
   - Fit exponential: I(t) = I0 × (1 - exp(-t/T1))
   - Extract T1 values

### Temperature Series

**Setup:**
- 8 spectra: 278 K, 283 K, 288 K, ..., 313 K
- 150 peaks expected
- Want to measure temperature coefficients

**Workflow:**
1. Load all 8 spectra
2. Configure first spectrum (278 K):
   - Detect peaks
   - Fix Positions: **OFF** (chemical shifts temperature-dependent)
   - Fix Linewidths: OFF (viscosity changes linewidths)
3. Process Series
4. Export CSV
5. Analysis:
   - Plot pos_f1 vs temperature for each peak
   - Calculate Δδ/ΔT (temperature coefficient in ppb/K)
   - Identify amide protons with high coefficients (solvent-exposed)

---

## Integration with Batch Processing

**Difference from Batch Processing:**

| Feature | Multi-Spectrum Viewer | Batch Processing |
|---------|----------------------|------------------|
| **Interaction** | GUI-based, manual configuration | CLI-based, automated |
| **Visualization** | Overlay plot, peak navigator | Minimal (console output) |
| **Peak Consistency** | Enforced (fix positions) | Optional |
| **Use Case** | Explore series interactively | Process large datasets overnight |
| **Export** | Series CSV (one file) | Per-spectrum CSV (many files) |

**When to Use Which:**
- **Multi-Spectrum Viewer**: Titrations, relaxation, small-medium series (<50 spectra)
- **Batch Processing**: High-throughput screening, large series (>50 spectra), automated pipelines

**See:** [ARCHITECTURE.md](ARCHITECTURE.md) section 6 for batch processing details.

---

## Code References

| Component | File | Key Functions |
|-----------|------|---------------|
| Multi-Spectrum Viewer GUI | `gui/multi_spectrum_viewer.py` | `MultiSpectrumViewer` class |
| Series Processing Logic | `processors/multi_spectrum_processor.py` | `process_series()` |
| Reference Parameter Transfer | `core/core_integrator.py` | `transfer_reference_params()` |
| Overlay Plotting | `gui/visualization.py` | `plot_overlay_spectra()` |

---

## Summary Checklist

**Before Processing Series:**
- [ ] All spectra loaded and visible in list
- [ ] First spectrum configured (detection + fitting)
- [ ] Peaks reviewed and edited if needed
- [ ] Fix Positions: **ON**
- [ ] Fix Linewidths: ON/OFF (based on experiment)
- [ ] Parallel Processing: ON (if > 100 peaks)
- [ ] **Set as Reference** clicked

**During Processing:**
- [ ] Monitor progress bar
- [ ] Watch console for errors
- [ ] Note any warning messages

**After Processing:**
- [ ] Review Peak Table for quality
- [ ] Check for outlier spectra (low R² across many peaks)
- [ ] Verify peak positions consistent (pos_f1, pos_f2 constant across rows)
- [ ] Export results as CSV or JSON

**Analysis:**
- [ ] Import CSV into analysis software
- [ ] Plot intensity vs condition for key peaks
- [ ] Fit binding/relaxation curves
- [ ] Report findings

---

**Next Steps:**
- Practice with sample data: `data_example/series_integration/`
- Explore [GUI_GUIDE.md](GUI_GUIDE.md) for general GUI features
- See [FITTING_LOGIC.md](FITTING_LOGIC.md) for technical details of series processing

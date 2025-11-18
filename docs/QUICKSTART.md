# 5-Minute Quickstart

**TL;DR**: Load spectrum → Detect Peaks → Fit All Peaks → Export. PS2D automatically handles overlaps. Default parameters work for most 15N-HSQC data.

---

## 1. Launch (30 seconds)

```bash
cd lunaNMR_v0o9
python3 launch_lunaNMR.py
```

Click **LunaNMR** in launcher dialog.

---

## 2. Load Data (1 minute)

**Option A: Spectrum + Peak List**
1. **File → Open Spectrum** → Select `.2ii` (Bruker), `.ft` (NMRPipe), or `.ucsf` (SPARKY)
2. **File → Open Peak List** → Select `.list` or `.tab` file
3. Spectrum displays with reference peaks as blue ×

**Option B: Spectrum Only**
1. **File → Open Spectrum**
2. Skip peak list → will detect all peaks

**Sample data**: Use `test_data/` folder if available, or your own 15N-HSQC.

---

## 3. Detect Peaks (30 seconds)

**Default settings work for most data.**

1. Check parameters (left panel):
   - `1H/15N (ppm)`: 0.03 / 0.2 ✓
   - `Noise Threshold`: 0.5 ✓
   - `Centroid Window X/Y`: 0.01 / 0.10 ✓

2. Click **Detect Peaks**

3. Console shows:
   ```
   Found 152 peaks total
   116/125 matched, 9 references retained
   ✅ Detection complete
   ```

4. Peak markers appear on spectrum (red dots)

---

## 4. Fit Peaks (2 minutes)

1. Click **Fit All Peaks**

2. Console shows progress:
   ```
   Cluster 1: 1 peak → 1D fitting
   Cluster 5: 3 peaks → PS2D 2D fitting
   ...
   ✅ 116 peaks fitted (97.9% success)
   ```

3. **Peak Navigator** (right panel) shows quality colors:
   - Green = Excellent/Good (R² ≥ 0.8)
   - Yellow = Fair (R² ≥ 0.5)
   - Red = Poor/Failed (R² < 0.5)

---

## 5. Inspect Results (1 minute)

**Option A: Navigate Peaks**
1. Click any peak in **Peak Navigator**
2. Spectrum auto-centers on peak
3. **Voigt Analysis** tab opens showing:
   - Experimental contours (top left)
   - Fitted peaks (top right)
   - Residual heatmap (bottom left)
   - Quality metrics (R², intensity, linewidths)

**Option B: Browse Table**
1. Main window shows peak table with columns:
   - Assignment
   - Position F1/F2
   - Height
   - Quality
2. Sort by clicking column headers
3. Right-click peak → **Show Voigt Analysis**

---

## 6. Export (30 seconds)

**File → Export Results** → Choose format:

| Format | Use Case |
|--------|----------|
| **CSV** | Peak table for Excel/R/Python |
| **JSON** | Full results (fitted surfaces, residuals) |
| **NMRPipe .tab** | Import into NMRDraw/NMRView |
| **SPARKY .list** | Import into SPARKY |

---

## Common Adjustments

### Too Many Peaks Detected
- Increase `Noise Threshold` (0.5 → 1.0)
- Increase `X×Y pixels` (3×1 → 5×2)

### Peaks at Wrong Positions
- Decrease `Centroid Window X/Y` (0.01/0.10 → 0.005/0.05)
- Check spectrum phasing in original software

### PS2D Fitting Fails
- Increase `Window Scale` slider (Simplified Mode)
- Check overlap ellipses (Tools → Show Ellipses)
- Manually adjust individual peak windows

### Need Consistent Positions Across Spectra
- Check **Fix Positions** before fitting
- Use **Analysis → Start Series Integration** for batches

---

## Keyboard Shortcuts

| Key | Action |
|-----|--------|
| `Ctrl+O` | Open Spectrum |
| `Ctrl+L` | Open Peak List |
| `Ctrl+D` | Detect Peaks |
| `Ctrl+F` | Fit All Peaks |
| `Ctrl+E` | Export Results |
| `←/→` | Navigate Peaks |

---

## Next Steps

- **Batch processing**: Tools → Batch Processing (for multiple spectra)
- **Series integration**: Analysis → Start Series Integration (relaxation/titration)
- **Parameter tuning**: See `docs/GUI_GUIDE.md` for advanced controls
- **Algorithm details**: See `docs/ALGORITHMS.md` for math/equations

---

## Troubleshooting

**"ModuleNotFoundError: No module named 'nmrglue'"**
→ `pip install -r requirements.txt`

**"Numba not installed" warning**
→ `pip install numba` (optional, 3-5× faster fitting)

**Matplotlib backend error**
→ Already handled in `launch_lunaNMR.py` (forces TkAgg)

**PS2D disabled warning**
→ Enable PS2D toggle in left panel (should be ON by default)

---

## Typical Processing Time

| Dataset | Peaks | Time (Sequential) | Time (Parallel) |
|---------|-------|-------------------|-----------------|
| Small | 50 | 30 sec | 20 sec |
| Medium | 125 | 2 min | 45 sec |
| Large | 300 | 8 min | 3 min |

(Tested on 6-core CPU with 15N-HSQC data)

---

**Done!** You've completed a full analysis workflow in <5 minutes.

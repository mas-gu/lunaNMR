# 5-Minute Quickstart

**TL;DR**: Load spectrum → Detect Peaks → Fit All Peaks → Export. PS2D automatically handles overlaps. Default parameters adapt automatically to 15N and 13C data.

---

## 1. Launch (30 seconds)

```bash
cd lunaNMR_v1o0
python launch_lunaNMR.py
```

Click **LunaNMR** in launcher dialog.

---

## 2. Load Data (1 minute)

**Option A: Spectrum + Peak List**
1. **Load data**
2. → Select `.2ii` (Bruker), `.ft` (NMRPipe), or `.ucsf` (SPARKY)
3. → Select `.txt` or `.tab` file
3. Spectrum displays with reference peaks as blue ×

**Option B: Spectrum Only**
1. **Load data**
2. → Select `.2ii` (Bruker), `.ft*` (NMRPipe), or `.ucsf` (SPARKY)
3. Skip peak list → will detect all peaks

---

## 3. Detect Peaks (30 seconds)

**Default settings work for most data.**

1. Click **Detect**

2. Console shows:
   ```
   x
   ```

3. Peak markers appear on spectrum (red dots)

---

## 4. Fit Peaks (2 minutes)

1. Click **Fit Spectrum**

2. Console shows progress:

3. **Peak Navigator** (right panel)

---

## 5. Inspect Results

** Navigate Peaks**
1. Click any peak in **Peak Navigator**
2. Spectrum auto-centers on peak
3. **3D Voigt Analysis** tab opens showing: (click on analysis symbol)
   - Experimental contours (top left)
   - Fitted peaks (top right)
   - Residual heatmap (bottom left)
   - Quality metrics (R², intensity, linewidths)

---


## Button Expert Mode

### All parameters for expect peak detection and 2D voigt fitting. Not needed for regular user.

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

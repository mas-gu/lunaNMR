# LunaNMR Quick Reference Card

**TL;DR**: One-page cheat sheet for common operations and parameters.

---

## Launch
```bash
cd lunaNMR_v1o0
python3 launch_lunaNMR.py
```

---

## Main Workflow

```
Load Spectrum → Detect Peaks → Fit All Peaks → Export Results

---

## Detection Parameters (Defaults)

| Parameter | 15N-HSQC | 13C-HSQC | Purpose |
|-----------|----------|----------|---------|
| 1H/15N (ppm) | 0.03 / 0.2 | 0.03 / 0.1 | Search window |
| Noise Threshold | 0.5 | 0.5 | S/N multiplier |
| Centroid X/Y (ppm) | 0.01 / 0.10 | 0.01 / 0.05 | Max position shift |

---

## Fitting Controls

| Toggle | Default | Effect |
|--------|---------|--------|
| PS2D multi-peak | ✓ ON | Auto 2D fitting for overlaps |
| Fix Positions | ☐ OFF | Zero drift (REQUIRED for series) |
| Fix Linewidths | ☐ OFF | Lock widths to reference |
| Parallel Processing | ☐ OFF | 2.7× faster (>100 peaks) |

---

## Overlap Thresholds (Auto-Clustering)

| Nucleus | F1 (ppm) | F2 (ppm) | Meaning |
|---------|----------|----------|---------|
| 15N | ±0.4 | ±0.04 | Peaks within this → cluster |
| 13C | ±0.1 | ±0.04 | Two-circle touching test |

---

## Quality Indicators

| Color | R² Range | Meaning |
|-------|----------|---------|
| 🟢 Green | ≥ 0.8 | Excellent/Good - Acceptable |
| 🟡 Yellow | 0.5-0.8 | Fair - Review recommended |
| 🔴 Red | < 0.5 | Poor/Failed - Needs attention |

---

## PS2D Fitting Stages

```
Stage 1: Linewidths + intensity      → Determine peak shapes (skip if fix_linewidths)
Stage 2: Positions + linewidths      → Refine centers (skip if fix_positions)
Stage 3: Global refinement           → Final polish
```

**Parameter Bounds:**
- Positions: ±1.5× FWHM (F1 max 0.1 ppm, F2 max 0.01 ppm)
- Linewidths: 0.5-5.0× median (F1), 0.5-2.0× median (F2)

---

## Series Integration (Titration/Relaxation)

**Quick Workflow:**
1. Fit Series
2. ✓ Fix Positions ON
3. Set as Reference
4. Process Series
5. Export Series Results → CSV

**Result:** One CSV row per peak per spectrum + table with Volume for all spectra together

---

## File Formats

| Format | Extension | Software |
|--------|-----------|----------|
| Bruker | `.2ii`, `.2rr` | TopSpin |
| NMRPipe | `.ft`, `.pipe` | NMRPipe |
| Varian | `.fid`, `.ft` | VnmrJ |
| SPARKY | `.ucsf` | SPARKY |

---
---

## Performance (150 peaks, 6-core CPU)

| Mode | Detection | Fitting | Total |
|------|-----------|---------|-------|
| Sequential | ~10 sec | ~2 min | ~2.5 min |
| Parallel | ~10 sec | ~45 sec | ~1 min |
| Parallel + Numba | ~10 sec | ~20 sec | ~30 sec |

**Recommendation:** Install Numba for 3-5× speedup
```bash
pip install numba
```

---

## Common Nucleus Configurations

### 15N-HSQC (600 MHz)
```
Detection window: 0.03/0.2 ppm
Overlap threshold: 0.04/0.4 ppm
Centroid window: 0.01/0.10 ppm
Typical linewidths: 0.02/0.10 ppm (1H/15N)
```

### 13C-HSQC (600 MHz)
```
Detection window: 0.03/0.1 ppm
Overlap threshold: 0.04/0.1 ppm
Centroid window: 0.01/0.05 ppm
Typical linewidths: 0.02/0.08 ppm (1H/13C)
```

---

## Advanced Parameters (Developer Mode)

**PS2D Radii** (elliptical data selection):
- 15N: radF1=0.4 ppm, radF2=0.04 ppm
- 13C: radF1=0.1 ppm, radF2=0.04 ppm

**LM Optimizer:**
- Initial λ = 0.001
- Max iterations = 250 per stage
- Convergence tolerance = 1e-6

**Cluster Limits:**
- Max cluster size: 15 peaks
- Max diameter: 2.0 ppm (F1) × 0.2 ppm (F2)

---

## Code References (Quick Lookup)

| Component | File |
|-----------|------|
| Peak Detection | `core/enhanced_peak_picker.py` |
| Clustering | `core/ps2d_exact_overlap_detector.py` |
| PS2D 2D Fitting | `core/ps2d_2d_fitter.py` |
| 1D Fallback | `core/enhanced_voigt_fitter.py` |
| Main Orchestrator | `core/core_integrator.py` |
| Multi-Spectrum | `gui/multi_spectrum_viewer.py` |
| Parallel Processing | `core/parallel_voigt_processor.py` |

---

## Documentation Index

| Topic | File |
|-------|------|
| Installation | [INSTALLATION.md](INSTALLATION.md) |
| 5-min tutorial | [QUICKSTART.md](QUICKSTART.md) |
| GUI features | [GUI_GUIDE.md](GUI_GUIDE.md) |
| Series integration | [MULTI_SPECTRUM_GUIDE.md](MULTI_SPECTRUM_GUIDE.md) |
| Algorithms | [ALGORITHMS.md](ALGORITHMS.md) |
| Fitting logic | [FITTING_LOGIC.md](FITTING_LOGIC.md) |
| Architecture | [ARCHITECTURE.md](ARCHITECTURE.md) |
| PS2D tuning | [advanced/PS2D_TUNING.md](advanced/PS2D_TUNING.md) |

---

## Support

**Documentation:** `/lunaNMR_v1o0/docs/`

**Issues:** Check console output for error messages

**Performance:** Install Numba (`pip install numba`) for 3-5× speedup

---

**Print this page** for desk reference while using LunaNMR.

# LunaNMR v0.9 Documentation

**TL;DR**: Start with [QUICKSTART.md](QUICKSTART.md) → Explore [GUI_GUIDE.md](GUI_GUIDE.md) for features → Consult [ALGORITHMS.md](ALGORITHMS.md) for technical details.

---

## Documentation Structure

### User Guides

| Document | Purpose | Audience | Time |
|----------|---------|----------|------|
| [INSTALLATION.md](INSTALLATION.md) | Setup and dependencies | New users | 5 min |
| [QUICKSTART.md](QUICKSTART.md) | 5-minute tutorial | New users | 5 min |
| [GUI_GUIDE.md](GUI_GUIDE.md) | Complete GUI reference | All users | 15 min |
| [MULTI_SPECTRUM_GUIDE.md](MULTI_SPECTRUM_GUIDE.md) | Series integration workflow | Advanced users | 10 min |

### Developer Guides

| Document | Purpose | Audience | Time |
|----------|---------|----------|------|
| [ARCHITECTURE.md](ARCHITECTURE.md) | System design and data flow | Developers | 20 min |
| [ALGORITHMS.md](ALGORITHMS.md) | Peak detection and PS2D fitting | Developers/Scientists | 15 min |
| [FITTING_LOGIC.md](FITTING_LOGIC.md) | Two-phase fitting workflow | Developers | 15 min |

### Advanced Topics

| Document | Purpose | Audience | Time |
|----------|---------|----------|------|
| [advanced/PS2D_TUNING.md](advanced/PS2D_TUNING.md) | Fine-tuning PS2D parameters | Experts | 30 min |
| [advanced/RADII_ANALYSIS.md](advanced/RADII_ANALYSIS.md) | Adaptive radii research | Researchers | 30 min |

---

## Quick Navigation

### I want to...

**...get started quickly**
→ [QUICKSTART.md](QUICKSTART.md)

**...understand the main GUI**
→ [GUI_GUIDE.md](GUI_GUIDE.md) sections 1-7

**...process multiple spectra (series integration)**
→ [MULTI_SPECTRUM_GUIDE.md](MULTI_SPECTRUM_GUIDE.md)

**...understand how peak detection works**
→ [ALGORITHMS.md](ALGORITHMS.md) section 1

**...understand how PS2D 2D fitting works**
→ [ALGORITHMS.md](ALGORITHMS.md) section 2
→ [FITTING_LOGIC.md](FITTING_LOGIC.md)

**...understand clustering and overlap detection**
→ [FITTING_LOGIC.md](FITTING_LOGIC.md) section 2

**...modify the system architecture**
→ [ARCHITECTURE.md](ARCHITECTURE.md)

**...tune PS2D fitting parameters**
→ [advanced/PS2D_TUNING.md](advanced/PS2D_TUNING.md)

---

## Key Concepts

### PS2D 2D Simultaneous Fitting

LunaNMR uses **PS2D** (Port from PINT's Simultaneous 2D fitting) as the core algorithm for resolving overlapping peaks. This is NOT an optional feature - it's the primary fitting engine.

**When it's used:**
- Automatically for all peaks separated by < 0.08 ppm (1H) or < 0.8 ppm (15N)
- Triggered by hierarchical overlap clustering
- Processes entire clusters (not individual peaks)

**How it works:**
- 5-stage Levenberg-Marquardt optimization
- Fits all peaks in a cluster simultaneously
- Uses elliptical data selection windows
- Respects GUI constraints (Fix Positions, Fix Linewidths)

See [FITTING_LOGIC.md](FITTING_LOGIC.md) for detailed workflow.

### Two-Phase Fitting Logic

**Phase 1: Peak Detection and Clustering**
- Prominence-based peak detection
- Top-contour centroid refinement
- Hierarchical graph-based clustering
- Creates disjoint clusters (each peak in exactly ONE cluster)

**Phase 2: Fitting Execution**
- Single-peak clusters → 1D cross-section fitting (fallback)
- Multi-peak clusters → PS2D 2D simultaneous fitting
- Parallel processing distributes clusters across cores
- Results consolidated by peak_number (integer matching)

See [FITTING_LOGIC.md](FITTING_LOGIC.md) for flowcharts.

### Quality Assessment

| Quality | R² Range | Color | Action |
|---------|----------|-------|--------|
| Excellent | ≥ 0.9 | Green | High confidence |
| Good | [0.8, 0.9) | Green | Acceptable |
| Fair | [0.5, 0.8) | Yellow | Review recommended |
| Poor | [0.2, 0.5) | Red | Manual inspection needed |
| Failed | < 0.2 | Red | Refit with adjusted parameters |

---

## File Format Support

| Format | Extension | Software | Status |
|--------|-----------|----------|--------|
| Bruker TopSpin | `.2ii`, `.2rr` | TopSpin | ✅ Full support |
| NMRPipe | `.ft`, `.pipe` | NMRPipe | ✅ Full support |
| Varian/Agilent | `.fid`, `.ft` | VnmrJ | ✅ Full support |
| SPARKY | `.ucsf` | SPARKY | ✅ Full support |

All formats loaded via `nmrglue` with automatic axis ordering detection.

---

## Common Workflows

### Single Spectrum Analysis
1. File → Open Spectrum
2. File → Open Peak List (optional)
3. Detect Peaks
4. Fit All Peaks
5. Export Results

### Series Integration (Titration/Relaxation)
1. Analysis → Start Series Integration
2. Select folder with multiple spectra
3. Configure parameters on first spectrum
4. Process Series
5. Export Series Results

See [MULTI_SPECTRUM_GUIDE.md](MULTI_SPECTRUM_GUIDE.md) for detailed workflow.

### Batch Processing
1. Tools → Batch Processing
2. Select folder
3. Configure detection/fitting parameters
4. Run batch
5. Review results table

---

## Troubleshooting

**Problem**: PS2D fitting fails for some peaks
→ Check overlap ellipses: Tools → Show Ellipses
→ Increase window scale slider
→ See [advanced/PS2D_TUNING.md](advanced/PS2D_TUNING.md) section 4

**Problem**: Peaks detected at wrong positions
→ Tighten centroid windows (decrease X/Y ppm values)
→ See [GUI_GUIDE.md](GUI_GUIDE.md) section 3.2

**Problem**: Slow performance with many peaks
→ Enable parallel processing (Use Parallel Processing checkbox)
→ Install Numba: `pip install numba` (3-5× speedup)
→ See [INSTALLATION.md](INSTALLATION.md) section 4

**Problem**: Inconsistent results across series
→ Enable "Fix Positions" before fitting
→ Use reference peak list from first spectrum
→ See [MULTI_SPECTRUM_GUIDE.md](MULTI_SPECTRUM_GUIDE.md) section 3

---

## Performance Characteristics

| Dataset Size | Detection | Fitting (Sequential) | Fitting (Parallel) | With Numba |
|--------------|-----------|---------------------|-------------------|------------|
| 50 peaks | ~5 sec | ~30 sec | ~20 sec | ~10 sec |
| 125 peaks | ~10 sec | ~2 min | ~45 sec | ~20 sec |
| 300 peaks | ~20 sec | ~8 min | ~3 min | ~90 sec |

**Recommendations:**
- Enable parallel processing for > 100 peaks
- Install Numba for production use (especially for large overlap clusters)
- Use sequential mode for < 50 peaks (parallel overhead not worth it)

---

## Contributing

When updating documentation:
- Keep it SHORT and FOCUSED
- Use TL;DR at top of each file
- Add concrete examples, not abstract descriptions
- Include file:line references for code
- Update this README when adding new docs

---

## Version History

- **v0.9** (2025): PS2D 2D simultaneous fitting, parallel processing, series integration
- **v0.8**: Enhanced peak detection with centroid refinement
- **v0.7**: Initial GUI implementation

---

**Next Steps**: Read [QUICKSTART.md](QUICKSTART.md) or dive into [GUI_GUIDE.md](GUI_GUIDE.md).

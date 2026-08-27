# LunaNMR v1.0 Documentation

New here? Start with [Getting Started](#getting-started) → then [ALGORITHMS.md](ALGORITHMS.md) for the technical details.

---

## Getting Started

Load spectrum → **Detect** → **Fit Spectrum** → inspect → export. PS2D handles overlaps automatically; defaults adapt to ¹⁵N and ¹³C.

1. **Launch** — `cd lunaNMR_v1o0 && python launch_lunaNMR.py`; the main window opens directly.
2. **Load data** — click **Load Data** and pick a spectrum (`.2ii`/`.2rr` Bruker, `.ft`/`.pipe` NMRPipe, `.ucsf` SPARKY, `.fid`/`.ft` Varian). The same dialog optionally loads a peak list; reference peaks show as blue ×. Without a list, all peaks are detected.
3. **Detect** — click **Detect**. Markers appear as red dots.
4. **Fit** — click **Fit Spectrum**. Progress logs to the console; fitted peaks populate the **Peak Navigator** (right panel).
5. **Inspect** — click a peak in the Navigator to auto-center. The **Voigt Analysis** tab shows contours, fitted peaks, the residual heatmap, and quality metrics (R², intensity, linewidths).

All detection/fitting parameters live behind the **Expert Mode** button (not needed for routine use). See [QUICK_REFERENCE.md](QUICK_REFERENCE.md) for the cheat-sheet.

---

## Documentation

| Document | Purpose |
|----------|---------|
| [INSTALLATION.md](INSTALLATION.md) | Setup and dependencies |
| [QUICK_REFERENCE.md](QUICK_REFERENCE.md) | Parameter cheat-sheet |
| [lunaNMR_guide.html](lunaNMR_guide.html) | Complete GUI/feature reference |
| [MULTI_SPECTRUM_GUIDE.md](MULTI_SPECTRUM_GUIDE.md) | Series integration workflow |
| [CLI.md](CLI.md) | Headless command-line reference |
| [CLI_AGENT.md](CLI_AGENT.md) | Machine contract for driving the CLI from an agent: JSON shapes, exit codes, silent-corruption gotchas |
| [CLI_AGENTS_DEEP/](CLI_AGENTS_DEEP/) | Long-form agent runbooks — relaxation and affinity, with the physical QC bands |
| [ARCHITECTURE.md](ARCHITECTURE.md) | System design and data flow |
| [ALGORITHMS.md](ALGORITHMS.md) | Peak detection, clustering, PS2D fitting |

---

## Key Concepts

### Analysis Modules (Modules menu)
- **DynamiXs Relaxation Analysis** — T₁/T₂, methyl T₂ (bi-exponential), spectral density, model-free, CPMG.
- **Kd / Titration Analysis** — CSP (1:1 quadratic isotherm) and intensity exponential-decay binding fits.

Both are scriptable headless (`python -m lunaNMR dynamixs …` / `kd`) — see [CLI.md](CLI.md).

### Project Bundles (`.lunaNMR`)
State-based save/load. Fit surfaces are stored slim and reconstructed on load. DynamiXs and Kd analyses persist as multiple named runs, reopenable from the Project Browser (File ▸ Project Contents…).

### PS2D 2D Simultaneous Fitting
The primary fitting engine for overlapping peaks — not optional.
- Triggered automatically by hierarchical overlap clustering for peaks separated by < 0.08 ppm (¹H) or < 0.8 ppm (¹⁵N).
- 5-stage Levenberg-Marquardt; fits every peak in a cluster simultaneously over elliptical data windows.
- Respects the Expert-Mode **Fix Positions** / **Fix Linewidths** flags.

### Two-Phase Fitting
- **Phase 1 (detect + cluster):** prominence detection → top-contour centroid refinement → hierarchical graph-based clustering into disjoint clusters (each peak in exactly one).
- **Phase 2 (fit):** single-peak clusters → 1D cross-section fitting; multi-peak clusters → PS2D. Parallel mode distributes whole clusters across cores; results are identical to sequential.

See [ALGORITHMS.md](ALGORITHMS.md) for detection, clustering, and fitting details.

### Quality Color Coding
Peak markers are colored by fit R² (thresholds vary by nucleus; GUI marker defaults below):

| Quality | R² | Color |
|---------|-----|-------|
| Excellent | ≥ 0.9 | green |
| Good | ≥ 0.8 | green/blue |
| Fair | ≥ 0.5 | orange |
| Poor / unfitted | < 0.5 | red |

---

## File Formats

Loaded via `nmrglue` with automatic axis-ordering detection.

| Software | Extension |
|----------|-----------|
| Bruker TopSpin | `.2ii`, `.2rr` |
| NMRPipe | `.ft`, `.pipe` |
| Varian/Agilent | `.fid`, `.ft` |
| SPARKY | `.ucsf` |

---

## Common Workflows

**Single spectrum:** Load Data → Detect → Fit Spectrum → File ▸ Export Results.

**Series (titration/relaxation):** click **Fit Series** → configure in the Series Integration dialog → process → results land in a `series_results_*/` folder (read by DynamiXs / Kd / overlay viewer). See [MULTI_SPECTRUM_GUIDE.md](MULTI_SPECTRUM_GUIDE.md).

**Batch (headless):** `python -m lunaNMR batch …` — see [CLI.md](CLI.md).

---

## Troubleshooting

| Problem | Try |
|---------|-----|
| PS2D fitting fails for some peaks | Expert Mode ▸ **Show Ellipses** to inspect windows; see [ALGORITHMS.md](ALGORITHMS.md) |
| Peaks at wrong positions | Tighten centroid windows (Expert Mode search X/Y ppm); see [lunaNMR_guide.html](lunaNMR_guide.html) |
| Slow with many peaks | Expert Mode ▸ **Use Parallel Processing**; ensure Numba is installed (see [INSTALLATION.md](INSTALLATION.md)) |
| Inconsistent series results | Enable **Fix Positions**; use a reference peak list from the first spectrum; see [MULTI_SPECTRUM_GUIDE.md](MULTI_SPECTRUM_GUIDE.md) |

Rough guide: sequential ~120 s / parallel ~45 s for 150 peaks (2.7× with 6 cores). Numba is required for large clusters — a 9-peak cluster is ~45–90 s with it, minutes without.

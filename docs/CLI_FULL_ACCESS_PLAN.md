# Plan: full CLI access for lunaNMR

Status: **implemented (2026-07-06).** Drafted 2026-06-17.

## Implemented — unified `python -m lunaNMR` (lunaNMR/cli.py + __main__.py)
All reasonable batch/analysis workflows are now headless subcommands, TDD'd in
`tests/test_cli.py` (24 tests):
- `kd` — Kd from a titration series (CSP quadratic isotherm / intensity decay).
- `series` — multi-spectrum / titration processing (`--mode time|titration`,
  `--peak-source reference|cascade|detected|independent`).
- `dynamixs t1t2` / `dynamixs methyl-t2` — relaxation fitting.
- `export kd` — headless CSP/intensity figures + summary.csv from a kd fit JSON.
- `project inventory` / `project remove` — inspect / prune a `.lunaNMR` bundle.
- `batch` — folds in the existing `batch_processing` CLI (detect + Voigt/PS2D fit).

Enabling change: `lunaNMR/__init__.py` made lazy (PEP 562 `__getattr__`) so the CLI
is headless (no eager GUI/core import; no `nmrglue` needed for `kd`/`export`/`project`).
`series`/`batch`/`dynamixs` still require the scientific stack (they read spectra / fit).

**Intentionally NOT exposed:** full headless project *save/load* of a live session —
a CLI has no in-memory session to serialize. `ProjectManager` already works with a
duck-typed session object (verified), so a composite "bundle these CLI outputs" verb
could be added later if a real need appears. The `pyproject.toml` / `lunanmr` console
script (vs `python -m lunaNMR`) is also still deferred.

The original plan text follows.

## Why this is tractable
The heavy analysis engines are already **Qt-free** and route output through
`utils/output_manager.py` (which falls back to `print()` when no GUI callback is set):
- `processors/multi_spectrum_processor.py` — series/titration processing (0 Qt imports)
- `core/core_integrator.py` — detection + Voigt/PS2D fitting (0 Qt imports)
- `dynamiXs_v2o0/dynamiXs_Kd/kd_fit.py` — Kd analysis; `run_kd_analysis_with_params`
  already takes a params dict and writes JSON + txt (0 Qt imports)

So "full CLI" is mostly **wrapping existing functions**, not rewriting algorithms. The one
real coupling is `utils/project_manager.py`, which reads session state off the GUI main
window (`self.main_window.{peaks, fit_results, saved_series, kd_analyses, dynamixs_state}`,
~6 references).

## Current CLI state (2026-06-17)
- **Accessible**: batch detect + Voigt/PS2D fitting over a folder
  (`python -m lunaNMR.batch_processing`); the whole ML pipeline (`scripts/`); some
  dynamiXs plotting scripts.
- **GUI-only**: multi-spectrum series/titration workflow (`process_nmr_series`, called only
  by `gui/dialogs/series_integration_dialog.py`); Kd analysis; dynamiXs T1/T2 & methyl-T2
  fitting; `.lunaNMR` project save/load; interactive viewers/editing/export.

## Scope
"Full CLI" = every batch/analysis workflow. NOT the inherently interactive features
(manual peak picking, real-time adjustment, click-to-exclude viewers) — those stay GUI.

## Architecture
- Single entry point `python -m lunaNMR` with **argparse subcommands**, each wrapping the
  Qt-free core. Console script `lunanmr = lunaNMR.cli:main`.
- Complex workflows take a **JSON/YAML config** (mirroring what the GUI dialogs collect),
  plus flags for common cases.
- Force matplotlib **Agg** backend for headless plotting; make the few top-level PySide6
  imports lazy so the CLI never needs a `QApplication`.

## Phases (quick wins first)
1. **`lunanmr kd` (LOW)** — wrap `run_kd_analysis_with_params` (already CLI-shaped):
   `--input tidy.csv --conc … --p0 --alpha --observable --intensity-from --out`. Proof of
   concept; validates the whole approach.
2. **`lunanmr series` (MED)** — wrap `MultiSpectrumProcessor.process_nmr_series`:
   `--spectra <folder/glob> --peaks ref.csv --out <dir> --mode time|titration
   --peak-source reference|cascade|detected`. Linchpin (currently only in the GUI dialog);
   produces the `series_results_*` outputs that feed Kd & dynamiXs.
3. **`lunanmr fit` / `detect` (LOW)** — re-expose the existing `batch_processing` CLI under
   the unified entry point.
4. **`lunanmr dynamixs` (MED)** — wrap the T1/T2 and methyl-T2 fitters (already mostly
   standalone in `dynamiXs_T1_T2/`): `dynamixs t1t2` / `dynamixs methyl-t2`.
5. **`lunanmr project` save/load (HIGH — the refactor)** — introduce a plain `SessionState`
   object holding peaks/fit_results/series/kd_analyses/dynamixs_state; have the GUI main
   window expose it and refactor `ProjectManager` to operate on `SessionState` instead of
   the `QMainWindow`. Then CLI can save/load `.lunaNMR` bundles headlessly. Needs careful
   TDD; only structural change in the plan.
6. **Headless export/report (MED)** — factor viewer plotting into pure functions
   (data → figure) so CLI can emit CSP/fit figures (SVG/PDF, editable text) + a summary
   report without the Qt viewers.

## Cross-cutting
- Headless guard: lazy PySide6 imports; Agg backend.
- One config schema describing a full pipeline (spectra → series → kd) for reproducible runs.
- End-to-end CLI tests per subcommand on synthetic data, asserting output files.
- README CLI section + `--help` epilogs.

## Effort / risk
Phases 1–4 are wrappers (days each, low risk — engines are ready). Phase 5 (project
persistence) is the real work and the only thing needing a refactor. Phase 6 is medium.

Recommended start: Phase 1 (`kd`) + the unified `python -m lunaNMR` skeleton.

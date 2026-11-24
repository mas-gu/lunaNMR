# Installation Guide – lunaNMR v1.0

Runs on macOS, Linux, Windows. All commands run from the repository root (the `lunaNMR/` folder `git clone` creates).

lunaNMR is used from its clone: run in place, or installed from it in **editable** mode (`pip install -e .`). A regular `pip install .` is not supported — see [section 5](#5-run-from-the-clone-or-install-it-in-editable-mode).

## 1. Prerequisites

- **Python** ≥ 3.9 (64‑bit).
- **Tk** — imported by the launcher, which uses it only for its startup error dialogs. Bundled with most Python builds; on Linux install `python3-tk` (Debian/Ubuntu) or `python3-tkinter` (Fedora/RHEL).
- **Git** — to clone the repo.

## 2. Get the Sources

```bash
git clone https://github.com/mas-gu/lunaNMR.git
cd lunaNMR
```

## 3. Virtual Environment (recommended)

```bash
python3 -m venv .venv
source .venv/bin/activate     # Linux/macOS
.venv\Scripts\activate        # Windows PowerShell
```

## 4. Install Dependencies

All runtime deps are in `requirements.txt`:

| Package | Notes |
|---|---|
| numpy ≥1.20, pandas ≥1.3, scipy ≥1.7 | core |
| lmfit ≥1.3.4 | DynamiXs T1/T2 and methyl T2 fitting |
| matplotlib ≥3.5 | backend forced to `QtAgg` |
| PySide6 ≥6.6 | Qt6 GUI — **required** |
| numba ≥0.57 | **required** — 3–5× faster 2D multi-peak fitting |
| nmrglue ≥0.9 | NMR file reading |

```bash
pip install -r requirements.txt
```

### Troubleshooting

- **`numba` wheel missing**: `pip install --upgrade pip`, retry. On Apple silicon use an arm64 Python.
- **`tkinter` missing**: install the OS package above.
- **Permission errors**: use a venv or `pip install --user`.

## 5. Run from the Clone, or Install It in Editable Mode

These are the two supported ways to use lunaNMR, and both keep the clone as the code that runs.

**From the clone, with no further step.** With the dependencies installed,
`python3 launch_lunaNMR.py` and `python -m lunaNMR <subcommand>` work from the repository root.

**Editable install**, to run it from any directory and to get the console scripts `lunanmr`
and `lunanmr-mcp` declared in `pyproject.toml`:

```bash
pip install -e .        # from the repository root; the clone stays the installed code
```

Without it `which lunanmr` finds nothing, and every `lunanmr …` line in the docs must be
read as `python -m lunaNMR …` run from the repository root. Keep the clone where it is
afterwards: an editable install points at it.

**A regular install is not supported.** `pip install .` (or a wheel built from the
repository) copies the `lunaNMR` package alone. DynamiXs, Kd / titration analysis and 1D
integration are not inside that package: they live in the `modules/` folder beside it and
are left behind. From such a copy `series`, `diagnose`, `batch` and `project` still run,
but the GUI does not open, and `dynamixs`, `kd`, `export kd`, `integrate-1d` (and
`project build` on a run holding a Kd fit) stop with one message naming the missing folder
and the editable install; under `--format json` it is the error envelope with `error.type`
`AnalysisModulesNotFound`. If you see that message, reinstall with `pip install -e .`
from the clone.

## 6. Verify

```bash
python3 lunaNMR/validation/verify_installation.py
```

Reports the Python version, the required packages, the internal modules, the analysis modules beside the package, whether the main GUI imports, and whether the fitting pipeline initialises. It exits 1 on a regular install, naming the missing `modules/` folder.

## 7. Launch

- **GUI**: `python3 launch_lunaNMR.py` — opens the LunaNMR main window. DynamiXs, Kd / Titration Analysis and 1D Integration are entries of its **Modules** menu.
- **Headless CLI**: `python -m lunaNMR <subcommand>`. See [CLI.md](CLI.md) for the subcommand list and flags.

## 8. Sample Data

Example spectra ship under `data_example/` (a 1D series in `1D/`, and an eight-point 2D T1 series split between `2DHSQC/`, which also holds the peak lists, and `more/`). Copy them into a separate workspace if you prefer not to modify the repo.

## 9. Update

```bash
git pull
pip install --upgrade -r requirements.txt
```

An editable install follows the clone, so there is nothing to reinstall. Re-run the verification script to confirm Numba still imports.

## 10. Cleanup

Deactivate and delete the venv, then remove the repo directory. No system-wide services.

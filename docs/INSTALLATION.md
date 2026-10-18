# Installation Guide – lunaNMR v1.0

Runs on macOS, Linux, Windows. All commands run from the `lunaNMR_v1o0/` directory.

## 1. Prerequisites

- **Python** ≥ 3.9 (64‑bit).
- **Tk** — imported by the launcher, which uses it only for its startup error dialogs. Bundled with most Python builds; on Linux install `python3-tk` (Debian/Ubuntu) or `python3-tkinter` (Fedora/RHEL).
- **Git** — to clone the repo.

## 2. Get the Sources

```bash
git clone https://github.com/mas-gu/lunaNMR.git
cd lunaNMR/lunaNMR_v1o0
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
| matplotlib ≥3.5 | backend forced to `QtAgg` |
| PySide6 ≥6.6 | Qt6 GUI — **required** |
| numba ≥0.57 | **required** — 3–5× faster 2D multi-peak fitting |
| networkx ≥2.5 | listed requirement; no module imports it |
| nmrglue ≥0.9 | NMR file reading |
| psutil ≥5.8 | optional — performance tracking |

```bash
pip install -r requirements.txt
```

### The `lunanmr` command (optional)

`python -m lunaNMR <subcommand>` works from the project directory with no extra step.
The shorter `lunanmr <subcommand>` is a console script declared in `pyproject.toml`, and
it exists only after an editable install:

```bash
pip install -e .        # from lunaNMR_v1o0/, adds `lunanmr` to PATH
```

Without it `which lunanmr` finds nothing, and every `lunanmr …` line in the docs must be
read as `python -m lunaNMR …`.

### Troubleshooting

- **`numba` wheel missing**: `pip install --upgrade pip`, retry. On Apple silicon use an arm64 Python.
- **`tkinter` missing**: install the OS package above.
- **Permission errors**: use a venv or `pip install --user`.

## 5. Verify

```bash
python3 lunaNMR/validation/verify_installation.py
```

Reports the Python version, the required packages (and optional `psutil`), the internal modules, whether the main GUI imports, and whether the fitting pipeline initialises.

## 6. Launch

- **GUI**: `python3 launch_lunaNMR.py` — opens the LunaNMR main window. DynamiXs, Kd / Titration Analysis and 1D Integration are entries of its **Modules** menu.
- **Headless CLI**: `python -m lunaNMR <subcommand>`. See [CLI.md](CLI.md) for the subcommand list and flags.

## 7. Sample Data

Example spectra ship under `data_example/` (a 1D series in `1D/`, and an eight-point 2D T1 series split between `2DHSQC/`, which also holds the peak lists, and `more/`). Copy them into a separate workspace if you prefer not to modify the repo.

## 8. Update

```bash
git pull
pip install --upgrade -r requirements.txt
```

Re-run the verification script to confirm Numba still imports.

## 9. Cleanup

Deactivate and delete the venv, then remove the repo directory. No system-wide services.

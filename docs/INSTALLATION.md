# Installation Guide – lunaNMR v0.9

Follow these steps to get lunaNMR running on macOS, Linux, or Windows with an activated Python environment and the correct native dependencies.

## 1. Prerequisites

- **Python**: 3.8 – 3.11 (64‑bit). Install from python.org or your package manager.
- **Git**: if you plan to clone from version control.
- **TK GUI toolkit**: included in most Python distributions. On some Linux distros install `python3-tk` (Ubuntu/Debian) or `python3-tkinter` (Fedora/RHEL).

## 2. Obtain the Sources

```bash
# Clone the repository
git clone <your fork or upstream>
cd lunaNMR_v0o9
```

If you downloaded an archive, extract it and open the extracted `lunaNMR_v0o9` folder in your terminal.

## 3. (Recommended) Create a Virtual Environment

```bash
python3 -m venv .venv
source .venv/bin/activate        # Linux/macOS
# or
.venv\Scripts\activate          # Windows PowerShell
```

Deactivate later with `deactivate`.

## 4. Install Python Dependencies

The project keeps all runtime dependencies in `requirements.txt`:

```text
numpy>=1.20.0
pandas>=1.3.0
matplotlib>=3.5.0
scipy>=1.7.0
scikit-learn>=1.0.0
numba>=0.57.0          # REQUIRED for PS2D 2D multi-peak fitting
networkx>=2.5
nmrglue>=0.9.0
psutil>=5.8.0          # Optional but used for performance dashboards
```

Install them with pip:

```bash
pip install -r requirements.txt
```

If you are behind a firewall, add `--proxy` as required. To install extras for development or testing, append packages such as `pytest`, `black`, or `flake8` manually.

### Troubleshooting dependency installs

- **`numba` wheel missing**: upgrade pip (`pip install --upgrade pip`) and retry. On Apple silicon, ensure you are using an arm64 Python build.
- **`tkinter` missing**: install OS packages noted above.
- **Permission errors**: use a virtual environment or `pip install --user`.

## 5. Verify the Installation

lunaNMR ships with an automated check:

```bash
python3 lunaNMR/validation/verify_installation.py
```

You should see a report confirming the Python version, required packages, and GUI availability.

## 6. Launching the Application

```bash
python3 launch_lunaNMR.py
```

This opens the launcher window. Choose **LunaNMR** for the main analysis suite (or **DynamiXs** if the optional module is installed).

## 7. Optional Data Packages

- Sample spectra reside under `data_example/` and `ml_training_data/`. Copy or mount your spectra into a separate workspace if you prefer not to modify the repository.
- Batch ML runs write to `batch_ml/` and `test_ml_output/`; ensure you have write permission.

## 8. Updating to a Newer Revision

```bash
git pull
pip install --upgrade -r requirements.txt
```

Re-run the verification script to confirm that compiled extensions (Numba) still import correctly after an upgrade.

## 9. Uninstallation / Cleanup

Deactivate and remove your virtual environment, then delete the `lunaNMR_v0o9` directory. There are no system-wide services to stop.

You are now ready to explore the GUI or the batch processing CLI.

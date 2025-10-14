# Repository Guidelines

## Project Structure & Module Organization
The entry point `launch_lunaNMR.py` spins up the GUI suite, while the core package lives under `lunaNMR/` (with `core/` algorithms, `processors/` workflows, `gui/` widgets, and `utils/` shared services). Batch automation scripts sit in `batch_ml/`, optional analysis modules ship in `modules/`, and reference datasets reside in `data_example/` and `ml_training_data/`. Validation utilities, including install checks, are in `lunaNMR/validation/`; keep new verification scripts there. Documentation straight to `docs/`, and generated reports or notebooks should be staged in `report/`.

## Build, Test, and Development Commands
- `python3 -m pip install -r requirements.txt` — install the runtime stack (GUI, scientific, and plotting deps).
- `python3 lunaNMR/validation/verify_installation.py` — quick smoke test to confirm the toolkit sees its dependencies and fonts.
- `python3 launch_lunaNMR.py` — launch the desktop selector for LunaNMR and DynamiXs modules.
- `python3 -m pytest tests` — run the automated suite; create the `tests/` package if you introduce new coverage.

## Coding Style & Naming Conventions
Follow PEP 8 with 4-space indentation and descriptive, snake_case module and function names. Favor type hints on public APIs and NumPy-style docstrings for processors and integrators. Format patches with `black` and lint with `flake8` (both pinned in `requirements.txt`). GUI resources should keep Tk widget IDs in CamelCase to match existing files; CLI helpers stay in snake_case.

## Testing Guidelines
Use `pytest` for all unit and functional checks. Mirror the runtime tree inside `tests/` (e.g., `tests/gui/test_main_gui.py`) and target the smallest reproducible fixtures—draw inputs from `data_example/` or craft synthetic spectra under `tests/fixtures/`. Name tests after the behavior under scrutiny (`test_detects_overlapping_peaks`). Maintain assertions for numerical tolerances and verify CLI entry points via `pytest -k cli` filters. Run `python3 -m pytest --cov=lunaNMR` before opening a pull request.

## Commit & Pull Request Guidelines
Keep commits concise with imperative subjects (`add simplified fit`). Group logically related changes and avoid bundling generated binaries. Pull requests should describe the scientific motivation, reference linked issues or datasets, note data requirements, and include before/after plots when GUI output changes. Attach logs from `verify_installation.py` or targeted pytest runs for regression-sensitive updates.

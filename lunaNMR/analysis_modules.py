# ABOUTME: Locates the analysis modules (DynamiXs, Kd, 1D integration) kept in `modules/`
# ABOUTME: beside the lunaNMR package, and says how to install lunaNMR when they are absent.

"""The analysis modules are not part of the lunaNMR package: they sit in `modules/`, next
to it in the clone, and import each other by bare name. `pip install .` copies the package
alone, so they are found only from the clone or from an editable install of it.
"""

import os
import sys

MODULES_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'modules'))

# Folder under `modules/` -> what it holds.
ANALYSIS_MODULES = {
    'dynamiXs_v2o0': 'DynamiXs relaxation and Kd / titration analysis',
    'integration_1d_v1o0': '1D integration',
}


class AnalysisModulesNotFound(ImportError):
    """An analysis module's folder is not beside the lunaNMR package."""


def module_folder(*parts):
    """The folder `modules/<parts>`; AnalysisModulesNotFound when it is not there."""
    path = os.path.join(MODULES_ROOT, *parts)
    if not os.path.isdir(path):
        raise AnalysisModulesNotFound(
            f"The lunaNMR analysis modules (DynamiXs, Kd, 1D integration) were not found: "
            f"{path} does not exist. They live in the `modules` folder beside the lunaNMR "
            f"package, which a plain `pip install .` does not copy. Install lunaNMR in "
            f"editable mode from its clone (`pip install -e .`), or run it from the clone.")
    return path


def add_to_import_path(*parts):
    """Put `modules/<parts>` first on sys.path, so its bare-name imports resolve."""
    path = module_folder(*parts)
    if path not in sys.path:
        sys.path.insert(0, path)
    return path

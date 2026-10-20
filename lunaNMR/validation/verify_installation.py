#!/usr/bin/env python3
# ABOUTME: Verifies that lunaNMR's dependencies, package modules and GUI entry point
# ABOUTME: all import, and reports what is missing when they do not.
"""Installation verification for lunaNMR.

Run from anywhere:  python3 lunaNMR/validation/verify_installation.py

Exits 0 when the checkout is usable and 1 otherwise, so it can gate CI as well
as answer "did my install work?".
"""

import importlib
import sys
from pathlib import Path

# Python puts this script's own directory on sys.path, not the caller's working
# directory, so an uninstalled checkout cannot import `lunaNMR` without help.
REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

LAUNCH_SCRIPT = "launch_lunaNMR.py"

# Required at runtime; mirrors the `dependencies` list in pyproject.toml.
REQUIRED_DEPENDENCIES = {
    'numpy': 'Numerical computing',
    'pandas': 'Data manipulation',
    'matplotlib': 'Plotting and visualization',
    'scipy': 'Scientific computing',
    'lmfit': 'Relaxation curve fitting (DynamiXs)',
    'nmrglue': 'NMR data file reading',
    'PySide6': 'Qt GUI framework',
    'numba': '2D multi-peak fitting acceleration',
}

# Absent installs degrade a feature rather than breaking the application.
OPTIONAL_DEPENDENCIES = {}

INTERNAL_MODULES = [
    ('lunaNMR.core.core_integrator', 'Core integration engine'),
    ('lunaNMR.core.enhanced_voigt_fitter', '1D Voigt profile fitting'),
    ('lunaNMR.core.enhanced_peak_picker', 'Peak detection'),
    ('lunaNMR.core.ps2d_2d_fitter', '2D simultaneous multi-peak fitting'),
    ('lunaNMR.core.parallel_voigt_processor', 'Parallel processing'),
    ('lunaNMR.processors.multi_spectrum_processor', 'Series processing'),
    ('lunaNMR.utils.file_manager', 'File handling'),
    ('lunaNMR.utils.config_manager', 'Configuration'),
    ('lunaNMR.utils.project_manager', 'Project save and load'),
    ('lunaNMR.gui.main_window', 'Main window'),
    ('lunaNMR.gui.components.spectrum_plotter', 'Spectrum plotting'),
]

def check_python_version():
    """Check Python version compatibility"""
    print(" Checking Python Version...")
    version = sys.version_info
    print(f"   Python {version.major}.{version.minor}.{version.micro}")

    if version.major < 3 or (version.major == 3 and version.minor < 9):
        print("   [FAIL] Python 3.9+ required")
        return False
    else:
        print("   [ok] Python version compatible")
        return True

def check_external_dependencies():
    """Check external package dependencies"""
    print("\nChecking External Dependencies...")

    missing = []
    available = []

    for package, description in REQUIRED_DEPENDENCIES.items():
        try:
            module = importlib.import_module(package)
            version = getattr(module, '__version__', 'unknown')
            print(f"   [ok] {package} {version} - {description}")
            available.append(package)
        except ImportError:
            print(f"   [FAIL] {package} - {description} (MISSING)")
            missing.append(package)

    for package, description in OPTIONAL_DEPENDENCIES.items():
        try:
            module = importlib.import_module(package)
            version = getattr(module, '__version__', 'unknown')
            print(f"   [ok] {package} {version} - {description} (optional)")
        except ImportError:
            print(f"   [--] {package} - {description} (optional, not installed)")

    return missing, available

def check_internal_modules():
    """Check internal module imports"""
    print("\nChecking Internal Modules...")

    missing = []
    available = []

    for module_name, description in INTERNAL_MODULES:
        try:
            importlib.import_module(module_name)
            print(f"   [ok] {module_name} - {description}")
            available.append(module_name)
        # A module that calls sys.exit() at import time raises SystemExit, which
        # is a BaseException and would otherwise abort the whole report.
        except (ImportError, SystemExit) as e:
            print(f"   [FAIL] {module_name} - {description} (ERROR: {e})")
            missing.append(module_name)

    return missing, available

def check_main_application():
    """Check main application can be imported"""
    print("\nChecking Main Application...")

    try:
        from lunaNMR.gui.main_window import LunaNMRMainWindow  # noqa: F401
        print("   [ok] Main GUI application can be imported")
        return True
    except (ImportError, SystemExit) as e:
        print(f"   [FAIL] Main GUI import failed: {e}")
        return False

def check_fitting_pipeline():
    """Check the integrator wires up a Voigt fitter with its quality methods"""
    print("\nChecking Fitting Pipeline...")

    try:
        from lunaNMR.core.core_integrator import EnhancedVoigtIntegrator
        integrator = EnhancedVoigtIntegrator()

        fitter = getattr(integrator, 'enhanced_fitter', None)
        if fitter is None:
            print("   [FAIL] Voigt fitter not initialized")
            return False
        print("   [ok] Voigt fitter initialized")

        required_methods = ('extract_local_peak_region', 'comprehensive_quality_assessment')
        for method in required_methods:
            if not hasattr(fitter, method):
                print(f"   [FAIL] {method} missing from fitter")
                return False
            print(f"   [ok] {method} available")

        return True

    except (Exception, SystemExit) as e:
        print(f"   [FAIL] Fitting pipeline check failed: {e}")
        return False

def generate_installation_report(missing_external, missing_internal, app_ok, fix_ok, python_ok):
    """Generate installation report"""
    print("\nINSTALLATION REPORT")
    print("=" * 50)

    if python_ok and not missing_external and not missing_internal and app_ok and fix_ok:
        print("INSTALLATION COMPLETE!")
        print("   [ok] All dependencies installed")
        print("   [ok] All modules available")
        print("   [ok] Main application ready")
        print("   [ok] Fitting pipeline ready")
        print(f"\nReady to run: python {LAUNCH_SCRIPT}")
        return True
    else:
        print("INSTALLATION INCOMPLETE!")

        if not python_ok:
            print("\n Python version is below the 3.9 minimum")

        if missing_external:
            print("\nMissing External Dependencies:")
            for pkg in missing_external:
                print(f"   - {pkg}")
            print(f"\nInstall with: pip install {' '.join(missing_external)}")
            print("   Or use: pip install -r requirements.txt")

        if missing_internal:
            print("\nMissing Internal Modules:")
            for mod in missing_internal:
                print(f"   - {mod}")
            print("   This indicates a problem with the installation")

        if not app_ok:
            print("\nMain application cannot be imported")
            print("   Check for import errors in lunaNMR.gui.main_window")

        if not fix_ok:
            print("\nFitting pipeline unavailable")
            print("   EnhancedVoigtIntegrator did not expose a usable Voigt fitter")

        return False

def main():
    """Main verification function"""
    print("NMR Peak Series Analysis - Installation Verification")
    print("=" * 65)

    python_ok = check_python_version()
    missing_external, _ = check_external_dependencies()
    missing_internal, _ = check_internal_modules()
    app_ok = check_main_application()
    fix_ok = check_fitting_pipeline()

    return generate_installation_report(
        missing_external, missing_internal, app_ok, fix_ok, python_ok)

if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)

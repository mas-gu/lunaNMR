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
    'sklearn': 'Machine learning (peak clustering)',
    'networkx': 'Graph clustering for peak detection',
    'nmrglue': 'NMR data file reading',
    'PySide6': 'Qt GUI framework',
    'numba': '2D multi-peak fitting acceleration',
}

# Absent installs degrade a feature rather than breaking the application.
OPTIONAL_DEPENDENCIES = {
    'torch': 'CNN peak classifier',
    'torchvision': 'CNN peak classifier',
    'psutil': 'Memory and performance monitoring',
}

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

    if version.major < 3 or (version.major == 3 and version.minor < 7):
        print("   ❌ Python 3.7+ required")
        return False
    else:
        print("   ✅ Python version compatible")
        return True

def check_external_dependencies():
    """Check external package dependencies"""
    print("\n📦 Checking External Dependencies...")

    missing = []
    available = []

    for package, description in REQUIRED_DEPENDENCIES.items():
        try:
            module = importlib.import_module(package)
            version = getattr(module, '__version__', 'unknown')
            print(f"   ✅ {package} {version} - {description}")
            available.append(package)
        except ImportError:
            print(f"   ❌ {package} - {description} (MISSING)")
            missing.append(package)

    for package, description in OPTIONAL_DEPENDENCIES.items():
        try:
            module = importlib.import_module(package)
            version = getattr(module, '__version__', 'unknown')
            print(f"   ✅ {package} {version} - {description} (optional)")
        except ImportError:
            print(f"   ➖ {package} - {description} (optional, not installed)")

    return missing, available

def check_internal_modules():
    """Check internal module imports"""
    print("\n🔧 Checking Internal Modules...")

    missing = []
    available = []

    for module_name, description in INTERNAL_MODULES:
        try:
            importlib.import_module(module_name)
            print(f"   ✅ {module_name} - {description}")
            available.append(module_name)
        # A module that calls sys.exit() at import time raises SystemExit, which
        # is a BaseException and would otherwise abort the whole report.
        except (ImportError, SystemExit) as e:
            print(f"   ❌ {module_name} - {description} (ERROR: {e})")
            missing.append(module_name)

    return missing, available

def check_main_application():
    """Check main application can be imported"""
    print("\n🎯 Checking Main Application...")

    try:
        from lunaNMR.gui.main_window import LunaNMRMainWindow  # noqa: F401
        print("   ✅ Main GUI application can be imported")
        return True
    except (ImportError, SystemExit) as e:
        print(f"   ❌ Main GUI import failed: {e}")
        return False

def check_fitting_pipeline():
    """Check the integrator wires up a Voigt fitter with its quality methods"""
    print("\n🔍 Checking Fitting Pipeline...")

    try:
        from lunaNMR.core.core_integrator import EnhancedVoigtIntegrator
        integrator = EnhancedVoigtIntegrator()

        fitter = getattr(integrator, 'enhanced_fitter', None)
        if fitter is None:
            print("   ❌ Voigt fitter not initialized")
            return False
        print("   ✅ Voigt fitter initialized")

        required_methods = ('extract_local_peak_region', 'comprehensive_quality_assessment')
        for method in required_methods:
            if not hasattr(fitter, method):
                print(f"   ❌ {method} missing from fitter")
                return False
            print(f"   ✅ {method} available")

        return True

    except (Exception, SystemExit) as e:
        print(f"   ❌ Fitting pipeline check failed: {e}")
        return False

def generate_installation_report(missing_external, missing_internal, app_ok, fix_ok, python_ok):
    """Generate installation report"""
    print("\n📋 INSTALLATION REPORT")
    print("=" * 50)

    if python_ok and not missing_external and not missing_internal and app_ok and fix_ok:
        print("🎉 INSTALLATION COMPLETE!")
        print("   ✅ All dependencies installed")
        print("   ✅ All modules available")
        print("   ✅ Main application ready")
        print("   ✅ Fitting pipeline ready")
        print(f"\n🚀 Ready to run: python {LAUNCH_SCRIPT}")
        return True
    else:
        print("❌ INSTALLATION INCOMPLETE!")

        if not python_ok:
            print("\n Python version is below the 3.7 minimum")

        if missing_external:
            print("\n📦 Missing External Dependencies:")
            for pkg in missing_external:
                print(f"   - {pkg}")
            print(f"\n💡 Install with: pip install {' '.join(missing_external)}")
            print("   Or use: pip install -r requirements.txt")

        if missing_internal:
            print("\n🔧 Missing Internal Modules:")
            for mod in missing_internal:
                print(f"   - {mod}")
            print("   This indicates a problem with the installation")

        if not app_ok:
            print("\n🎯 Main application cannot be imported")
            print("   Check for import errors in lunaNMR.gui.main_window")

        if not fix_ok:
            print("\n🔍 Fitting pipeline unavailable")
            print("   EnhancedVoigtIntegrator did not expose a usable Voigt fitter")

        return False

def main():
    """Main verification function"""
    print("🔍 NMR Peak Series Analysis - Installation Verification")
    print("=" * 65)

    # Run all checks
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

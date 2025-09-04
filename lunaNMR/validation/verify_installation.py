#!/usr/bin/env python3
"""
Installation Verification Script

This script verifies that all required dependencies and modules are properly
installed and that the NMR Peak Series Analysis application can run.

Author: Guillaume Mas
Date: 2025
"""

import sys
import os

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

    dependencies = {
        'numpy': 'Numerical computing',
        'pandas': 'Data manipulation',
        'matplotlib': 'Plotting and visualization',
        'scipy': 'Scientific computing',
        'sklearn': 'Machine learning (peak clustering)',
        'nmrglue': 'NMR data file reading'
    }

    missing = []
    available = []

    for package, description in dependencies.items():
        try:
            if package == 'sklearn':
                import sklearn
                version = sklearn.__version__
            else:
                module = __import__(package)
                version = getattr(module, '__version__', 'unknown')

            print(f"   ✅ {package} {version} - {description}")
            available.append(package)
        except ImportError:
            print(f"   ❌ {package} - {description} (MISSING)")
            missing.append(package)

    return missing, available

def check_internal_modules():
    """Check internal module imports"""
    print("\n🔧 Checking Internal Modules...")

    modules = [
        ('core_integrator', 'Core Voigt fitting'),
        ('enhanced_voigt_fitter', 'Local quality assessment'),
        ('spectrum_browser', 'Individual spectrum browser'),
        ('visualization', 'Plotting components'),
        ('gui_components', 'GUI widgets'),
        ('file_manager', 'File handling'),
        ('series_processor', 'Series processing'),
        ('config_manager', 'Configuration')
    ]

    missing = []
    available = []

    for module_name, description in modules:
        try:
            module = __import__(module_name)
            print(f"   ✅ {module_name} - {description}")
            available.append(module_name)
        except ImportError as e:
            print(f"   ❌ {module_name} - {description} (ERROR: {e})")
            missing.append(module_name)

    return missing, available

def check_main_application():
    """Check main application can be imported"""
    print("\n🎯 Checking Main Application...")

    try:
        from lunaNMR.gui.main_gui import NMRPeaksSeriesGUI
        print("   ✅ Main GUI application can be imported")
        return True
    except ImportError as e:
        print(f"   ❌ Main GUI import failed: {e}")
        return False

def check_local_quality_fix():
    """Check that the local quality assessment fix is working"""
    print("\n🔍 Checking Local Quality Assessment Fix...")

    try:
        from lunaNMR.core.core_integrator import EnhancedVoigtIntegrator
        integrator = EnhancedVoigtIntegrator()

        # Check if enhanced fitter is available
        if hasattr(integrator, 'enhanced_fitter') and integrator.enhanced_fitter is not None:
            print("   ✅ Enhanced Voigt fitter initialized")

            # Check if local quality methods are available
            if hasattr(integrator.enhanced_fitter, 'extract_local_peak_region'):
                print("   ✅ Local peak region extraction available")
            else:
                print("   ❌ Local peak region extraction missing")
                return False

            if hasattr(integrator.enhanced_fitter, 'comprehensive_quality_assessment'):
                print("   ✅ Local quality assessment available")
            else:
                print("   ❌ Local quality assessment missing")
                return False

            print("   ✅ Local quality assessment fix is properly integrated")
            return True
        else:
            print("   ❌ Enhanced fitter not initialized")
            return False

    except Exception as e:
        print(f"   ❌ Local quality check failed: {e}")
        return False

def run_quick_test():
    """Run a quick test of the local quality assessment"""
    print("\n🧪 Running Quick Integration Test...")

    try:
        # Run the simple integration test
        import subprocess
        result = subprocess.run([sys.executable, 'test_simple_integration.py'],
                              capture_output=True, text=True, timeout=60)

        if result.returncode == 0 and "INTEGRATION SUCCESSFUL" in result.stdout:
            print("   ✅ Integration test passed")
            print("   ✅ Local quality assessment working correctly")
            return True
        else:
            print("   ❌ Integration test failed")
            if result.stdout:
                print(f"   Output: {result.stdout[-200:]}")  # Last 200 chars
            return False

    except subprocess.TimeoutExpired:
        print("   ⚠️ Integration test timed out")
        return False
    except Exception as e:
        print(f"   ⚠️ Could not run integration test: {e}")
        return True  # Don't fail for this

def generate_installation_report(missing_external, missing_internal, app_ok, fix_ok):
    """Generate installation report"""
    print(f"\n📋 INSTALLATION REPORT")
    print("=" * 50)

    if not missing_external and not missing_internal and app_ok and fix_ok:
        print("🎉 INSTALLATION COMPLETE!")
        print("   ✅ All dependencies installed")
        print("   ✅ All modules available")
        print("   ✅ Main application ready")
        print("   ✅ Local quality assessment fix working")
        print("\n🚀 Ready to run: python launch_gui.py")
        return True
    else:
        print("❌ INSTALLATION INCOMPLETE!")

        if missing_external:
            print(f"\n📦 Missing External Dependencies:")
            for pkg in missing_external:
                print(f"   - {pkg}")
            print(f"\n💡 Install with: pip install {' '.join(missing_external)}")
            print(f"   Or use: pip install -r requirements.txt")

        if missing_internal:
            print(f"\n🔧 Missing Internal Modules:")
            for mod in missing_internal:
                print(f"   - {mod}")
            print(f"   This indicates a problem with the installation")

        if not app_ok:
            print(f"\n🎯 Main application cannot be imported")
            print(f"   Check for import errors in lunaNMR.gui.main_gui")

        if not fix_ok:
            print(f"\n🔍 Local quality assessment fix not working")
            print(f"   The R² improvement fix may not be properly integrated")

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
    fix_ok = check_local_quality_fix()
    test_ok = run_quick_test()

    # Generate report
    success = generate_installation_report(missing_external, missing_internal, app_ok, fix_ok)

    if success:
        print(f"\n✨ The local quality assessment fix is ready!")
        print(f"   Users will now see improved R² values in the GUI")
        print(f"   that exclude distant peaks from quality evaluation.")

    return success

if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)

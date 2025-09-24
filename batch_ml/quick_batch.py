#!/usr/bin/env python3
"""
Quick Batch Processing Script - Minimal Setup

This is the simplest possible script for batch processing NMR spectra.
It provides a minimal interface with sensible defaults for quick processing.

Usage:
    python quick_batch.py /path/to/spectra
    python quick_batch.py /path/to/spectra 15N1H  # specify nucleus type
    python quick_batch.py /path/to/spectra auto   # enable optimization

Example for your 100 spectra folder:
    python quick_batch.py /path/to/your/100/spectra auto
"""

import sys
from pathlib import Path

# Add the lunaNMR path to Python path (adjusted for new location)
script_dir = Path(__file__).parent
lunanmr_path = script_dir.parent  # Go up to lunaNMR_v0o9 folder
sys.path.insert(0, str(lunanmr_path))

def print_usage():
    """Print usage information."""
    print("Quick Batch NMR Processing")
    print("=" * 40)
    print("Usage:")
    print("  python quick_batch.py FOLDER")
    print("  python quick_batch.py FOLDER NUCLEUS")
    print("  python quick_batch.py FOLDER auto")
    print()
    print("Arguments:")
    print("  FOLDER    - Path to folder containing NMR spectra")
    print("  NUCLEUS   - Optional: 15N1H, 13C1H, 1H (auto-detect if not specified)")
    print("  auto      - Enable automatic parameter optimization")
    print()
    print("Examples:")
    print("  python quick_batch.py /data/spectra")
    print("  python quick_batch.py /data/spectra 15N1H")
    print("  python quick_batch.py /data/spectra auto")

def main():
    """Main entry point for quick batch processing."""
    # Check arguments
    if len(sys.argv) < 2 or sys.argv[1] in ['-h', '--help']:
        print_usage()
        return 0

    folder = sys.argv[1]

    # Parse second argument
    nucleus_type = None
    optimize = False

    if len(sys.argv) > 2:
        arg2 = sys.argv[2]
        if arg2 == 'auto':
            optimize = True
        elif arg2 in ['15N1H', '13C1H', '1H']:
            nucleus_type = arg2
        else:
            print(f"Warning: Unknown second argument '{arg2}', ignoring")

    # Check folder exists
    if not Path(folder).exists():
        print(f"Error: Folder '{folder}' does not exist")
        return 1

    try:
        from lunaNMR.batch_processing import BatchProcessor

        print(f"🚀 Starting batch processing...")
        print(f"   Folder: {folder}")
        print(f"   Nucleus: {nucleus_type or 'auto-detect'}")
        print(f"   Optimization: {'enabled' if optimize else 'disabled'}")
        print()

        # Create processor with sensible defaults
        config = {
            'processing': {
                'skip_on_error': True,  # Continue on errors
                'quality_threshold': 0.8,
                'max_fitting_attempts': 3
            },
            'logging': {
                'level': 'INFO',
                'detailed_logging': False,
                'progress_interval': 10  # Report every 10 files
            }
        }

        processor = BatchProcessor(config)

        # Process folder
        results = processor.process_folder(
            folder_path=folder,
            nucleus_type=nucleus_type,
            auto_optimize=optimize
        )

        # Print results
        print("\n" + "="*50)
        print("✅ BATCH PROCESSING COMPLETED")
        print("="*50)
        print(f"📁 Total files: {results['total_files']}")
        print(f"✅ Successfully processed: {results['processed_files']}")
        print(f"❌ Failed: {results['failed_files']}")
        print(f"🧬 Peaks detected: {results['total_peaks_detected']}")
        print(f"📊 ML samples collected: {results['total_ml_samples']}")

        success_rate = (results['processed_files'] / results['total_files'] * 100) if results['total_files'] > 0 else 0
        print(f"📈 Success rate: {success_rate:.1f}%")

        if results['total_ml_samples'] > 0:
            print(f"\n🎉 SUCCESS: {results['total_ml_samples']} ML training samples collected!")
            print("   Your training dataset is ready for Phase 2 ML model development")
        else:
            print("\n⚠️  Warning: No ML samples collected - check spectrum quality")

        return 0

    except ImportError as e:
        print(f"❌ Error: Could not import batch processing components")
        print(f"   Details: {e}")
        print("\n   Please ensure lunaNMR is properly set up")
        return 1
    except Exception as e:
        print(f"❌ Error during processing: {e}")
        return 1

if __name__ == "__main__":
    sys.exit(main())
#!/usr/bin/env python3
# ABOUTME: CLI script to extract Bruker 2D processed data (2rr) from TopSpin directories.
# ABOUTME: Copies 2rr + procs + proc2s to flat output structure for use with lunaNMR.
"""
Extract Bruker 2D Processed Data

This script scans a Bruker TopSpin dataset directory and extracts all 2D processed
spectra (experiments with 2rr files) to a flat output directory structure.

Usage:
    python extract_bruker_2d.py /path/to/bruker/dataset -o /path/to/output

Example:
    python extract_bruker_2d.py 20171009_CH1_domain -o extracted_spectra

    Input structure:
        20171009_CH1_domain/
        ├── 1/pdata/1/2rr, procs, proc2s
        ├── 2/pdata/1/2rr, procs, proc2s
        └── ...

    Output structure:
        extracted_spectra/
        ├── 1/2rr, procs, proc2s
        ├── 2/2rr, procs, proc2s
        └── ...
"""

import argparse
import os
import shutil
import sys
from pathlib import Path


def find_2d_experiments(bruker_dir: Path) -> list[tuple[str, Path]]:
    """Find all experiments with 2rr files in pdata/1/.

    Args:
        bruker_dir: Path to Bruker dataset directory

    Returns:
        List of (experiment_number, pdata_path) tuples, sorted by experiment number
    """
    experiments = []

    for item in bruker_dir.iterdir():
        if not item.is_dir():
            continue

        # Check if this looks like an experiment number (numeric folder name)
        try:
            exp_num = int(item.name)
        except ValueError:
            continue

        # Check for 2rr in pdata/1/
        pdata_path = item / "pdata" / "1"
        if (pdata_path / "2rr").exists():
            experiments.append((item.name, pdata_path))

    # Sort by experiment number
    experiments.sort(key=lambda x: int(x[0]))
    return experiments


def extract_experiment(exp_name: str, pdata_path: Path, output_dir: Path,
                       dry_run: bool = False, verbose: bool = False) -> bool:
    """Extract a single experiment's 2D data.

    Args:
        exp_name: Experiment number/name
        pdata_path: Path to pdata/1/ directory
        output_dir: Base output directory
        dry_run: If True, don't actually copy files
        verbose: If True, print detailed info

    Returns:
        True if successful, False otherwise
    """
    # Required files for lunaNMR
    required_files = ["2rr", "procs", "proc2s"]

    # Check all required files exist
    missing = [f for f in required_files if not (pdata_path / f).exists()]
    if missing:
        print(f"  Warning: Experiment {exp_name} missing files: {missing}")
        return False

    # Create output directory
    exp_output = output_dir / exp_name

    if dry_run:
        print(f"  [DRY RUN] Would create: {exp_output}")
        for f in required_files:
            print(f"    Copy: {pdata_path / f}")
        return True

    # Create directory and copy files
    exp_output.mkdir(parents=True, exist_ok=True)

    for f in required_files:
        src = pdata_path / f
        dst = exp_output / f
        shutil.copy2(src, dst)
        if verbose:
            print(f"    Copied: {f}")

    return True


def main():
    parser = argparse.ArgumentParser(
        description="Extract Bruker 2D processed data for lunaNMR",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s 20171009_CH1_domain -o extracted
  %(prog)s /data/bruker/T1_series -o ~/nmr/T1_spectra
  %(prog)s . -o output --dry-run  # Preview what would be extracted
        """
    )

    parser.add_argument(
        "bruker_dir",
        type=Path,
        help="Path to Bruker dataset directory (contains numbered experiment folders)"
    )

    parser.add_argument(
        "-o", "--output",
        type=Path,
        required=True,
        help="Output directory for extracted spectra"
    )

    parser.add_argument(
        "-n", "--dry-run",
        action="store_true",
        help="Show what would be extracted without copying files"
    )

    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="Print detailed progress information"
    )

    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite existing output directory"
    )

    args = parser.parse_args()

    # Validate input directory
    if not args.bruker_dir.exists():
        print(f"Error: Input directory does not exist: {args.bruker_dir}")
        sys.exit(1)

    if not args.bruker_dir.is_dir():
        print(f"Error: Input is not a directory: {args.bruker_dir}")
        sys.exit(1)

    # Check output directory
    if args.output.exists() and not args.dry_run:
        if not args.overwrite:
            print(f"Error: Output directory already exists: {args.output}")
            print("Use --overwrite to replace existing output")
            sys.exit(1)
        else:
            if args.verbose:
                print(f"Removing existing output directory: {args.output}")
            shutil.rmtree(args.output)

    # Find 2D experiments
    print(f"Scanning: {args.bruker_dir.resolve()}")
    experiments = find_2d_experiments(args.bruker_dir)

    if not experiments:
        print("No 2D experiments found (no 2rr files in pdata/1/)")
        sys.exit(0)

    print(f"Found {len(experiments)} 2D experiment(s): {', '.join(e[0] for e in experiments)}")

    if args.dry_run:
        print(f"\n[DRY RUN] Would extract to: {args.output.resolve()}")
    else:
        print(f"\nExtracting to: {args.output.resolve()}")

    # Extract each experiment
    success_count = 0
    for exp_name, pdata_path in experiments:
        print(f"\nExperiment {exp_name}:")
        if extract_experiment(exp_name, pdata_path, args.output,
                             dry_run=args.dry_run, verbose=args.verbose):
            success_count += 1

    # Summary
    print(f"\n{'=' * 40}")
    if args.dry_run:
        print(f"[DRY RUN] Would extract {success_count}/{len(experiments)} experiments")
    else:
        print(f"Successfully extracted {success_count}/{len(experiments)} experiments")
        print(f"Output: {args.output.resolve()}")


if __name__ == "__main__":
    main()

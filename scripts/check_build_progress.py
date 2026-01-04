# ABOUTME: Monitors progress of build_training_data.py by comparing checkpoint vs data folder.
# ABOUTME: Reports percentage completion based on spectrum/peaklist file pairs.

"""
Progress checker for build_training_data.py

Compares build_checkpoint.json against the input data folder to calculate
how far along the training data collection is.

Usage:
    python scripts/check_build_progress.py
    python scripts/check_build_progress.py --checkpoint path/to/checkpoint.json
    python scripts/check_build_progress.py --watch  # Continuous monitoring
"""

import argparse
import json
import sys
import time
from pathlib import Path
from typing import Set, Tuple


SPECTRUM_EXTENSIONS = {'.ft', '.ft2', '.ft1', '.2rr', '.2ii', '.fid', '.ucsf', '.pipe'}
PEAK_LIST_EXTENSIONS = {'.txt', '.csv'}


def find_spectrum_files(data_dir: Path) -> Set[Path]:
    """Find all spectrum files in directory tree."""
    spectrum_files = set()
    for ext in SPECTRUM_EXTENSIONS:
        spectrum_files.update(data_dir.glob(f'*{ext}'))
        spectrum_files.update(data_dir.glob(f'**/*{ext}'))
    return spectrum_files


def find_peak_list_for_spectrum(spectrum_path: Path) -> Path | None:
    """Find matching peak list for a spectrum file."""
    parent = spectrum_path.parent
    stem = spectrum_path.stem

    # Try specific patterns first
    candidates = [
        parent / f"{stem}_peaks.txt",
        parent / f"{stem}_peaklist.txt",
        parent / f"{stem}.txt",
        parent / f"{stem}_peaks.csv",
        parent / f"{stem}.csv",
    ]

    for candidate in candidates:
        if candidate.exists():
            return candidate

    # Check for any .txt file in same directory (common pattern)
    txt_files = list(parent.glob('*.txt'))
    if len(txt_files) == 1:
        return txt_files[0]

    return None


def count_valid_pairs(data_dir: Path) -> Tuple[int, int, Set[Path]]:
    """
    Count spectrum files that have matching peak lists.

    Returns: (total_spectra, valid_pairs, valid_spectrum_paths)
    """
    spectrum_files = find_spectrum_files(data_dir)
    valid_spectra = set()

    for spec in spectrum_files:
        if find_peak_list_for_spectrum(spec):
            valid_spectra.add(spec)

    return len(spectrum_files), len(valid_spectra), valid_spectra


def load_checkpoint(checkpoint_path: Path) -> dict:
    """Load checkpoint file."""
    if not checkpoint_path.exists():
        return {}

    with open(checkpoint_path, 'r') as f:
        return json.load(f)


def get_progress(checkpoint_path: Path) -> dict:
    """Calculate progress from checkpoint."""
    checkpoint = load_checkpoint(checkpoint_path)

    if not checkpoint:
        return {
            'status': 'no_checkpoint',
            'message': f'No checkpoint file found at {checkpoint_path}'
        }

    input_dir = Path(checkpoint.get('input_dir', ''))
    processed = checkpoint.get('processed', [])

    if not input_dir.exists():
        return {
            'status': 'error',
            'message': f'Input directory not found: {input_dir}'
        }

    # Count files
    total_spectra, valid_pairs, valid_paths = count_valid_pairs(input_dir)
    processed_count = len(processed)

    # Calculate percentage
    if valid_pairs > 0:
        percentage = (processed_count / valid_pairs) * 100
    else:
        percentage = 0.0

    # Find remaining files
    processed_set = {Path(p) for p in processed}
    remaining = valid_paths - processed_set

    return {
        'status': 'ok',
        'input_dir': str(input_dir),
        'total_spectrum_files': total_spectra,
        'valid_pairs': valid_pairs,
        'processed': processed_count,
        'remaining': len(remaining),
        'percentage': percentage,
        'timestamp': checkpoint.get('timestamp', 'unknown'),
        'mode': checkpoint.get('mode', 'unknown'),
    }


def format_progress(progress: dict, verbose: bool = False) -> str:
    """Format progress for display."""
    if progress['status'] != 'ok':
        return f"❌ {progress['message']}"

    bar_width = 40
    filled = int(bar_width * progress['percentage'] / 100)
    bar = '█' * filled + '░' * (bar_width - filled)

    lines = [
        f"",
        f"📊 Build Training Data Progress",
        f"{'─' * 50}",
        f"Input:     {progress['input_dir']}",
        f"Mode:      {progress['mode']}",
        f"",
        f"[{bar}] {progress['percentage']:.1f}%",
        f"",
        f"  Processed:  {progress['processed']:>5} files",
        f"  Remaining:  {progress['remaining']:>5} files",
        f"  Total:      {progress['valid_pairs']:>5} valid pairs",
        f"",
        f"Last update: {progress['timestamp']}",
    ]

    if verbose:
        lines.append(f"  (Total spectrum files: {progress['total_spectrum_files']})")

    return '\n'.join(lines)


def watch_progress(checkpoint_path: Path, interval: int = 10):
    """Continuously monitor progress."""
    print(f"Watching {checkpoint_path} (Ctrl+C to stop)\n")

    last_processed = -1

    try:
        while True:
            progress = get_progress(checkpoint_path)

            # Clear screen and print
            print('\033[2J\033[H', end='')  # Clear screen, move cursor to top
            print(format_progress(progress))

            if progress['status'] == 'ok':
                current = progress['processed']
                if last_processed >= 0 and current > last_processed:
                    rate = (current - last_processed) / interval
                    eta_files = progress['remaining']
                    if rate > 0:
                        eta_seconds = eta_files / rate
                        eta_min = int(eta_seconds / 60)
                        eta_sec = int(eta_seconds % 60)
                        print(f"\n  Rate: {rate:.1f} files/sec")
                        print(f"  ETA:  {eta_min}m {eta_sec}s")
                last_processed = current

            time.sleep(interval)

    except KeyboardInterrupt:
        print("\n\nStopped watching.")


def main():
    parser = argparse.ArgumentParser(
        description='Check progress of build_training_data.py',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    parser.add_argument(
        'checkpoint',
        nargs='?',
        type=Path,
        default=Path('ml_training_data/build_checkpoint.json'),
        help='Path to checkpoint file (default: ml_training_data/build_checkpoint.json)'
    )
    parser.add_argument(
        '--watch', '-w',
        action='store_true',
        help='Continuously monitor progress'
    )
    parser.add_argument(
        '--interval', '-i',
        type=int,
        default=10,
        help='Watch interval in seconds (default: 10)'
    )
    parser.add_argument(
        '--verbose', '-v',
        action='store_true',
        help='Show additional details'
    )
    parser.add_argument(
        '--json',
        action='store_true',
        help='Output as JSON'
    )

    args = parser.parse_args()

    if args.watch:
        watch_progress(args.checkpoint, args.interval)
    else:
        progress = get_progress(args.checkpoint)

        if args.json:
            print(json.dumps(progress, indent=2))
        else:
            print(format_progress(progress, args.verbose))

        return 0 if progress['status'] == 'ok' else 1


if __name__ == '__main__':
    sys.exit(main())

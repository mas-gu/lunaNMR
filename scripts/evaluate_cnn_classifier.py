#!/usr/bin/env python3
# ABOUTME: Evaluation script for CNN-based peak classifier.
# ABOUTME: Compares CNN detection against S/N baseline and RF models.

"""
Evaluate CNN Peak Classifier

This script evaluates a trained CNN peak classifier against expert peak lists.
Supports hybrid mode (S/N detection + CNN filter) for comparison with RF models.

Usage:
    python evaluate_cnn_classifier.py \
        --model lunaNMR/ml/pretrained/cnn_peak_classifier.pt \
        --spectra-dirs /path/to/spectra \
        --hybrid

Features:
- Evaluates on spectra with matching .txt peak lists
- Supports --test-only to evaluate only on held-out test spectra
- Hybrid mode: S/N detection → CNN filter (for direct comparison with RF hybrid)
- Computes precision, recall, F1
"""

import argparse
import csv
import logging
import sys
import time
from pathlib import Path
from typing import Dict, List, Tuple, Optional

import numpy as np

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

try:
    import torch
except ImportError:
    print("ERROR: PyTorch is required. Install with: pip install torch torchvision")
    sys.exit(1)

try:
    import nmrglue as ng
except ImportError:
    print("ERROR: nmrglue is required. Install with: pip install nmrglue")
    sys.exit(1)

from scipy.ndimage import maximum_filter

from lunaNMR.ml.cnn_peak_classifier import CNNPeakClassifier
from lunaNMR.core.core_integrator import EnhancedVoigtIntegrator

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def find_spectrum_peaklist_pairs(spectra_dirs: List[Path]) -> List[Tuple[Path, Path]]:
    """Find all spectrum files that have matching .txt peak list files."""
    nmr_extensions = ['.ft', '.ft2', '.ft3', '.pipe', '.ucsf']
    pairs = []

    for spectra_dir in spectra_dirs:
        for ext in nmr_extensions:
            for spectrum_path in spectra_dir.rglob(f"*{ext}"):
                txt_path = spectrum_path.with_suffix('.txt')
                if txt_path.exists():
                    pairs.append((spectrum_path, txt_path))

    return pairs


def load_peaklist_txt(txt_path: Path) -> List[Dict]:
    """Load peak list from .txt file (CSV format)."""
    peaks = []

    with open(txt_path, 'r') as f:
        reader = csv.reader(f)
        header = next(reader, None)

        for row in reader:
            if len(row) >= 3:
                try:
                    pos_f2 = float(row[1].strip())  # 1H
                    pos_f1 = float(row[2].strip())  # 15N/13C
                    peaks.append({
                        'assignment': row[0].strip(),
                        'pos_f1': pos_f1,
                        'pos_f2': pos_f2,
                    })
                except (ValueError, IndexError):
                    continue

    return peaks


def load_spectrum(spectrum_path: Path) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Load NMR spectrum file."""
    suffix = spectrum_path.suffix.lower()

    if suffix in ['.ft2', '.pipe', '.ft', '']:
        try:
            nmr_dict, nmr_data = ng.pipe.read(str(spectrum_path))
            uc_x = ng.pipe.make_uc(nmr_dict, nmr_data, dim=1)
            uc_y = ng.pipe.make_uc(nmr_dict, nmr_data, dim=0)
        except:
            nmr_dict, nmr_data = ng.bruker.read(str(spectrum_path.parent))
            udic = ng.bruker.guess_udic(nmr_dict, nmr_data, strip_fake=True)
            uc_x = ng.fileiobase.uc_from_udic(udic, dim=1)
            uc_y = ng.fileiobase.uc_from_udic(udic, dim=0)
    elif suffix == '.ucsf':
        nmr_dict, nmr_data = ng.sparky.read(str(spectrum_path))
        uc_x = ng.sparky.make_uc(nmr_dict, nmr_data, dim=1)
        uc_y = ng.sparky.make_uc(nmr_dict, nmr_data, dim=0)
    else:
        raise ValueError(f"Unknown format: {suffix}")

    ppm_x = uc_x.ppm_scale()
    ppm_y = uc_y.ppm_scale()

    return nmr_data, ppm_x, ppm_y


def ppm_to_index(ppm_axis: np.ndarray, ppm_value: float) -> int:
    """Convert PPM value to nearest index."""
    return int(np.argmin(np.abs(ppm_axis - ppm_value)))


def match_peaks(
    detected: List[Dict],
    reference: List[Dict],
    ppm_x: np.ndarray,
    ppm_y: np.ndarray,
    tolerance_f1: float = 0.3,
    tolerance_f2: float = 0.03,
) -> Tuple[int, int, int]:
    """Match detected peaks against reference peak list."""
    detected_ppm = []
    for d in detected:
        y_idx = d['y_idx']
        x_idx = d['x_idx']
        pos_f1 = ppm_y[y_idx]
        pos_f2 = ppm_x[x_idx]
        detected_ppm.append((pos_f1, pos_f2))

    ref_ppm = [(r['pos_f1'], r['pos_f2']) for r in reference]

    matched_detected = set()
    matched_ref = set()

    for i, (ref_f1, ref_f2) in enumerate(ref_ppm):
        best_j = None
        best_dist = float('inf')

        for j, (det_f1, det_f2) in enumerate(detected_ppm):
            if j in matched_detected:
                continue

            d_f1 = abs(det_f1 - ref_f1) / tolerance_f1
            d_f2 = abs(det_f2 - ref_f2) / tolerance_f2

            if d_f1 <= 1.0 and d_f2 <= 1.0:
                dist = np.sqrt(d_f1**2 + d_f2**2)
                if dist < best_dist:
                    best_dist = dist
                    best_j = j

        if best_j is not None:
            matched_ref.add(i)
            matched_detected.add(best_j)

    tp = len(matched_ref)
    fp = len(detected_ppm) - len(matched_detected)
    fn = len(ref_ppm) - len(matched_ref)

    return tp, fp, fn


def match_peaks_lunaNMR(
    detected: List[Dict],
    reference: List[Dict],
    tolerance_f1: float = 0.3,
    tolerance_f2: float = 0.03,
) -> Tuple[int, int, int]:
    """Match lunaNMR-detected peaks against reference peak list."""
    detected_ppm = [(d['ppm_y'], d['ppm_x']) for d in detected]
    ref_ppm = [(r['pos_f1'], r['pos_f2']) for r in reference]

    matched_detected = set()
    matched_ref = set()

    for i, (ref_f1, ref_f2) in enumerate(ref_ppm):
        best_j = None
        best_dist = float('inf')

        for j, (det_f1, det_f2) in enumerate(detected_ppm):
            if j in matched_detected:
                continue

            d_f1 = abs(det_f1 - ref_f1) / tolerance_f1
            d_f2 = abs(det_f2 - ref_f2) / tolerance_f2

            if d_f1 <= 1.0 and d_f2 <= 1.0:
                dist = np.sqrt(d_f1**2 + d_f2**2)
                if dist < best_dist:
                    best_dist = dist
                    best_j = j

        if best_j is not None:
            matched_ref.add(i)
            matched_detected.add(best_j)

    tp = len(matched_ref)
    fp = len(detected_ppm) - len(matched_detected)
    fn = len(ref_ppm) - len(matched_ref)

    return tp, fp, fn


def evaluate_cnn_standalone(
    classifier: CNNPeakClassifier,
    pairs: List[Tuple[Path, Path]],
    threshold: float = 0.5,
    min_snr: float = 3.0,
) -> Dict:
    """
    Evaluate CNN in standalone mode (local maxima → CNN classification).

    Args:
        classifier: Trained CNN classifier
        pairs: List of (spectrum_path, peaklist_path) tuples
        threshold: CNN classification threshold
        min_snr: Minimum S/N for candidate pre-filtering

    Returns:
        Metrics dict
    """
    total_tp = 0
    total_fp = 0
    total_fn = 0
    spectra_evaluated = 0
    total_time = 0.0

    for i, (spectrum_path, txt_path) in enumerate(pairs):
        reference_peaks = load_peaklist_txt(txt_path)
        if len(reference_peaks) < 3:
            continue

        try:
            spectrum, ppm_x, ppm_y = load_spectrum(spectrum_path)
        except Exception as e:
            logger.warning(f"Error loading {spectrum_path.name}: {e}")
            continue

        # Estimate noise level
        ny, nx = spectrum.shape
        corners = [
            spectrum[:20, :20],
            spectrum[:20, -20:],
            spectrum[-20:, :20],
            spectrum[-20:, -20:],
        ]
        noise_level = np.std(np.concatenate([c.flatten() for c in corners]))

        t_start = time.perf_counter()

        # Find local maxima candidates
        local_max = maximum_filter(spectrum, size=3)
        is_max = (spectrum == local_max) & (spectrum > min_snr * noise_level)
        y_indices, x_indices = np.where(is_max)

        candidates = [
            {'y_idx': int(y), 'x_idx': int(x)}
            for y, x in zip(y_indices, x_indices)
        ]

        if len(candidates) == 0:
            continue

        # CNN classification
        is_peak, confidence = classifier.classify_spectrum(
            spectrum, candidates, threshold=threshold
        )

        # Filter detected peaks
        detected = [c for c, p in zip(candidates, is_peak) if p]

        t_elapsed = time.perf_counter() - t_start
        total_time += t_elapsed

        # Match against reference
        tp, fp, fn = match_peaks(detected, reference_peaks, ppm_x, ppm_y)

        total_tp += tp
        total_fp += fp
        total_fn += fn
        spectra_evaluated += 1

        logger.info(f"[{i+1}/{len(pairs)}] {spectrum_path.name} (CNN): "
                   f"candidates={len(candidates)}, detected={len(detected)}, ref={len(reference_peaks)}, "
                   f"TP={tp}, FP={fp}, FN={fn}")

    precision = total_tp / (total_tp + total_fp) if (total_tp + total_fp) > 0 else 0
    recall = total_tp / (total_tp + total_fn) if (total_tp + total_fn) > 0 else 0
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0

    return {
        'method': f'CNN Standalone (threshold={threshold})',
        'threshold': threshold,
        'spectra_evaluated': spectra_evaluated,
        'total_tp': total_tp,
        'total_fp': total_fp,
        'total_fn': total_fn,
        'precision': precision,
        'recall': recall,
        'f1': f1,
        'total_time': total_time,
        'avg_time_per_spectrum': total_time / spectra_evaluated if spectra_evaluated > 0 else 0,
    }


def evaluate_cnn_hybrid(
    classifier: CNNPeakClassifier,
    pairs: List[Tuple[Path, Path]],
    sn_threshold: float = 10.0,
    expected_peak_count: int = 250,
    cnn_threshold: float = 0.3,
) -> Dict:
    """
    Evaluate CNN in hybrid mode (S/N detection → CNN filter).

    This is directly comparable to RF hybrid mode.

    Args:
        classifier: Trained CNN classifier
        pairs: List of (spectrum_path, peaklist_path) tuples
        sn_threshold: S/N threshold for initial detection
        expected_peak_count: Max peaks from S/N detection
        cnn_threshold: CNN threshold for filtering

    Returns:
        Metrics dict
    """
    total_tp = 0
    total_fp = 0
    total_fn = 0
    spectra_evaluated = 0
    total_time = 0.0
    total_sn_candidates = 0
    total_after_cnn = 0

    for i, (spectrum_path, txt_path) in enumerate(pairs):
        reference_peaks = load_peaklist_txt(txt_path)
        if len(reference_peaks) < 3:
            continue

        try:
            t_start = time.perf_counter()

            # Step 1: lunaNMR S/N detection
            integrator = EnhancedVoigtIntegrator()
            integrator.load_nmr_file(str(spectrum_path))
            integrator.sn_threshold = sn_threshold
            integrator.expected_peak_count = expected_peak_count

            sn_detected = integrator.detect_peaks_sn_native()

            if not sn_detected:
                continue

            total_sn_candidates += len(sn_detected)

            # Get spectrum data
            spectrum = integrator.nmr_data
            ppm_x = integrator.ppm_x_axis
            ppm_y = integrator.ppm_y_axis

            # Step 2: Convert to candidates for CNN
            candidates = []
            for peak in sn_detected:
                pos_f2 = peak['ppm_x']
                pos_f1 = peak['ppm_y']
                y_idx = ppm_to_index(ppm_y, pos_f1)
                x_idx = ppm_to_index(ppm_x, pos_f2)
                candidates.append({
                    'y_idx': y_idx,
                    'x_idx': x_idx,
                    'ppm_x': pos_f2,
                    'ppm_y': pos_f1,
                })

            # Step 3: CNN classification
            is_peak, confidence = classifier.classify_spectrum(
                spectrum, candidates, threshold=cnn_threshold
            )

            # Filter to keep only CNN-approved peaks
            filtered = [c for c, passed in zip(candidates, is_peak) if passed]

            t_elapsed = time.perf_counter() - t_start
            total_time += t_elapsed
            total_after_cnn += len(filtered)

        except Exception as e:
            logger.warning(f"Error with hybrid detection for {spectrum_path.name}: {e}")
            continue

        # Match against reference
        tp, fp, fn = match_peaks_lunaNMR(filtered, reference_peaks)

        total_tp += tp
        total_fp += fp
        total_fn += fn
        spectra_evaluated += 1

        logger.info(f"[{i+1}/{len(pairs)}] {spectrum_path.name} (Hybrid S/N={sn_threshold}+CNN={cnn_threshold}): "
                   f"S/N={len(sn_detected)}, after CNN={len(filtered)}, ref={len(reference_peaks)}, "
                   f"TP={tp}, FP={fp}, FN={fn}")

    precision = total_tp / (total_tp + total_fp) if (total_tp + total_fp) > 0 else 0
    recall = total_tp / (total_tp + total_fn) if (total_tp + total_fn) > 0 else 0
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0
    filter_rate = 1 - (total_after_cnn / total_sn_candidates) if total_sn_candidates > 0 else 0

    return {
        'method': f'Hybrid S/N+CNN (S/N={sn_threshold}, CNN={cnn_threshold})',
        'sn_threshold': sn_threshold,
        'cnn_threshold': cnn_threshold,
        'spectra_evaluated': spectra_evaluated,
        'total_tp': total_tp,
        'total_fp': total_fp,
        'total_fn': total_fn,
        'precision': precision,
        'recall': recall,
        'f1': f1,
        'total_time': total_time,
        'avg_time_per_spectrum': total_time / spectra_evaluated if spectra_evaluated > 0 else 0,
        'total_sn_candidates': total_sn_candidates,
        'total_after_cnn': total_after_cnn,
        'filter_rate': filter_rate,
    }


def evaluate_sn_baseline(
    pairs: List[Tuple[Path, Path]],
    sn_threshold: float = 10.0,
    expected_peak_count: int = 250,
) -> Dict:
    """Evaluate lunaNMR S/N detection as baseline."""
    total_tp = 0
    total_fp = 0
    total_fn = 0
    spectra_evaluated = 0
    total_time = 0.0

    for i, (spectrum_path, txt_path) in enumerate(pairs):
        reference_peaks = load_peaklist_txt(txt_path)
        if len(reference_peaks) < 3:
            continue

        try:
            integrator = EnhancedVoigtIntegrator()
            integrator.load_nmr_file(str(spectrum_path))
            integrator.sn_threshold = sn_threshold
            integrator.expected_peak_count = expected_peak_count

            t_start = time.perf_counter()
            detected = integrator.detect_peaks_sn_native()
            t_elapsed = time.perf_counter() - t_start
            total_time += t_elapsed

            if not detected:
                continue

        except Exception as e:
            logger.warning(f"Error with S/N detection for {spectrum_path.name}: {e}")
            continue

        tp, fp, fn = match_peaks_lunaNMR(detected, reference_peaks)

        total_tp += tp
        total_fp += fp
        total_fn += fn
        spectra_evaluated += 1

        logger.info(f"[{i+1}/{len(pairs)}] {spectrum_path.name} (S/N={sn_threshold}): "
                   f"detected={len(detected)}, ref={len(reference_peaks)}, "
                   f"TP={tp}, FP={fp}, FN={fn}")

    precision = total_tp / (total_tp + total_fp) if (total_tp + total_fp) > 0 else 0
    recall = total_tp / (total_tp + total_fn) if (total_tp + total_fn) > 0 else 0
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0

    return {
        'method': f'lunaNMR S/N (threshold={sn_threshold})',
        'sn_threshold': sn_threshold,
        'spectra_evaluated': spectra_evaluated,
        'total_tp': total_tp,
        'total_fp': total_fp,
        'total_fn': total_fn,
        'precision': precision,
        'recall': recall,
        'f1': f1,
        'total_time': total_time,
        'avg_time_per_spectrum': total_time / spectra_evaluated if spectra_evaluated > 0 else 0,
    }


def print_comparison(cnn_results: Dict, sn_results: Dict):
    """Print comparison table."""
    print("\n" + "=" * 85)
    print("CNN HYBRID vs STANDALONE S/N COMPARISON")
    print("=" * 85)
    print(f"{'Metric':<25} {'CNN Hybrid':<28} {'Standalone S/N':<20}")
    print("-" * 85)
    print(f"{'Spectra Evaluated':<25} {cnn_results['spectra_evaluated']:<28} {sn_results['spectra_evaluated']:<20}")
    print(f"{'True Positives':<25} {cnn_results['total_tp']:<28} {sn_results['total_tp']:<20}")
    print(f"{'False Positives':<25} {cnn_results['total_fp']:<28} {sn_results['total_fp']:<20}")
    print(f"{'False Negatives':<25} {cnn_results['total_fn']:<28} {sn_results['total_fn']:<20}")
    print("-" * 85)
    print(f"{'Precision':<25} {cnn_results['precision']:.4f}{'':<24} {sn_results['precision']:.4f}")
    print(f"{'Recall':<25} {cnn_results['recall']:.4f}{'':<24} {sn_results['recall']:.4f}")
    print(f"{'F1 Score':<25} {cnn_results['f1']:.4f}{'':<24} {sn_results['f1']:.4f}")
    print("-" * 85)

    if 'total_sn_candidates' in cnn_results:
        filter_pct = cnn_results['filter_rate'] * 100
        print(f"{'S/N Candidates':<25} {cnn_results['total_sn_candidates']:<28}")
        print(f"{'After CNN Filter':<25} {cnn_results['total_after_cnn']:<28}")
        print(f"{'Filter Rate':<25} {filter_pct:.1f}% rejected{'':<16}")
    print("-" * 85)

    cnn_time = cnn_results.get('avg_time_per_spectrum', 0)
    sn_time = sn_results.get('avg_time_per_spectrum', 0)
    print(f"{'Total Time (s)':<25} {cnn_results.get('total_time', 0):<28.2f} {sn_results.get('total_time', 0):<20.2f}")
    print(f"{'Avg Time/Spectrum (ms)':<25} {cnn_time*1000:<28.2f} {sn_time*1000:<20.2f}")
    print("=" * 85)

    # Analysis
    print("\nANALYSIS (CNN Hybrid vs Standalone S/N):")

    if cnn_results['precision'] > sn_results['precision']:
        diff = cnn_results['precision'] - sn_results['precision']
        print(f"  + CNN has {diff*100:.1f}% better precision")
    else:
        diff = sn_results['precision'] - cnn_results['precision']
        print(f"  - CNN has {diff*100:.1f}% worse precision")

    recall_loss = sn_results['recall'] - cnn_results['recall']
    if recall_loss < 0.01:
        print(f"  + Recall preserved ({cnn_results['recall']:.1%} vs {sn_results['recall']:.1%})")
    else:
        print(f"  ! Recall dropped {recall_loss*100:.1f}%")

    if cnn_results['f1'] > sn_results['f1']:
        diff = cnn_results['f1'] - sn_results['f1']
        print(f"  + CNN is {diff*100:.1f}% better overall (F1: {cnn_results['f1']:.4f} vs {sn_results['f1']:.4f})")
    else:
        diff = sn_results['f1'] - cnn_results['f1']
        print(f"  - CNN is {diff*100:.1f}% worse overall (F1: {cnn_results['f1']:.4f} vs {sn_results['f1']:.4f})")


def main():
    parser = argparse.ArgumentParser(description='Evaluate CNN peak classifier')
    parser.add_argument('--model', type=Path, required=True,
                       help='Path to trained CNN model (.pt)')
    parser.add_argument('--spectra-dirs', type=Path, nargs='+', required=True,
                       help='Directories to search for spectrum files with matching .txt peak lists')
    parser.add_argument('--threshold', type=float, default=0.5,
                       help='CNN classification threshold (default: 0.5)')
    parser.add_argument('--sn-threshold', type=float, default=10.0,
                       help='S/N threshold for baseline/hybrid detection (default: 10.0)')
    parser.add_argument('--expected-peaks', type=int, default=250,
                       help='Expected peak count for S/N detection (default: 250)')
    parser.add_argument('--test-only', action='store_true',
                       help='Only evaluate on test spectra from model metadata')
    parser.add_argument('--hybrid', action='store_true',
                       help='Evaluate hybrid mode (S/N + CNN filter)')
    parser.add_argument('--hybrid-threshold', type=float, default=0.3,
                       help='CNN threshold for hybrid filtering (default: 0.3)')

    args = parser.parse_args()

    # Load CNN model
    logger.info(f"Loading CNN model from {args.model}")
    classifier = CNNPeakClassifier(model_path=args.model)

    # Find spectrum/peak-list pairs
    logger.info("Searching for spectrum files with matching .txt peak lists...")
    pairs = find_spectrum_peaklist_pairs(args.spectra_dirs)
    logger.info(f"Found {len(pairs)} spectrum/peak-list pairs")

    if len(pairs) == 0:
        logger.error("No spectrum/peak-list pairs found!")
        sys.exit(1)

    # Filter to test-only if requested
    if args.test_only:
        checkpoint = torch.load(args.model, map_location='cpu')
        metadata = checkpoint.get('metadata', {})
        test_spectra = set(metadata.get('test_spectra', []))

        if test_spectra:
            original_count = len(pairs)
            pairs = [(sp, tp) for sp, tp in pairs
                     if sp.name in test_spectra or sp.stem in test_spectra]
            logger.info(f"Filtered to {len(pairs)} test-only spectra (from {original_count})")

            if len(pairs) == 0:
                logger.error("No test spectra found!")
                sys.exit(1)
        else:
            logger.warning("No test spectra in model metadata - evaluating on all")

    # Show found pairs
    for spectrum_path, txt_path in pairs[:5]:
        logger.info(f"  {spectrum_path.name} -> {txt_path.name}")
    if len(pairs) > 5:
        logger.info(f"  ... and {len(pairs) - 5} more")

    # Evaluate S/N baseline
    logger.info("\n" + "=" * 50)
    logger.info(f"Evaluating lunaNMR S/N Detection (S/N={args.sn_threshold})...")
    logger.info("=" * 50)
    sn_results = evaluate_sn_baseline(
        pairs,
        sn_threshold=args.sn_threshold,
        expected_peak_count=args.expected_peaks,
    )

    if args.hybrid:
        # Hybrid mode: S/N + CNN filter
        logger.info("\n" + "=" * 50)
        logger.info(f"Evaluating CNN Hybrid (S/N={args.sn_threshold} + CNN={args.hybrid_threshold})...")
        logger.info("=" * 50)
        cnn_results = evaluate_cnn_hybrid(
            classifier, pairs,
            sn_threshold=args.sn_threshold,
            expected_peak_count=args.expected_peaks,
            cnn_threshold=args.hybrid_threshold,
        )
    else:
        # Standalone CNN
        logger.info("\n" + "=" * 50)
        logger.info(f"Evaluating CNN Standalone (threshold={args.threshold})...")
        logger.info("=" * 50)
        cnn_results = evaluate_cnn_standalone(
            classifier, pairs,
            threshold=args.threshold,
        )

    # Print comparison
    print_comparison(cnn_results, sn_results)


if __name__ == '__main__':
    main()

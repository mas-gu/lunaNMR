#!/usr/bin/env python3
# ABOUTME: Script to evaluate trained peak classifier against expert peak lists.
# ABOUTME: Compares ML detection with S/N baseline and computes precision/recall metrics.

"""
Evaluate Peak Classifier

This script evaluates a trained peak classifier against expert-validated peak lists:
1. Loads trained model
2. For each test spectrum with a matching .txt peak list:
   - Runs ML peak detection (local maxima + classification)
   - Compares detected peaks vs expert peak list (.txt)
   - Computes precision, recall, F1
3. Compares performance vs traditional S/N threshold detection
4. Reports feature importance analysis

Usage:
    python evaluate_peak_classifier.py --model lunaNMR/ml/pretrained/peak_classifier.joblib \
                                       --spectra-dirs /path/to/spectra

The script finds all spectrum files (.ft, .ft2) that have matching .txt peak lists
in the same directory and uses those as ground truth.
"""

import argparse
import json
import logging
import sys
import time
from pathlib import Path
from typing import Dict, List, Tuple, Optional
import numpy as np

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

try:
    import nmrglue as ng
except ImportError:
    print("ERROR: nmrglue is required. Install with: pip install nmrglue")
    sys.exit(1)

from scipy.ndimage import maximum_filter
from lunaNMR.ml.peak_feature_extractor import PeakFeatureExtractor
from lunaNMR.ml.peak_classifier import PeakClassifier, MLPeakDetector
from lunaNMR.ml.feature_extractor import FeatureExtractor
from lunaNMR import EnhancedVoigtIntegrator
import csv
import joblib

# Nucleus detector instance for reuse
_nucleus_detector = FeatureExtractor()

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


class ModelWrapper:
    """Wrapper to handle both old (PeakClassifier) and new (v2.0.0+) model formats."""

    def __init__(self, model_path: Path):
        """Load model from file, auto-detecting format."""
        self.model_path = model_path
        model_data = joblib.load(model_path)

        if isinstance(model_data, dict) and 'version' in model_data:
            # New format (v2.0.0)
            self.scaler = model_data['scaler']
            self.model = model_data['model']
            self.model_type = model_data.get('model_type', 'rf')
            self.feature_names = model_data.get('feature_names', [])
            self.test_spectra = model_data.get('test_spectra', [])
            self.version = model_data.get('version', '2.0.0')
            logger.info(f"Loaded {self.model_type.upper()} model (v{self.version})")
        else:
            # Old format (PeakClassifier)
            self.scaler = model_data.get('scaler')
            self.model = model_data.get('model')
            self.model_type = 'rf'
            self.feature_names = model_data.get('feature_names', [])
            self.test_spectra = []
            self.version = '1.0.0'
            logger.info(f"Loaded legacy PeakClassifier model (v{self.version})")

    def classify(self, X: np.ndarray, threshold: float = 0.5) -> Tuple[np.ndarray, np.ndarray]:
        """Classify features, returning (is_peak, confidence)."""
        X_scaled = self.scaler.transform(X)
        probs = self.model.predict_proba(X_scaled)[:, 1]
        is_peak = probs >= threshold
        return is_peak, probs

    def get_feature_importance(self) -> Dict[str, float]:
        """Get feature importance from model."""
        importances = self.model.feature_importances_
        if self.feature_names:
            return dict(zip(self.feature_names, importances.tolist()))
        else:
            return {f"feature_{i}": imp for i, imp in enumerate(importances)}


class NucleusModelSelector:
    """
    Selector for nucleus-specific models.

    Loads and manages models for 15N, 13C, and optionally a fallback model.
    Automatically selects the appropriate model based on detected nucleus type.
    """

    def __init__(
        self,
        model_path: Optional[Path] = None,
        model_15N_path: Optional[Path] = None,
        model_13C_path: Optional[Path] = None,
    ):
        """
        Initialize with available model paths.

        Args:
            model_path: Fallback model for any nucleus type
            model_15N_path: Model trained specifically on 15N spectra
            model_13C_path: Model trained specifically on 13C spectra
        """
        self.models = {}
        self.fallback_model = None

        if model_15N_path:
            logger.info(f"Loading 15N-specific model: {model_15N_path}")
            self.models['15N'] = ModelWrapper(model_15N_path)

        if model_13C_path:
            logger.info(f"Loading 13C-specific model: {model_13C_path}")
            self.models['13C'] = ModelWrapper(model_13C_path)

        if model_path:
            logger.info(f"Loading fallback model: {model_path}")
            self.fallback_model = ModelWrapper(model_path)

        if not self.models and not self.fallback_model:
            raise ValueError("Must provide at least one model")

    def get_model(self, nucleus_type: str) -> ModelWrapper:
        """
        Get the appropriate model for the given nucleus type.

        Args:
            nucleus_type: '15N', '13C', or other

        Returns:
            ModelWrapper for the appropriate model
        """
        if nucleus_type in self.models:
            return self.models[nucleus_type]
        elif self.fallback_model:
            return self.fallback_model
        else:
            # Return any available model
            available = list(self.models.values())
            if available:
                logger.warning(f"No model for {nucleus_type}, using {list(self.models.keys())[0]} model")
                return available[0]
            raise ValueError(f"No model available for nucleus type: {nucleus_type}")

    def get_any_model(self) -> ModelWrapper:
        """Get any available model (for feature importance display)."""
        if self.fallback_model:
            return self.fallback_model
        return list(self.models.values())[0]


def find_spectrum_peaklist_pairs(spectra_dirs: List[Path]) -> List[Tuple[Path, Path]]:
    """
    Find all spectrum files that have matching .txt peak list files.

    Returns list of (spectrum_path, peaklist_path) tuples.
    """
    nmr_extensions = ['.ft', '.ft2', '.ft3', '.pipe', '.ucsf']
    pairs = []

    for spectra_dir in spectra_dirs:
        # Search recursively
        for ext in nmr_extensions:
            for spectrum_path in spectra_dir.rglob(f"*{ext}"):
                # Check for matching .txt file
                txt_path = spectrum_path.with_suffix('.txt')
                if txt_path.exists():
                    pairs.append((spectrum_path, txt_path))

    return pairs


def load_peaklist_txt(txt_path: Path) -> List[Dict]:
    """
    Load peak list from .txt file (CSV format).

    Expected format: Assignment, Position_X, Position_Y, Height
    Position_X = 1H (ppm), Position_Y = 15N/13C (ppm)

    Returns list of dicts with pos_f1 (Y) and pos_f2 (X).
    """
    peaks = []

    with open(txt_path, 'r') as f:
        reader = csv.reader(f)
        header = next(reader, None)  # Skip header

        for row in reader:
            if len(row) >= 3:
                try:
                    # Position_X = 1H = F2, Position_Y = 15N = F1
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


def estimate_noise_level(spectrum: np.ndarray) -> float:
    """Estimate noise level from spectrum corners."""
    ny, nx = spectrum.shape
    corner_size_y = max(10, ny // 20)
    corner_size_x = max(10, nx // 20)

    corners = [
        spectrum[:corner_size_y, :corner_size_x],
        spectrum[:corner_size_y, -corner_size_x:],
        spectrum[-corner_size_y:, :corner_size_x],
        spectrum[-corner_size_y:, -corner_size_x:],
    ]

    corner_data = np.concatenate([c.flatten() for c in corners])
    return float(np.std(corner_data))


def load_spectrum(spectrum_path: Path) -> Tuple[np.ndarray, np.ndarray, np.ndarray, Dict]:
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

    return nmr_data, ppm_x, ppm_y, nmr_dict


def find_spectrum_file(spectrum_name: str, spectra_dirs: List[Path]) -> Optional[Path]:
    """Find spectrum file by name in search directories."""
    # NMR file extensions (in priority order)
    nmr_extensions = ['.ft', '.ft2', '.ft3', '.pipe', '.ucsf']
    base_name = Path(spectrum_name).stem

    for spectra_dir in spectra_dirs:
        # First pass: try exact filename matches
        for ext in [''] + nmr_extensions:
            path = spectra_dir / (spectrum_name + ext)
            if path.exists() and path.suffix.lower() in nmr_extensions + ['']:
                return path
            path = spectra_dir / (base_name + ext)
            if path.exists() and path.suffix.lower() in nmr_extensions + ['']:
                return path

        # Second pass: recursive search with NMR extensions only
        for ext in nmr_extensions:
            matches = list(spectra_dir.rglob(f"{base_name}{ext}"))
            if matches:
                return matches[0]

    return None


def ppm_to_index(ppm_axis: np.ndarray, ppm_value: float) -> int:
    """Convert PPM value to nearest index."""
    return int(np.argmin(np.abs(ppm_axis - ppm_value)))


def find_local_maxima(
    spectrum: np.ndarray,
    min_snr: float = 0.0,
    noise_level: float = 1.0,
) -> List[Dict]:
    """
    Find all local maxima in spectrum.

    Returns list of candidates with y_idx, x_idx, intensity.
    """
    # Local maximum filter (3x3 neighborhood)
    local_max = maximum_filter(spectrum, size=3)
    is_max = (spectrum == local_max) & (spectrum > 0)

    # Optional S/N filter
    if min_snr > 0:
        is_max = is_max & (spectrum > min_snr * noise_level)

    # Get coordinates
    y_indices, x_indices = np.where(is_max)

    candidates = []
    for y_idx, x_idx in zip(y_indices, x_indices):
        candidates.append({
            'y_idx': int(y_idx),
            'x_idx': int(x_idx),
            'intensity': float(spectrum[y_idx, x_idx]),
        })

    return candidates


def match_peaks(
    detected: List[Dict],
    reference: List[Dict],
    ppm_x: np.ndarray,
    ppm_y: np.ndarray,
    tolerance_f1: float = 0.3,
    tolerance_f2: float = 0.03,
) -> Tuple[int, int, int]:
    """
    Match detected peaks against reference peak list.

    Args:
        detected: List of detected peaks (with y_idx, x_idx)
        reference: List of reference peaks (with pos_f1, pos_f2)
        ppm_x, ppm_y: PPM axes
        tolerance_f1: Matching tolerance in F1 (ppm)
        tolerance_f2: Matching tolerance in F2 (ppm)

    Returns:
        Tuple of (true_positives, false_positives, false_negatives)
    """
    # Convert detected to PPM coordinates
    detected_ppm = []
    for d in detected:
        y_idx = d['y_idx']
        x_idx = d['x_idx']
        pos_f1 = ppm_y[y_idx]
        pos_f2 = ppm_x[x_idx]
        detected_ppm.append((pos_f1, pos_f2))

    # Reference PPM coordinates
    ref_ppm = [(r['pos_f1'], r['pos_f2']) for r in reference]

    # Match each reference peak to nearest detected peak
    matched_detected = set()
    matched_ref = set()

    for i, (ref_f1, ref_f2) in enumerate(ref_ppm):
        best_j = None
        best_dist = float('inf')

        for j, (det_f1, det_f2) in enumerate(detected_ppm):
            if j in matched_detected:
                continue

            # Normalized distance (scaled by tolerance)
            d_f1 = abs(det_f1 - ref_f1) / tolerance_f1
            d_f2 = abs(det_f2 - ref_f2) / tolerance_f2

            # Within tolerance?
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


def evaluate_ml_detector(
    model_path: Path,
    test_data: Dict,
    spectra_dirs: List[Path],
    threshold: float = 0.4,
    min_snr_prefilter: float = 0.0,
    min_r2: float = 0.85,
) -> Dict:
    """
    Evaluate ML peak detector on test data.

    Args:
        threshold: ML classification confidence threshold
        min_snr_prefilter: Minimum S/N to consider as candidate (pre-filter)

    Returns metrics dict.
    """
    classifier = PeakClassifier.load(model_path)

    total_tp = 0
    total_fp = 0
    total_fn = 0
    spectra_evaluated = 0
    total_time = 0.0

    samples = test_data.get('samples', [])

    for i, sample in enumerate(samples):
        spectrum_name = sample.get('spectrum_name', '')
        if not spectrum_name:
            continue

        spectrum_path = find_spectrum_file(spectrum_name, spectra_dirs)
        if spectrum_path is None:
            logger.warning(f"Spectrum not found: {spectrum_name}")
            continue

        try:
            spectrum, ppm_x, ppm_y, _ = load_spectrum(spectrum_path)
        except Exception as e:
            logger.warning(f"Error loading {spectrum_name}: {e}")
            continue

        noise_level = sample.get('noise_level', 1.0)
        nucleus_type = sample.get('nucleus_type', '15N')

        # Reference peaks (high quality only)
        peaks = sample.get('peaks', [])
        reference_peaks = [
            p for p in peaks
            if p.get('r_squared', 0) >= min_r2 and p.get('pos_f1', 0) > 0
        ]

        if len(reference_peaks) < 5:
            continue

        # ML detection with optional S/N pre-filter
        t_start = time.perf_counter()

        candidates = find_local_maxima(spectrum, min_snr=min_snr_prefilter, noise_level=noise_level)

        if len(candidates) == 0:
            continue

        # Extract features with nucleus-specific templates
        ppm_per_point_f1 = abs(ppm_y[1] - ppm_y[0]) if len(ppm_y) > 1 else 0.1
        ppm_per_point_f2 = abs(ppm_x[1] - ppm_x[0]) if len(ppm_x) > 1 else 0.01
        extractor = PeakFeatureExtractor(
            noise_level=noise_level,
            nucleus_type=nucleus_type,
            ppm_per_point_f1=ppm_per_point_f1,
            ppm_per_point_f2=ppm_per_point_f2,
        )
        features = extractor.extract_features_batch(spectrum, candidates)

        # Classify
        is_peak, confidence = classifier.classify(features, threshold)

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

        logger.info(f"[{i+1}/{len(samples)}] {spectrum_name} (ML, prefilter={min_snr_prefilter}): "
                   f"candidates={len(candidates)}, detected={len(detected)}, ref={len(reference_peaks)}, "
                   f"TP={tp}, FP={fp}, FN={fn}")

    # Compute metrics
    precision = total_tp / (total_tp + total_fp) if (total_tp + total_fp) > 0 else 0
    recall = total_tp / (total_tp + total_fn) if (total_tp + total_fn) > 0 else 0
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0

    avg_time = total_time / spectra_evaluated if spectra_evaluated > 0 else 0

    return {
        'method': f'ML Detector (prefilter={min_snr_prefilter})',
        'threshold': threshold,
        'min_snr_prefilter': min_snr_prefilter,
        'spectra_evaluated': spectra_evaluated,
        'total_tp': total_tp,
        'total_fp': total_fp,
        'total_fn': total_fn,
        'precision': precision,
        'recall': recall,
        'f1': f1,
        'total_time': total_time,
        'avg_time_per_spectrum': avg_time,
    }


def match_peaks_lunaNMR(
    detected: List[Dict],
    reference: List[Dict],
    tolerance_f1: float = 0.3,
    tolerance_f2: float = 0.03,
) -> Tuple[int, int, int]:
    """
    Match lunaNMR-detected peaks against reference peak list.

    Args:
        detected: List of detected peaks from lunaNMR (with ppm_x, ppm_y)
        reference: List of reference peaks (with pos_f1, pos_f2)
        tolerance_f1: Matching tolerance in F1 (ppm)
        tolerance_f2: Matching tolerance in F2 (ppm)

    Returns:
        Tuple of (true_positives, false_positives, false_negatives)
    """
    # lunaNMR uses ppm_x (F2) and ppm_y (F1)
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


def evaluate_sn_detector(
    test_data: Dict,
    spectra_dirs: List[Path],
    sn_threshold: float = 10.0,
    expected_peak_count: int = 250,
    min_r2: float = 0.85,
) -> Dict:
    """
    Evaluate lunaNMR S/N threshold detection using actual detection routines.

    Uses EnhancedVoigtIntegrator.detect_peaks_sn_native() with GUI default
    parameters for fair comparison.

    Args:
        sn_threshold: S/N threshold multiplier (GUI default: 10.0)
        expected_peak_count: Max peaks to return (default: 250)

    Returns metrics dict.
    """
    total_tp = 0
    total_fp = 0
    total_fn = 0
    spectra_evaluated = 0
    total_time = 0.0

    samples = test_data.get('samples', [])

    for i, sample in enumerate(samples):
        spectrum_name = sample.get('spectrum_name', '')
        if not spectrum_name:
            continue

        spectrum_path = find_spectrum_file(spectrum_name, spectra_dirs)
        if spectrum_path is None:
            continue

        # Reference peaks
        peaks = sample.get('peaks', [])
        reference_peaks = [
            p for p in peaks
            if p.get('r_squared', 0) >= min_r2 and p.get('pos_f1', 0) > 0
        ]

        if len(reference_peaks) < 5:
            continue

        try:
            # Use lunaNMR's actual detection
            integrator = EnhancedVoigtIntegrator()
            integrator.load_nmr_file(str(spectrum_path))

            # Configure with GUI default parameters
            integrator.sn_threshold = sn_threshold
            integrator.expected_peak_count = expected_peak_count

            # Time only the detection step (not file loading)
            t_start = time.perf_counter()
            detected = integrator.detect_peaks_sn_native()
            t_elapsed = time.perf_counter() - t_start
            total_time += t_elapsed

            if not detected:
                continue

        except Exception as e:
            logger.warning(f"Error with lunaNMR detection for {spectrum_name}: {e}")
            continue

        # Match against reference
        tp, fp, fn = match_peaks_lunaNMR(detected, reference_peaks)

        total_tp += tp
        total_fp += fp
        total_fn += fn
        spectra_evaluated += 1

        logger.info(f"[{i+1}/{len(samples)}] {spectrum_name} (lunaNMR S/N={sn_threshold}): "
                   f"detected={len(detected)}, ref={len(reference_peaks)}, "
                   f"TP={tp}, FP={fp}, FN={fn}")

    precision = total_tp / (total_tp + total_fp) if (total_tp + total_fp) > 0 else 0
    recall = total_tp / (total_tp + total_fn) if (total_tp + total_fn) > 0 else 0
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0

    avg_time = total_time / spectra_evaluated if spectra_evaluated > 0 else 0

    return {
        'method': f'lunaNMR S/N (threshold={sn_threshold})',
        'sn_threshold': sn_threshold,
        'expected_peak_count': expected_peak_count,
        'spectra_evaluated': spectra_evaluated,
        'total_tp': total_tp,
        'total_fp': total_fp,
        'total_fn': total_fn,
        'precision': precision,
        'recall': recall,
        'f1': f1,
        'total_time': total_time,
        'avg_time_per_spectrum': avg_time,
    }


def print_comparison(ml_results: Dict, sn_results: Dict):
    """Print comparison table."""
    print("\n" + "=" * 70)
    print("DETECTION METHOD COMPARISON")
    print("=" * 70)
    print(f"{'Metric':<25} {'ML Detector':<20} {'lunaNMR S/N':<20}")
    print("-" * 70)
    print(f"{'Spectra Evaluated':<25} {ml_results['spectra_evaluated']:<20} {sn_results['spectra_evaluated']:<20}")
    print(f"{'True Positives':<25} {ml_results['total_tp']:<20} {sn_results['total_tp']:<20}")
    print(f"{'False Positives':<25} {ml_results['total_fp']:<20} {sn_results['total_fp']:<20}")
    print(f"{'False Negatives':<25} {ml_results['total_fn']:<20} {sn_results['total_fn']:<20}")
    print("-" * 70)
    print(f"{'Precision':<25} {ml_results['precision']:.4f}{'':<16} {sn_results['precision']:.4f}")
    print(f"{'Recall':<25} {ml_results['recall']:.4f}{'':<16} {sn_results['recall']:.4f}")
    print(f"{'F1 Score':<25} {ml_results['f1']:.4f}{'':<16} {sn_results['f1']:.4f}")
    print("-" * 70)
    # Timing info
    ml_time = ml_results.get('avg_time_per_spectrum', 0)
    sn_time = sn_results.get('avg_time_per_spectrum', 0)
    print(f"{'Total Time (s)':<25} {ml_results.get('total_time', 0):<20.2f} {sn_results.get('total_time', 0):<20.2f}")
    print(f"{'Avg Time/Spectrum (ms)':<25} {ml_time*1000:<20.2f} {sn_time*1000:<20.2f}")
    print("=" * 70)

    # Analysis
    print("\nANALYSIS:")
    if ml_results['recall'] > sn_results['recall']:
        diff = ml_results['recall'] - sn_results['recall']
        print(f"  ✓ ML catches {diff*100:.1f}% more peaks (better recall)")
    elif sn_results['recall'] > ml_results['recall']:
        diff = sn_results['recall'] - ml_results['recall']
        print(f"  ✗ S/N catches {diff*100:.1f}% more peaks")

    if ml_results['precision'] > sn_results['precision']:
        diff = ml_results['precision'] - sn_results['precision']
        print(f"  ✓ ML has {diff*100:.1f}% fewer false positives (better precision)")
    elif sn_results['precision'] > ml_results['precision']:
        diff = sn_results['precision'] - ml_results['precision']
        print(f"  ✗ S/N has {diff*100:.1f}% fewer false positives")

    if ml_results['f1'] > sn_results['f1']:
        diff = ml_results['f1'] - sn_results['f1']
        print(f"  ✓ ML is {diff*100:.1f}% better overall (F1 score)")
    elif sn_results['f1'] > ml_results['f1']:
        diff = sn_results['f1'] - ml_results['f1']
        print(f"  ✗ S/N is {diff*100:.1f}% better overall (F1 score)")

    # Speed comparison
    ml_time = ml_results.get('avg_time_per_spectrum', 0)
    sn_time = sn_results.get('avg_time_per_spectrum', 0)
    if ml_time > 0 and sn_time > 0:
        if ml_time < sn_time:
            ratio = sn_time / ml_time
            print(f"  ✓ ML is {ratio:.1f}x faster ({ml_time*1000:.1f}ms vs {sn_time*1000:.1f}ms per spectrum)")
        else:
            ratio = ml_time / sn_time
            print(f"  ✗ S/N is {ratio:.1f}x faster ({sn_time*1000:.1f}ms vs {ml_time*1000:.1f}ms per spectrum)")


def print_hybrid_comparison(hybrid_results: Dict, sn_results: Dict):
    """Print comparison table for hybrid vs S/N."""
    print("\n" + "=" * 85)
    print("HYBRID (S/N + ML FILTER) vs STANDALONE S/N COMPARISON")
    print("=" * 85)
    print(f"{'Metric':<25} {'Hybrid S/N+ML':<28} {'Standalone S/N':<20}")
    print("-" * 85)
    print(f"{'Spectra Evaluated':<25} {hybrid_results['spectra_evaluated']:<28} {sn_results['spectra_evaluated']:<20}")
    print(f"{'True Positives':<25} {hybrid_results['total_tp']:<28} {sn_results['total_tp']:<20}")
    print(f"{'False Positives':<25} {hybrid_results['total_fp']:<28} {sn_results['total_fp']:<20}")
    print(f"{'False Negatives':<25} {hybrid_results['total_fn']:<28} {sn_results['total_fn']:<20}")
    print("-" * 85)
    print(f"{'Precision':<25} {hybrid_results['precision']:.4f}{'':<24} {sn_results['precision']:.4f}")
    print(f"{'Recall':<25} {hybrid_results['recall']:.4f}{'':<24} {sn_results['recall']:.4f}")
    print(f"{'F1 Score':<25} {hybrid_results['f1']:.4f}{'':<24} {sn_results['f1']:.4f}")
    print("-" * 85)

    # Hybrid-specific stats
    if 'total_sn_candidates' in hybrid_results:
        filter_pct = hybrid_results['filter_rate'] * 100
        print(f"{'S/N Candidates':<25} {hybrid_results['total_sn_candidates']:<28}")
        print(f"{'After ML Filter':<25} {hybrid_results['total_after_ml']:<28}")
        print(f"{'Filter Rate':<25} {filter_pct:.1f}% rejected{'':<16}")
    print("-" * 85)

    # Timing
    hybrid_time = hybrid_results.get('avg_time_per_spectrum', 0)
    sn_time = sn_results.get('avg_time_per_spectrum', 0)
    print(f"{'Total Time (s)':<25} {hybrid_results.get('total_time', 0):<28.2f} {sn_results.get('total_time', 0):<20.2f}")
    print(f"{'Avg Time/Spectrum (ms)':<25} {hybrid_time*1000:<28.2f} {sn_time*1000:<20.2f}")
    print("=" * 85)

    # Analysis
    print("\nANALYSIS (Hybrid vs Standalone S/N):")

    # Precision improvement
    if hybrid_results['precision'] > sn_results['precision']:
        diff = hybrid_results['precision'] - sn_results['precision']
        print(f"  ✓ Hybrid has {diff*100:.1f}% better precision (fewer false positives)")
    else:
        diff = sn_results['precision'] - hybrid_results['precision']
        print(f"  ✗ Hybrid has {diff*100:.1f}% worse precision")

    # Recall preservation
    recall_loss = sn_results['recall'] - hybrid_results['recall']
    if recall_loss < 0.01:
        print(f"  ✓ Recall preserved ({hybrid_results['recall']:.1%} vs {sn_results['recall']:.1%})")
    else:
        print(f"  ⚠ Recall dropped {recall_loss*100:.1f}% ({hybrid_results['recall']:.1%} vs {sn_results['recall']:.1%})")

    # F1 comparison
    if hybrid_results['f1'] > sn_results['f1']:
        diff = hybrid_results['f1'] - sn_results['f1']
        print(f"  ✓ Hybrid is {diff*100:.1f}% better overall (F1: {hybrid_results['f1']:.4f} vs {sn_results['f1']:.4f})")
    else:
        diff = sn_results['f1'] - hybrid_results['f1']
        print(f"  ✗ Hybrid is {diff*100:.1f}% worse overall (F1: {hybrid_results['f1']:.4f} vs {sn_results['f1']:.4f})")

    # Speed
    if hybrid_time > 0 and sn_time > 0:
        overhead = hybrid_time - sn_time
        ratio = hybrid_time / sn_time
        print(f"  ⏱ ML filtering adds {overhead*1000:.0f}ms overhead ({ratio:.1f}x slower)")


def evaluate_ml_detector_txt(
    model_path: Path,
    pairs: List[Tuple[Path, Path]],
    threshold: float = 0.4,
    min_snr_prefilter: float = 0.0,
) -> Tuple[Dict, ModelWrapper]:
    """
    Evaluate ML peak detector using .txt peak lists as ground truth.

    Args:
        model_path: Path to trained model
        pairs: List of (spectrum_path, peaklist_path) tuples
        threshold: ML classification confidence threshold
        min_snr_prefilter: Minimum S/N to consider as candidate (pre-filter)

    Returns:
        Tuple of (metrics dict, model wrapper)
    """
    classifier = ModelWrapper(model_path)

    total_tp = 0
    total_fp = 0
    total_fn = 0
    spectra_evaluated = 0
    total_time = 0.0

    for i, (spectrum_path, txt_path) in enumerate(pairs):
        try:
            spectrum, ppm_x, ppm_y, _ = load_spectrum(spectrum_path)
        except Exception as e:
            logger.warning(f"Error loading {spectrum_path.name}: {e}")
            continue

        # Load reference peaks from .txt
        reference_peaks = load_peaklist_txt(txt_path)
        if len(reference_peaks) < 3:
            logger.warning(f"Too few peaks in {txt_path.name}: {len(reference_peaks)}")
            continue

        noise_level = estimate_noise_level(spectrum)

        # Detect nucleus type using lunaNMR's standard method
        nucleus_type = _nucleus_detector._detect_nucleus_type(ppm_y)

        # ML detection with timing
        t_start = time.perf_counter()

        candidates = find_local_maxima(spectrum, min_snr=min_snr_prefilter, noise_level=noise_level)

        if len(candidates) == 0:
            logger.warning(f"No candidates in {spectrum_path.name}")
            continue

        # Extract features with nucleus-specific templates
        ppm_per_point_f1 = abs(ppm_y[1] - ppm_y[0]) if len(ppm_y) > 1 else 0.1
        ppm_per_point_f2 = abs(ppm_x[1] - ppm_x[0]) if len(ppm_x) > 1 else 0.01
        extractor = PeakFeatureExtractor(
            noise_level=noise_level,
            nucleus_type=nucleus_type,
            ppm_per_point_f1=ppm_per_point_f1,
            ppm_per_point_f2=ppm_per_point_f2,
        )
        features = extractor.extract_features_batch(spectrum, candidates)

        # Classify
        is_peak, confidence = classifier.classify(features, threshold)

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

        logger.info(f"[{i+1}/{len(pairs)}] {spectrum_path.name} (ML): "
                   f"candidates={len(candidates)}, detected={len(detected)}, ref={len(reference_peaks)}, "
                   f"TP={tp}, FP={fp}, FN={fn}, time={t_elapsed*1000:.1f}ms")

    # Compute metrics
    precision = total_tp / (total_tp + total_fp) if (total_tp + total_fp) > 0 else 0
    recall = total_tp / (total_tp + total_fn) if (total_tp + total_fn) > 0 else 0
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0

    avg_time = total_time / spectra_evaluated if spectra_evaluated > 0 else 0

    results = {
        'method': f'ML Detector (prefilter={min_snr_prefilter})',
        'threshold': threshold,
        'min_snr_prefilter': min_snr_prefilter,
        'spectra_evaluated': spectra_evaluated,
        'total_tp': total_tp,
        'total_fp': total_fp,
        'total_fn': total_fn,
        'precision': precision,
        'recall': recall,
        'f1': f1,
        'total_time': total_time,
        'avg_time_per_spectrum': avg_time,
    }
    return results, classifier


def evaluate_hybrid_detector_txt(
    model_selector: NucleusModelSelector,
    pairs: List[Tuple[Path, Path]],
    sn_threshold: float = 10.0,
    expected_peak_count: int = 250,
    ml_threshold: float = 0.3,
) -> Tuple[Dict, NucleusModelSelector]:
    """
    Evaluate hybrid S/N + ML filter detection.

    Strategy:
    1. Use lunaNMR S/N detection to find candidates (high recall)
    2. Extract ML features at each S/N candidate position
    3. Apply ML classifier to filter out false positives (nucleus-specific)
    4. Return filtered peaks

    This combines S/N's high recall with ML's precision improvement.

    Args:
        model_selector: NucleusModelSelector with nucleus-specific models
        pairs: List of (spectrum_path, peaklist_path) tuples
        sn_threshold: S/N threshold for initial detection
        expected_peak_count: Max peaks from S/N
        ml_threshold: ML confidence threshold for filtering (lower = keep more)

    Returns:
        Tuple of (metrics dict, model selector)
    """

    total_tp = 0
    total_fp = 0
    total_fn = 0
    spectra_evaluated = 0
    total_time = 0.0
    total_sn_candidates = 0
    total_after_ml = 0

    # Per-nucleus tracking
    per_nucleus = {
        '15N': {'tp': 0, 'fp': 0, 'fn': 0, 'spectra': 0, 'sn_candidates': 0, 'after_ml': 0},
        '13C': {'tp': 0, 'fp': 0, 'fn': 0, 'spectra': 0, 'sn_candidates': 0, 'after_ml': 0},
    }

    for i, (spectrum_path, txt_path) in enumerate(pairs):
        # Load reference peaks from .txt
        reference_peaks = load_peaklist_txt(txt_path)
        if len(reference_peaks) < 3:
            continue

        try:
            # Step 1: Use lunaNMR S/N detection (high recall)
            t_start = time.perf_counter()

            integrator = EnhancedVoigtIntegrator()
            integrator.load_nmr_file(str(spectrum_path))
            integrator.sn_threshold = sn_threshold
            integrator.expected_peak_count = expected_peak_count

            sn_detected = integrator.detect_peaks_sn_native()

            if not sn_detected:
                continue

            total_sn_candidates += len(sn_detected)

            # Get spectrum data for feature extraction
            spectrum = integrator.nmr_data
            ppm_x = integrator.ppm_x_axis
            ppm_y = integrator.ppm_y_axis
            noise_level = integrator.noise_level

            # Detect nucleus type
            nucleus_type = _nucleus_detector._detect_nucleus_type(ppm_y)

            # Step 2: Extract ML features at S/N candidate positions
            ppm_per_point_f1 = abs(ppm_y[1] - ppm_y[0]) if len(ppm_y) > 1 else 0.1
            ppm_per_point_f2 = abs(ppm_x[1] - ppm_x[0]) if len(ppm_x) > 1 else 0.01
            extractor = PeakFeatureExtractor(
                noise_level=noise_level,
                nucleus_type=nucleus_type,
                ppm_per_point_f1=ppm_per_point_f1,
                ppm_per_point_f2=ppm_per_point_f2,
            )

            # Convert S/N detections to index-based candidates for feature extraction
            candidates = []
            for peak in sn_detected:
                # lunaNMR uses ppm_x (F2) and ppm_y (F1)
                pos_f2 = peak['ppm_x']  # 1H
                pos_f1 = peak['ppm_y']  # 15N/13C
                y_idx = ppm_to_index(ppm_y, pos_f1)
                x_idx = ppm_to_index(ppm_x, pos_f2)
                candidates.append({
                    'y_idx': y_idx,
                    'x_idx': x_idx,
                    'ppm_x': pos_f2,
                    'ppm_y': pos_f1,
                })

            # Extract features batch
            features = extractor.extract_features_batch(spectrum, candidates)

            # Step 3: ML classification to filter (using nucleus-specific model)
            classifier = model_selector.get_model(nucleus_type)
            is_peak, confidence = classifier.classify(features, ml_threshold)

            # Filter to keep only ML-approved peaks
            filtered = [c for c, passed in zip(candidates, is_peak) if passed]

            t_elapsed = time.perf_counter() - t_start
            total_time += t_elapsed
            total_after_ml += len(filtered)

        except Exception as e:
            logger.warning(f"Error with hybrid detection for {spectrum_path.name}: {e}")
            continue

        # Match against reference (using filtered candidates in lunaNMR format)
        tp, fp, fn = match_peaks_lunaNMR(filtered, reference_peaks)

        total_tp += tp
        total_fp += fp
        total_fn += fn
        spectra_evaluated += 1

        # Update per-nucleus tracking
        if nucleus_type in per_nucleus:
            per_nucleus[nucleus_type]['tp'] += tp
            per_nucleus[nucleus_type]['fp'] += fp
            per_nucleus[nucleus_type]['fn'] += fn
            per_nucleus[nucleus_type]['spectra'] += 1
            per_nucleus[nucleus_type]['sn_candidates'] += len(sn_detected)
            per_nucleus[nucleus_type]['after_ml'] += len(filtered)

        logger.info(f"[{i+1}/{len(pairs)}] {spectrum_path.name} ({nucleus_type}, Hybrid S/N={sn_threshold}+ML={ml_threshold}): "
                   f"S/N={len(sn_detected)}, after ML={len(filtered)}, ref={len(reference_peaks)}, "
                   f"TP={tp}, FP={fp}, FN={fn}")

    precision = total_tp / (total_tp + total_fp) if (total_tp + total_fp) > 0 else 0
    recall = total_tp / (total_tp + total_fn) if (total_tp + total_fn) > 0 else 0
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0

    avg_time = total_time / spectra_evaluated if spectra_evaluated > 0 else 0
    filter_rate = 1 - (total_after_ml / total_sn_candidates) if total_sn_candidates > 0 else 0

    results = {
        'method': f'Hybrid S/N+ML (S/N={sn_threshold}, ML={ml_threshold})',
        'sn_threshold': sn_threshold,
        'ml_threshold': ml_threshold,
        'expected_peak_count': expected_peak_count,
        'spectra_evaluated': spectra_evaluated,
        'total_tp': total_tp,
        'total_fp': total_fp,
        'total_fn': total_fn,
        'precision': precision,
        'recall': recall,
        'f1': f1,
        'total_time': total_time,
        'avg_time_per_spectrum': avg_time,
        'total_sn_candidates': total_sn_candidates,
        'total_after_ml': total_after_ml,
        'filter_rate': filter_rate,
        'per_nucleus': per_nucleus,
    }
    return results, model_selector


def print_per_nucleus_metrics(hybrid_results: Dict):
    """Print separate metrics for 15N and 13C spectra."""
    per_nucleus = hybrid_results.get('per_nucleus', {})

    if not per_nucleus:
        print("\nNo per-nucleus data available.")
        return

    print("\n" + "=" * 85)
    print("PER-NUCLEUS METRICS BREAKDOWN")
    print("=" * 85)

    print(f"{'Metric':<25} {'15N Spectra':<28} {'13C Spectra':<28}")
    print("-" * 85)

    # Get data for each nucleus
    data_15N = per_nucleus.get('15N', {})
    data_13C = per_nucleus.get('13C', {})

    # Spectra count
    spectra_15N = data_15N.get('spectra', 0)
    spectra_13C = data_13C.get('spectra', 0)
    print(f"{'Spectra Evaluated':<25} {spectra_15N:<28} {spectra_13C:<28}")

    # TP/FP/FN
    tp_15N = data_15N.get('tp', 0)
    fp_15N = data_15N.get('fp', 0)
    fn_15N = data_15N.get('fn', 0)

    tp_13C = data_13C.get('tp', 0)
    fp_13C = data_13C.get('fp', 0)
    fn_13C = data_13C.get('fn', 0)

    print(f"{'True Positives':<25} {tp_15N:<28} {tp_13C:<28}")
    print(f"{'False Positives':<25} {fp_15N:<28} {fp_13C:<28}")
    print(f"{'False Negatives':<25} {fn_15N:<28} {fn_13C:<28}")
    print("-" * 85)

    # Precision/Recall/F1 for 15N
    prec_15N = tp_15N / (tp_15N + fp_15N) if (tp_15N + fp_15N) > 0 else 0
    rec_15N = tp_15N / (tp_15N + fn_15N) if (tp_15N + fn_15N) > 0 else 0
    f1_15N = 2 * prec_15N * rec_15N / (prec_15N + rec_15N) if (prec_15N + rec_15N) > 0 else 0

    # Precision/Recall/F1 for 13C
    prec_13C = tp_13C / (tp_13C + fp_13C) if (tp_13C + fp_13C) > 0 else 0
    rec_13C = tp_13C / (tp_13C + fn_13C) if (tp_13C + fn_13C) > 0 else 0
    f1_13C = 2 * prec_13C * rec_13C / (prec_13C + rec_13C) if (prec_13C + rec_13C) > 0 else 0

    print(f"{'Precision':<25} {prec_15N:.4f}{'':<24} {prec_13C:.4f}")
    print(f"{'Recall':<25} {rec_15N:.4f}{'':<24} {rec_13C:.4f}")
    print(f"{'F1 Score':<25} {f1_15N:.4f}{'':<24} {f1_13C:.4f}")
    print("-" * 85)

    # S/N candidates and filter rate
    sn_15N = data_15N.get('sn_candidates', 0)
    sn_13C = data_13C.get('sn_candidates', 0)
    after_15N = data_15N.get('after_ml', 0)
    after_13C = data_13C.get('after_ml', 0)

    filter_15N = 1 - (after_15N / sn_15N) if sn_15N > 0 else 0
    filter_13C = 1 - (after_13C / sn_13C) if sn_13C > 0 else 0

    print(f"{'S/N Candidates':<25} {sn_15N:<28} {sn_13C:<28}")
    print(f"{'After ML Filter':<25} {after_15N:<28} {after_13C:<28}")
    print(f"{'Filter Rate':<25} {filter_15N*100:.1f}% rejected{'':<16} {filter_13C*100:.1f}% rejected")
    print("=" * 85)

    # Analysis
    print("\nPER-NUCLEUS ANALYSIS:")

    if spectra_15N > 0:
        print(f"\n  15N Spectra ({spectra_15N} total):")
        print(f"    F1 Score: {f1_15N:.2%} (Precision: {prec_15N:.2%}, Recall: {rec_15N:.2%})")
        if sn_15N > 0:
            print(f"    ML filter rejected {filter_15N*100:.1f}% of S/N candidates")
    else:
        print("\n  15N Spectra: No data")

    if spectra_13C > 0:
        print(f"\n  13C Spectra ({spectra_13C} total):")
        print(f"    F1 Score: {f1_13C:.2%} (Precision: {prec_13C:.2%}, Recall: {rec_13C:.2%})")
        if sn_13C > 0:
            print(f"    ML filter rejected {filter_13C*100:.1f}% of S/N candidates")
    else:
        print("\n  13C Spectra: No data")

    # Compare if both have data
    if spectra_15N > 0 and spectra_13C > 0:
        print("\n  COMPARISON:")
        if f1_15N > f1_13C:
            diff = f1_15N - f1_13C
            print(f"    15N performs {diff*100:.1f}% better (F1: {f1_15N:.2%} vs {f1_13C:.2%})")
        elif f1_13C > f1_15N:
            diff = f1_13C - f1_15N
            print(f"    13C performs {diff*100:.1f}% better (F1: {f1_13C:.2%} vs {f1_15N:.2%})")
        else:
            print(f"    Both nuclei perform equally (F1: {f1_15N:.2%})")


def evaluate_sn_detector_txt(
    pairs: List[Tuple[Path, Path]],
    sn_threshold: float = 10.0,
    expected_peak_count: int = 250,
) -> Dict:
    """
    Evaluate lunaNMR S/N detection using .txt peak lists as ground truth.

    Args:
        pairs: List of (spectrum_path, peaklist_path) tuples
        sn_threshold: S/N threshold multiplier
        expected_peak_count: Max peaks to return

    Returns metrics dict.
    """
    total_tp = 0
    total_fp = 0
    total_fn = 0
    spectra_evaluated = 0
    total_time = 0.0

    for i, (spectrum_path, txt_path) in enumerate(pairs):
        # Load reference peaks from .txt
        reference_peaks = load_peaklist_txt(txt_path)
        if len(reference_peaks) < 3:
            continue

        try:
            # Use lunaNMR's actual detection
            integrator = EnhancedVoigtIntegrator()
            integrator.load_nmr_file(str(spectrum_path))

            # Configure parameters
            integrator.sn_threshold = sn_threshold
            integrator.expected_peak_count = expected_peak_count

            # Time only detection
            t_start = time.perf_counter()
            detected = integrator.detect_peaks_sn_native()
            t_elapsed = time.perf_counter() - t_start
            total_time += t_elapsed

            if not detected:
                continue

        except Exception as e:
            logger.warning(f"Error with lunaNMR detection for {spectrum_path.name}: {e}")
            continue

        # Match against reference
        tp, fp, fn = match_peaks_lunaNMR(detected, reference_peaks)

        total_tp += tp
        total_fp += fp
        total_fn += fn
        spectra_evaluated += 1

        logger.info(f"[{i+1}/{len(pairs)}] {spectrum_path.name} (lunaNMR S/N={sn_threshold}): "
                   f"detected={len(detected)}, ref={len(reference_peaks)}, "
                   f"TP={tp}, FP={fp}, FN={fn}, time={t_elapsed*1000:.1f}ms")

    precision = total_tp / (total_tp + total_fp) if (total_tp + total_fp) > 0 else 0
    recall = total_tp / (total_tp + total_fn) if (total_tp + total_fn) > 0 else 0
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0

    avg_time = total_time / spectra_evaluated if spectra_evaluated > 0 else 0

    return {
        'method': f'lunaNMR S/N (threshold={sn_threshold})',
        'sn_threshold': sn_threshold,
        'expected_peak_count': expected_peak_count,
        'spectra_evaluated': spectra_evaluated,
        'total_tp': total_tp,
        'total_fp': total_fp,
        'total_fn': total_fn,
        'precision': precision,
        'recall': recall,
        'f1': f1,
        'total_time': total_time,
        'avg_time_per_spectrum': avg_time,
    }


def main():
    parser = argparse.ArgumentParser(description='Evaluate peak classifier against expert peak lists')
    parser.add_argument('--model', type=Path, required=False, default=None,
                       help='Path to trained model (.joblib) - used for all spectra')
    parser.add_argument('--model-15N', type=Path, required=False, default=None,
                       help='Path to 15N-specific model (.joblib)')
    parser.add_argument('--model-13C', type=Path, required=False, default=None,
                       help='Path to 13C-specific model (.joblib)')
    parser.add_argument('--spectra-dirs', type=Path, nargs='+', required=True,
                       help='Directories to search for spectrum files with matching .txt peak lists')
    parser.add_argument('--threshold', type=float, default=0.4,
                       help='ML classification threshold (default: 0.4)')
    parser.add_argument('--ml-prefilter', type=float, default=3.0,
                       help='Minimum S/N pre-filter for ML candidates (default: 3.0)')
    parser.add_argument('--sn-threshold', type=float, default=10.0,
                       help='S/N threshold for lunaNMR detection (default: 10.0)')
    parser.add_argument('--expected-peaks', type=int, default=250,
                       help='Expected peak count for lunaNMR detection (default: 250)')
    parser.add_argument('--test-only', action='store_true',
                       help='Only evaluate on test spectra (from model v2.0.0 train/test split)')
    parser.add_argument('--hybrid', action='store_true',
                       help='Evaluate hybrid mode: S/N detection + ML filter (instead of standalone ML)')
    parser.add_argument('--hybrid-ml-threshold', type=float, default=0.3,
                       help='ML threshold for hybrid filtering (default: 0.3, lower = keep more)')
    parser.add_argument('--per-nucleus', action='store_true',
                       help='Report metrics separately for 15N and 13C spectra')

    args = parser.parse_args()

    # Validate model arguments
    if not args.model and not args.model_15N and not args.model_13C:
        parser.error("Must provide --model or at least one of --model-15N / --model-13C")

    # Determine if using nucleus-specific models
    use_nucleus_specific = args.model_15N or args.model_13C
    if use_nucleus_specific:
        logger.info("Using NUCLEUS-SPECIFIC models:")
        if args.model_15N:
            logger.info(f"  15N model: {args.model_15N}")
        if args.model_13C:
            logger.info(f"  13C model: {args.model_13C}")
        if args.model:
            logger.info(f"  Fallback model: {args.model}")
    else:
        logger.info(f"Using single model for all spectra: {args.model}")

    # Find spectrum/peak-list pairs
    logger.info(f"Searching for spectrum files with matching .txt peak lists...")
    pairs = find_spectrum_peaklist_pairs(args.spectra_dirs)
    logger.info(f"Found {len(pairs)} spectrum/peak-list pairs")

    if len(pairs) == 0:
        logger.error("No spectrum/peak-list pairs found!")
        sys.exit(1)

    # Filter to test-only spectra if requested
    if args.test_only:
        # Load model(s) to get test spectra list - use any available model
        model_to_load = args.model or args.model_15N or args.model_13C
        try:
            model_data = joblib.load(model_to_load)
            if isinstance(model_data, dict) and 'test_spectra' in model_data:
                test_spectra = set(model_data['test_spectra'])
                if test_spectra:
                    original_count = len(pairs)
                    # Match by name (with ext) or stem (without ext)
                    pairs = [(sp, tp) for sp, tp in pairs
                             if sp.name in test_spectra or sp.stem in test_spectra]
                    logger.info(f"Filtered to {len(pairs)} test-only spectra (from {original_count})")
                    if len(pairs) == 0:
                        logger.error("No test spectra found in the provided directories!")
                        logger.error(f"Test spectra expected: {list(test_spectra)[:5]}...")
                        sys.exit(1)
                else:
                    logger.warning("Model has no test_spectra list - evaluating on all spectra")
            else:
                logger.warning("Model doesn't contain test_spectra (older format) - evaluating on all spectra")
        except Exception as e:
            logger.warning(f"Could not read test spectra from model: {e}")
            logger.warning("Evaluating on all spectra")

    # Show found pairs
    for spectrum_path, txt_path in pairs[:5]:
        logger.info(f"  {spectrum_path.name} -> {txt_path.name}")
    if len(pairs) > 5:
        logger.info(f"  ... and {len(pairs) - 5} more")

    # Evaluate lunaNMR S/N baseline (always needed)
    logger.info("\n" + "=" * 50)
    logger.info(f"Evaluating lunaNMR S/N Detection (S/N={args.sn_threshold}, peaks={args.expected_peaks})...")
    logger.info("=" * 50)
    sn_results = evaluate_sn_detector_txt(
        pairs,
        sn_threshold=args.sn_threshold,
        expected_peak_count=args.expected_peaks,
    )

    if args.hybrid:
        # Hybrid mode: S/N detection + ML filter
        logger.info("\n" + "=" * 50)
        logger.info(f"Evaluating HYBRID Detector (S/N={args.sn_threshold} + ML filter threshold={args.hybrid_ml_threshold})...")
        logger.info("=" * 50)

        # Create model selector for nucleus-specific models
        model_selector = NucleusModelSelector(
            model_path=args.model,
            model_15N_path=args.model_15N,
            model_13C_path=args.model_13C,
        )

        hybrid_results, model_selector = evaluate_hybrid_detector_txt(
            model_selector, pairs,
            sn_threshold=args.sn_threshold,
            expected_peak_count=args.expected_peaks,
            ml_threshold=args.hybrid_ml_threshold,
        )

        # Print hybrid comparison
        print_hybrid_comparison(hybrid_results, sn_results)

        # Print per-nucleus metrics if requested
        if args.per_nucleus:
            print_per_nucleus_metrics(hybrid_results)

        # Feature importance (use any model from selector)
        try:
            classifier = model_selector.get_any_model()
            importance = classifier.get_feature_importance()
            print("\nFEATURE IMPORTANCE:")
            for name, imp in sorted(importance.items(), key=lambda x: -x[1]):
                bar = "█" * int(imp * 50)
                print(f"  {name:<20} {imp:.4f} {bar}")
        except Exception as e:
            logger.warning(f"Could not get feature importance: {e}")

    else:
        # Standard mode: standalone ML detector (use single model)
        if not args.model:
            logger.error("--model is required for non-hybrid mode")
            sys.exit(1)

        logger.info("\n" + "=" * 50)
        logger.info(f"Evaluating ML Peak Detector (threshold={args.threshold}, prefilter={args.ml_prefilter})...")
        logger.info("=" * 50)
        ml_results, classifier = evaluate_ml_detector_txt(
            args.model, pairs,
            threshold=args.threshold,
            min_snr_prefilter=args.ml_prefilter,
        )

        # Print comparison
        print_comparison(ml_results, sn_results)

        # Feature importance
        try:
            importance = classifier.get_feature_importance()
            print("\nFEATURE IMPORTANCE:")
            for name, imp in sorted(importance.items(), key=lambda x: -x[1]):
                bar = "█" * int(imp * 50)
                print(f"  {name:<20} {imp:.4f} {bar}")
        except Exception as e:
            logger.warning(f"Could not get feature importance: {e}")


if __name__ == '__main__':
    main()

#!/usr/bin/env python3
# ABOUTME: Script to train peak classifier from collected training data and NMR spectra.
# ABOUTME: Supports Random Forest and XGBoost with spectrum-level train/test splitting.

"""
Train Peak Classifier

This script:
1. Loads training JSON (peak positions from collected data)
2. Splits spectra into train/test sets (spectrum-level, not sample-level)
3. Loads spectrum files and extracts features
4. Extracts features at peak positions (positive samples)
5. Extracts features from hard negatives or noise regions (negative samples)
6. Trains Random Forest or XGBoost classifier
7. Saves model and test spectrum list

Usage:
    python train_peak_classifier.py --training-dir ml_training_data/adaptativev2 \
                                   --spectra-dirs /path/to/spectra \
                                   --output lunaNMR/ml/pretrained/peak_classifier.joblib \
                                   --model-type xgb \
                                   --hard-negatives
"""

import argparse
import json
import logging
import sys
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
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import StandardScaler

# XGBoost and LightGBM imported lazily when needed to avoid library load errors
XGBClassifier = None
LGBMClassifier = None

def get_xgboost_classifier():
    """Lazily import and return XGBClassifier."""
    global XGBClassifier
    if XGBClassifier is None:
        try:
            from xgboost import XGBClassifier as XGB
            XGBClassifier = XGB
        except Exception as e:
            raise ImportError(f"XGBoost not available: {e}. Install with: pip install xgboost")
    return XGBClassifier

def get_lightgbm_classifier():
    """Lazily import and return LGBMClassifier."""
    global LGBMClassifier
    if LGBMClassifier is None:
        try:
            from lightgbm import LGBMClassifier as LGBM
            LGBMClassifier = LGBM
        except Exception as e:
            raise ImportError(f"LightGBM not available: {e}. Install with: pip install lightgbm")
    return LGBMClassifier

from lunaNMR.ml.peak_feature_extractor import PeakFeatureExtractor, PeakFeatures
from lunaNMR.ml.peak_classifier import PeakClassifier
import joblib

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def load_training_data(training_dir: Path) -> Dict:
    """Load training data from JSON file."""
    training_file = training_dir / "training_data.json"
    if not training_file.exists():
        raise FileNotFoundError(f"Training data not found: {training_file}")

    with open(training_file, 'r') as f:
        data = json.load(f)

    logger.info(f"Loaded {data.get('n_samples', len(data.get('samples', [])))} training samples")
    return data


def find_spectrum_file(spectrum_name: str, spectra_dirs: List[Path]) -> Optional[Path]:
    """
    Find spectrum file by name in search directories.

    Searches for common NMR file formats.
    """
    # NMR file extensions (in priority order)
    nmr_extensions = ['.ft', '.ft2', '.ft3', '.pipe', '.ucsf']

    # Also try without extension
    base_name = Path(spectrum_name).stem

    for spectra_dir in spectra_dirs:
        # First pass: try exact filename matches
        for ext in [''] + nmr_extensions:
            # Try exact name
            path = spectra_dir / (spectrum_name + ext)
            if path.exists() and path.suffix.lower() in nmr_extensions + ['']:
                return path

            # Try base name with extension
            path = spectra_dir / (base_name + ext)
            if path.exists() and path.suffix.lower() in nmr_extensions + ['']:
                return path

        # Second pass: recursive search with NMR extensions only
        for ext in nmr_extensions:
            matches = list(spectra_dir.rglob(f"{base_name}{ext}"))
            if matches:
                return matches[0]
            # Try with wildcard prefix (for directories)
            matches = list(spectra_dir.rglob(f"*/{base_name}{ext}"))
            if matches:
                return matches[0]

    return None


def load_spectrum(spectrum_path: Path) -> Tuple[np.ndarray, np.ndarray, np.ndarray, Dict]:
    """
    Load NMR spectrum file.

    Returns:
        Tuple of (spectrum_data, ppm_x_axis, ppm_y_axis, metadata)
    """
    # Determine format and load
    suffix = spectrum_path.suffix.lower()

    if suffix in ['.ft', '.ft2', '.ft3', '.pipe', '']:
        # NMRPipe format
        try:
            nmr_dict, nmr_data = ng.pipe.read(str(spectrum_path))
            uc_x = ng.pipe.make_uc(nmr_dict, nmr_data, dim=1)
            uc_y = ng.pipe.make_uc(nmr_dict, nmr_data, dim=0)
        except:
            # Try Bruker format
            nmr_dict, nmr_data = ng.bruker.read(str(spectrum_path.parent))
            udic = ng.bruker.guess_udic(nmr_dict, nmr_data, strip_fake=True)
            uc_x = ng.fileiobase.uc_from_udic(udic, dim=1)
            uc_y = ng.fileiobase.uc_from_udic(udic, dim=0)
    elif suffix == '.ucsf':
        # SPARKY format
        nmr_dict, nmr_data = ng.sparky.read(str(spectrum_path))
        uc_x = ng.sparky.make_uc(nmr_dict, nmr_data, dim=1)
        uc_y = ng.sparky.make_uc(nmr_dict, nmr_data, dim=0)
    else:
        raise ValueError(f"Unknown format: {suffix}")

    # Generate PPM axes
    ppm_x = uc_x.ppm_scale()
    ppm_y = uc_y.ppm_scale()

    return nmr_data, ppm_x, ppm_y, nmr_dict


def ppm_to_index(ppm_axis: np.ndarray, ppm_value: float) -> int:
    """Convert PPM value to nearest index."""
    return int(np.argmin(np.abs(ppm_axis - ppm_value)))


def get_valid_noise_region(
    ppm_x: np.ndarray,
    ppm_y: np.ndarray,
    nucleus_type: str,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Get valid noise sampling region based on nucleus type.

    Uses same constraints as core_integrator.py:3712-3720.

    Returns:
        Tuple of (x_valid_mask, y_valid_mask)
    """
    # Nucleus-adaptive 1H constraint
    if nucleus_type == '15N':
        # 15N-HSQC: 1H >= 5.8 ppm (exclude water/aliphatic, keep amide)
        x_valid = ppm_x >= 5.8
    elif nucleus_type == '13C':
        # 13C-HSQC: 1H <= 2.75 ppm (exclude downfield, keep aliphatic)
        x_valid = ppm_x <= 2.75
    else:
        x_valid = np.ones(len(ppm_x), dtype=bool)

    # Exclude Y-dimension edges (0.5 ppm margin)
    Y_EDGE_MARGIN = 0.5
    y_min_edge = ppm_y.min() + Y_EDGE_MARGIN
    y_max_edge = ppm_y.max() - Y_EDGE_MARGIN
    y_valid = (ppm_y >= y_min_edge) & (ppm_y <= y_max_edge)

    return x_valid, y_valid


def sample_noise_locations(
    spectrum: np.ndarray,
    ppm_x: np.ndarray,
    ppm_y: np.ndarray,
    peak_positions: List[Tuple[float, float]],
    nucleus_type: str,
    n_samples: int,
    min_distance_ppm: float = 0.5,
) -> List[Tuple[int, int]]:
    """
    Sample random noise locations from valid regions.

    Args:
        spectrum: 2D spectrum data
        ppm_x: X-axis PPM values
        ppm_y: Y-axis PPM values
        peak_positions: List of (pos_f1, pos_f2) peak positions
        nucleus_type: '15N' or '13C'
        n_samples: Number of noise samples to generate
        min_distance_ppm: Minimum distance from peaks (in F1 dimension)

    Returns:
        List of (y_idx, x_idx) noise sample locations
    """
    # Get valid region masks
    x_valid, y_valid = get_valid_noise_region(ppm_x, ppm_y, nucleus_type)

    # Get indices of valid region
    valid_x_indices = np.where(x_valid)[0]
    valid_y_indices = np.where(y_valid)[0]

    if len(valid_x_indices) == 0 or len(valid_y_indices) == 0:
        logger.warning(f"No valid noise region for nucleus type {nucleus_type}")
        return []

    # Convert peak positions to indices
    peak_indices = []
    for pos_f1, pos_f2 in peak_positions:
        y_idx = ppm_to_index(ppm_y, pos_f1)
        x_idx = ppm_to_index(ppm_x, pos_f2)
        peak_indices.append((y_idx, x_idx))

    # Sample random locations
    noise_locations = []
    max_attempts = n_samples * 10
    attempts = 0

    while len(noise_locations) < n_samples and attempts < max_attempts:
        # Random indices from valid regions
        y_idx = np.random.choice(valid_y_indices)
        x_idx = np.random.choice(valid_x_indices)

        # Check distance from all peaks
        too_close = False
        for py, px in peak_indices:
            # Use F1 (Y) dimension for distance (larger PPM range)
            dy = abs(ppm_y[y_idx] - ppm_y[py])
            if dy < min_distance_ppm:
                too_close = True
                break

        if not too_close:
            noise_locations.append((y_idx, x_idx))

        attempts += 1

    return noise_locations


def find_local_maxima(
    spectrum: np.ndarray,
    min_snr: float = 3.0,
    noise_level: float = 1.0,
) -> List[Tuple[int, int]]:
    """
    Find all local maxima in spectrum above S/N threshold.

    Returns list of (y_idx, x_idx) tuples.
    """
    # Local maximum filter (3x3 neighborhood)
    local_max = maximum_filter(spectrum, size=3)
    is_max = (spectrum == local_max) & (spectrum > 0)

    # Apply S/N filter
    if min_snr > 0:
        is_max = is_max & (spectrum > min_snr * noise_level)

    # Get coordinates
    y_indices, x_indices = np.where(is_max)

    return list(zip(y_indices.tolist(), x_indices.tolist()))


def find_hard_negatives(
    spectrum: np.ndarray,
    ppm_x: np.ndarray,
    ppm_y: np.ndarray,
    peak_positions: List[Tuple[float, float]],
    nucleus_type: str,
    noise_level: float,
    min_snr_prefilter: float = 3.0,
    tolerance_f1: float = 0.3,
    tolerance_f2: float = 0.03,
) -> List[Tuple[int, int]]:
    """
    Find hard negatives: local maxima that are NOT real peaks.

    These are the challenging cases that look like peaks but aren't.

    Args:
        spectrum: 2D spectrum data
        ppm_x: X-axis PPM values
        ppm_y: Y-axis PPM values
        peak_positions: List of (pos_f1, pos_f2) real peak positions
        nucleus_type: '15N' or '13C'
        noise_level: Noise level for S/N filtering
        min_snr_prefilter: Minimum S/N for candidates
        tolerance_f1: Matching tolerance in F1 (ppm)
        tolerance_f2: Matching tolerance in F2 (ppm)

    Returns:
        List of (y_idx, x_idx) hard negative locations
    """
    # Get valid detection region
    x_valid, y_valid = get_valid_noise_region(ppm_x, ppm_y, nucleus_type)

    # Find all local maxima
    all_maxima = find_local_maxima(spectrum, min_snr=min_snr_prefilter, noise_level=noise_level)

    # Filter to valid region and check against real peaks
    hard_negatives = []

    for y_idx, x_idx in all_maxima:
        # Check if in valid region
        if not x_valid[x_idx] or not y_valid[y_idx]:
            continue

        # Get PPM coordinates of this maximum
        max_f1 = ppm_y[y_idx]
        max_f2 = ppm_x[x_idx]

        # Check if this matches any real peak
        is_real_peak = False
        for peak_f1, peak_f2 in peak_positions:
            d_f1 = abs(max_f1 - peak_f1)
            d_f2 = abs(max_f2 - peak_f2)

            if d_f1 <= tolerance_f1 and d_f2 <= tolerance_f2:
                is_real_peak = True
                break

        if not is_real_peak:
            hard_negatives.append((y_idx, x_idx))

    return hard_negatives


def extract_training_features(
    sample: Dict,
    spectrum: np.ndarray,
    ppm_x: np.ndarray,
    ppm_y: np.ndarray,
    min_r2: float = 0.85,
    use_hard_negatives: bool = True,
    hard_neg_prefilter: float = 3.0,
) -> Tuple[List[np.ndarray], List[int], int, int]:
    """
    Extract features from one training sample.

    Args:
        sample: Training sample dict with peaks
        spectrum: 2D spectrum data
        ppm_x, ppm_y: PPM axes
        min_r2: Minimum R² for peaks to include
        use_hard_negatives: If True, use local maxima as negatives (hard negatives)
                           If False, use random noise locations
        hard_neg_prefilter: S/N prefilter for hard negative candidates

    Returns:
        Tuple of (features_list, labels_list, n_peaks, n_noise)
    """
    features_list = []
    labels_list = []

    noise_level = sample.get('noise_level', 1.0)
    nucleus_type = sample.get('nucleus_type', '15N')

    # Compute ppm per point from axes
    ppm_per_point_f1 = abs(ppm_y[1] - ppm_y[0]) if len(ppm_y) > 1 else 0.1
    ppm_per_point_f2 = abs(ppm_x[1] - ppm_x[0]) if len(ppm_x) > 1 else 0.01

    extractor = PeakFeatureExtractor(
        noise_level=noise_level,
        nucleus_type=nucleus_type,
        ppm_per_point_f1=ppm_per_point_f1,
        ppm_per_point_f2=ppm_per_point_f2,
    )

    # Get peak positions
    peaks = sample.get('peaks', [])
    peak_positions = []

    # Extract positive samples (peaks with good R²)
    for peak in peaks:
        r_squared = peak.get('r_squared', 0)
        if r_squared < min_r2:
            continue

        pos_f1 = peak.get('pos_f1', 0)
        pos_f2 = peak.get('pos_f2', 0)

        if pos_f1 <= 0 or pos_f2 <= 0:
            continue

        peak_positions.append((pos_f1, pos_f2))

        # Convert to indices
        y_idx = ppm_to_index(ppm_y, pos_f1)
        x_idx = ppm_to_index(ppm_x, pos_f2)

        # Extract features
        features = extractor.extract_features(spectrum, y_idx, x_idx)
        features_list.append(features.to_array())
        labels_list.append(1)  # Peak

    n_peaks = len(labels_list)

    # Get negative samples
    if use_hard_negatives:
        # Hard negatives: local maxima that aren't real peaks
        noise_locations = find_hard_negatives(
            spectrum, ppm_x, ppm_y,
            peak_positions, nucleus_type, noise_level,
            min_snr_prefilter=hard_neg_prefilter,
        )
    else:
        # Random noise: sample from baseline regions
        noise_locations = sample_noise_locations(
            spectrum, ppm_x, ppm_y,
            peak_positions, nucleus_type,
            n_samples=n_peaks,
        )

    # Extract negative samples (noise/hard negatives)
    for y_idx, x_idx in noise_locations:
        features = extractor.extract_features(spectrum, y_idx, x_idx)
        features_list.append(features.to_array())
        labels_list.append(0)  # Noise

    n_noise = len(noise_locations)

    return features_list, labels_list, n_peaks, n_noise


def main():
    parser = argparse.ArgumentParser(description='Train peak classifier from collected data')
    parser.add_argument('--training-dir', type=Path, required=True,
                       help='Directory containing training_data.json')
    parser.add_argument('--spectra-dirs', type=Path, nargs='+', required=True,
                       help='Directories to search for spectrum files')
    parser.add_argument('--output', type=Path, required=True,
                       help='Output path for trained model (.joblib)')
    parser.add_argument('--min-r2', type=float, default=0.85,
                       help='Minimum R² for peaks to include (default: 0.85)')
    parser.add_argument('--test-size', type=float, default=0.15,
                       help='Fraction of SPECTRA for testing (default: 0.15)')
    parser.add_argument('--val-size', type=float, default=0.15,
                       help='Fraction of SPECTRA for validation (default: 0.15)')
    parser.add_argument('--cv', type=int, default=5,
                       help='Cross-validation folds (default: 5)')
    parser.add_argument('--dry-run', action='store_true',
                       help='Show what would be done without training')
    # Hard negative sampling options
    parser.add_argument('--hard-negatives', action='store_true', default=True,
                       help='Use hard negative sampling (local maxima that are not peaks)')
    parser.add_argument('--random-negatives', action='store_true',
                       help='Use random noise sampling instead of hard negatives')
    parser.add_argument('--hard-neg-prefilter', type=float, default=3.0,
                       help='S/N prefilter for hard negative candidates (default: 3.0)')
    # Model type
    parser.add_argument('--model-type', type=str, default='rf', choices=['rf', 'xgb', 'lgb'],
                       help='Model type: rf (Random Forest), xgb (XGBoost), or lgb (LightGBM)')
    # Nucleus-specific training
    parser.add_argument('--nucleus-type', type=str, default=None, choices=['15N', '13C'],
                       help='Train nucleus-specific model (15N or 13C). If not set, train on all data.')

    args = parser.parse_args()

    # Check XGBoost/LightGBM availability (lazy import)
    if args.model_type == 'xgb':
        try:
            get_xgboost_classifier()
        except ImportError as e:
            logger.error(str(e))
            sys.exit(1)
    elif args.model_type == 'lgb':
        try:
            get_lightgbm_classifier()
        except ImportError as e:
            logger.error(str(e))
            sys.exit(1)

    # Determine negative sampling strategy
    use_hard_negatives = args.hard_negatives and not args.random_negatives

    # Load training data
    logger.info(f"Loading training data from {args.training_dir}")
    logger.info(f"Model type: {args.model_type.upper()}")
    logger.info(f"Negative sampling: {'HARD NEGATIVES' if use_hard_negatives else 'RANDOM NOISE'}")
    if use_hard_negatives:
        logger.info(f"  Hard negative S/N prefilter: {args.hard_neg_prefilter}")
    training_data = load_training_data(args.training_dir)
    samples = training_data.get('samples', [])

    if not samples:
        logger.error("No training samples found")
        sys.exit(1)

    # First pass: find all available spectra
    logger.info("\n=== Phase 1: Finding available spectra ===")
    available_samples = []
    for sample in samples:
        spectrum_name = sample.get('spectrum_name', '')
        if not spectrum_name:
            continue
        spectrum_path = find_spectrum_file(spectrum_name, args.spectra_dirs)
        if spectrum_path is not None:
            available_samples.append((sample, spectrum_path))

    logger.info(f"Found {len(available_samples)} spectra (out of {len(samples)} samples)")

    # Filter by nucleus type if specified
    if args.nucleus_type:
        original_count = len(available_samples)
        available_samples = [
            (sample, path) for sample, path in available_samples
            if sample.get('nucleus_type', '15N') == args.nucleus_type
        ]
        logger.info(f"Filtered to {len(available_samples)} {args.nucleus_type} spectra (from {original_count})")

    if len(available_samples) < 7:
        logger.error("Need at least 7 spectra for spectrum-level splitting")
        sys.exit(1)

    # Spectrum-level train/val/test split
    logger.info("\n=== Phase 2: Spectrum-level train/val/test split ===")
    np.random.seed(42)
    indices = np.random.permutation(len(available_samples))

    n_test = max(1, int(len(available_samples) * args.test_size))
    n_val = max(1, int(len(available_samples) * args.val_size))
    n_train = len(available_samples) - n_test - n_val

    test_indices = indices[:n_test]
    val_indices = indices[n_test:n_test + n_val]
    train_indices = indices[n_test + n_val:]

    train_samples = [available_samples[i] for i in train_indices]
    val_samples = [available_samples[i] for i in val_indices]
    test_samples = [available_samples[i] for i in test_indices]

    logger.info(f"  Train spectra: {len(train_samples)}")
    logger.info(f"  Validation spectra: {len(val_samples)}")
    logger.info(f"  Test spectra: {len(test_samples)}")

    # Save test spectra list
    test_spectra_file = args.output.parent / (args.output.stem + "_test_spectra.json")
    test_spectra_names = [s[0].get('spectrum_name', '') for s in test_samples]
    with open(test_spectra_file, 'w') as f:
        json.dump({'test_spectra': test_spectra_names}, f, indent=2)
    logger.info(f"  Test spectra list saved to: {test_spectra_file}")

    if args.dry_run:
        logger.info("\nDry run - would train on these spectra:")
        for sample, path in train_samples[:5]:
            logger.info(f"  {sample.get('spectrum_name', '')}")
        if len(train_samples) > 5:
            logger.info(f"  ... and {len(train_samples) - 5} more")
        return

    # Extract features from each split
    def extract_from_samples(sample_list, split_name):
        features_list = []
        labels_list = []
        total_peaks = 0
        total_neg = 0

        for i, (sample, spectrum_path) in enumerate(sample_list):
            try:
                spectrum, ppm_x, ppm_y, _ = load_spectrum(spectrum_path)
                features, labels, n_peaks, n_noise = extract_training_features(
                    sample, spectrum, ppm_x, ppm_y,
                    min_r2=args.min_r2,
                    use_hard_negatives=use_hard_negatives,
                    hard_neg_prefilter=args.hard_neg_prefilter,
                )
                features_list.extend(features)
                labels_list.extend(labels)
                total_peaks += n_peaks
                total_neg += n_noise

                neg_type = "hard neg" if use_hard_negatives else "noise"
                logger.info(f"  [{split_name}] [{i+1}/{len(sample_list)}] {sample.get('spectrum_name', '')}: "
                           f"{n_peaks} peaks, {n_noise} {neg_type}")

            except Exception as e:
                logger.error(f"  Error processing {sample.get('spectrum_name', '')}: {e}")
                continue

        return np.array(features_list), np.array(labels_list), total_peaks, total_neg

    logger.info("\n=== Phase 3: Extracting features ===")
    X_train, y_train, train_peaks, train_neg = extract_from_samples(train_samples, "TRAIN")
    X_val, y_val, val_peaks, val_neg = extract_from_samples(val_samples, "VAL")
    X_test, y_test, test_peaks, test_neg = extract_from_samples(test_samples, "TEST")

    logger.info(f"\nFeature extraction complete:")
    logger.info(f"  Train: {len(y_train)} samples ({train_peaks} peaks, {train_neg} negatives)")
    logger.info(f"  Val: {len(y_val)} samples ({val_peaks} peaks, {val_neg} negatives)")
    logger.info(f"  Test: {len(y_test)} samples ({test_peaks} peaks, {test_neg} negatives)")

    if len(y_train) < 100:
        logger.error("Insufficient training data (need at least 100 samples)")
        sys.exit(1)

    # Scale features
    scaler = StandardScaler()
    X_train_scaled = scaler.fit_transform(X_train)
    X_val_scaled = scaler.transform(X_val)
    X_test_scaled = scaler.transform(X_test)

    # Calculate class weight for imbalanced data
    n_neg = np.sum(y_train == 0)
    n_pos = np.sum(y_train == 1)
    scale_pos_weight = n_neg / n_pos if n_pos > 0 else 1.0
    logger.info(f"\nClass balance: {n_pos} peaks, {n_neg} negatives (ratio: {scale_pos_weight:.2f})")

    # Train model
    logger.info(f"\n=== Phase 4: Training {args.model_type.upper()} classifier ===")

    if args.model_type == 'rf':
        model = RandomForestClassifier(
            n_estimators=200,
            max_depth=15,
            min_samples_leaf=3,
            class_weight='balanced',
            n_jobs=-1,
            random_state=42,
        )
    elif args.model_type == 'xgb':
        XGB = get_xgboost_classifier()
        model = XGB(
            n_estimators=200,
            max_depth=8,
            learning_rate=0.1,
            scale_pos_weight=scale_pos_weight,
            n_jobs=-1,
            random_state=42,
            eval_metric='logloss',
        )
    else:  # lgb
        LGBM = get_lightgbm_classifier()
        model = LGBM(
            n_estimators=200,
            max_depth=8,
            learning_rate=0.1,
            scale_pos_weight=scale_pos_weight,
            n_jobs=-1,
            random_state=42,
            verbose=-1,  # Suppress LightGBM output
        )

    model.fit(X_train_scaled, y_train)

    # Evaluate on validation set
    y_val_pred = model.predict(X_val_scaled)
    y_val_proba = model.predict_proba(X_val_scaled)[:, 1]

    val_tp = np.sum((y_val_pred == 1) & (y_val == 1))
    val_fp = np.sum((y_val_pred == 1) & (y_val == 0))
    val_fn = np.sum((y_val_pred == 0) & (y_val == 1))
    val_precision = val_tp / (val_tp + val_fp) if (val_tp + val_fp) > 0 else 0
    val_recall = val_tp / (val_tp + val_fn) if (val_tp + val_fn) > 0 else 0
    val_f1 = 2 * val_precision * val_recall / (val_precision + val_recall) if (val_precision + val_recall) > 0 else 0

    logger.info(f"\nValidation results:")
    logger.info(f"  Precision: {val_precision:.4f}")
    logger.info(f"  Recall: {val_recall:.4f}")
    logger.info(f"  F1: {val_f1:.4f}")

    # Evaluate on test set (held-out spectra)
    y_test_proba = model.predict_proba(X_test_scaled)[:, 1]

    for threshold in [0.3, 0.4, 0.5]:
        y_test_pred = (y_test_proba >= threshold).astype(int)
        test_tp = np.sum((y_test_pred == 1) & (y_test == 1))
        test_fp = np.sum((y_test_pred == 1) & (y_test == 0))
        test_fn = np.sum((y_test_pred == 0) & (y_test == 1))
        test_precision = test_tp / (test_tp + test_fp) if (test_tp + test_fp) > 0 else 0
        test_recall = test_tp / (test_tp + test_fn) if (test_tp + test_fn) > 0 else 0
        test_f1 = 2 * test_precision * test_recall / (test_precision + test_recall) if (test_precision + test_recall) > 0 else 0

        logger.info(f"\nTest results (threshold={threshold}) - HELD-OUT SPECTRA:")
        logger.info(f"  Precision: {test_precision:.4f}")
        logger.info(f"  Recall: {test_recall:.4f}")
        logger.info(f"  F1: {test_f1:.4f}")

    # Feature importance
    feature_names = PeakFeatures.feature_names()
    if args.model_type == 'rf':
        importances = model.feature_importances_
    else:
        importances = model.feature_importances_

    logger.info("\nFeature importance:")
    importance_dict = dict(zip(feature_names, importances))
    for name, imp in sorted(importance_dict.items(), key=lambda x: -x[1]):
        logger.info(f"  {name}: {imp:.4f}")

    # Save model (scaler + model together)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    model_data = {
        'scaler': scaler,
        'model': model,
        'model_type': args.model_type,
        'feature_names': feature_names,
        'test_spectra': test_spectra_names,
        'nucleus_type': args.nucleus_type,  # None = all nuclei, '15N' or '13C' = specific
        'version': '2.1.0',  # Bumped for nucleus-specific support
    }
    joblib.dump(model_data, args.output)
    logger.info(f"\nModel saved to {args.output}")
    if args.nucleus_type:
        logger.info(f"  Nucleus-specific model: {args.nucleus_type}")


if __name__ == '__main__':
    main()

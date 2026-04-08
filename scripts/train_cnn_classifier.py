#!/usr/bin/env python3
# ABOUTME: Training script for CNN-based peak classifier.
# ABOUTME: Trains a small CNN on 2D spectrum patches to distinguish peaks from noise.

"""
Train CNN Peak Classifier

This script trains a CNN to classify 2D spectrum patches as peaks or noise.
Uses spectrum-level train/val/test split to prevent data leakage.

Usage:
    python train_cnn_classifier.py \
        --training-dir ml_training_data/adaptativev2 \
        --spectra-dirs /path/to/spectra \
        --output lunaNMR/ml/pretrained/cnn_peak_classifier.pt

Features:
- Spectrum-level splitting (entire spectra held out for validation/test)
- Hard negative mining (local maxima that aren't real peaks)
- Data augmentation (rotation, flip, noise, intensity scaling)
- Early stopping on validation loss
- Class weighting for imbalanced data
"""

import argparse
import json
import logging
import random
import sys
import time
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

try:
    import torch
    import torch.nn as nn
    import torch.optim as optim
    from torch.utils.data import Dataset, DataLoader
except ImportError:
    print("ERROR: PyTorch is required. Install with: pip install torch torchvision")
    sys.exit(1)

try:
    import nmrglue as ng
except ImportError:
    print("ERROR: nmrglue is required. Install with: pip install nmrglue")
    sys.exit(1)

from scipy.ndimage import maximum_filter

from lunaNMR.ml.cnn_peak_classifier import create_peak_cnn, CNNConfig, count_parameters

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


class PeakPatchDataset(Dataset):
    """
    PyTorch Dataset for peak patches.

    Stores pre-extracted patches with labels.
    Applies data augmentation during training.
    """

    def __init__(
        self,
        patches: np.ndarray,
        labels: np.ndarray,
        augment: bool = False,
    ):
        """
        Args:
            patches: Array of shape (n, patch_size, patch_size)
            labels: Array of shape (n,) with 0/1 labels
            augment: Whether to apply data augmentation
        """
        self.patches = patches.astype(np.float32)
        self.labels = labels.astype(np.float32)
        self.augment = augment

    def __len__(self):
        return len(self.labels)

    def __getitem__(self, idx):
        patch = self.patches[idx].copy()
        label = self.labels[idx]

        if self.augment:
            patch = self._augment(patch)

        # Normalize (z-score per patch)
        mean = patch.mean()
        std = patch.std()
        if std > 1e-10:
            patch = (patch - mean) / std
        else:
            patch = patch - mean

        # Add channel dimension: (H, W) -> (1, H, W)
        patch = patch[np.newaxis, :, :]

        # Use torch.tensor() instead of torch.from_numpy() for numpy 2.x compatibility
        return torch.tensor(patch, dtype=torch.float32), torch.tensor(label, dtype=torch.float32)

    def _augment(self, patch: np.ndarray) -> np.ndarray:
        """Apply random augmentation."""
        # Random rotation (0, 90, 180, 270 degrees)
        k = random.randint(0, 3)
        patch = np.rot90(patch, k)

        # Random horizontal flip
        if random.random() > 0.5:
            patch = np.fliplr(patch).copy()

        # Random vertical flip
        if random.random() > 0.5:
            patch = np.flipud(patch).copy()

        # Random noise injection (10% of cases)
        if random.random() > 0.9:
            noise = np.random.normal(0, 0.1, patch.shape)
            patch = patch + noise * patch.std()

        # Random intensity scaling (0.8-1.2x)
        if random.random() > 0.5:
            scale = random.uniform(0.8, 1.2)
            patch = patch * scale

        return patch


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


def find_spectrum_file(spectrum_name: str, spectra_dirs: List[Path]) -> Optional[Path]:
    """Find spectrum file by name in search directories."""
    nmr_extensions = ['.ft', '.ft2', '.ft3', '.pipe', '.ucsf']
    base_name = Path(spectrum_name).stem

    for spectra_dir in spectra_dirs:
        for ext in [''] + nmr_extensions:
            path = spectra_dir / (spectrum_name + ext)
            if path.exists():
                return path
            path = spectra_dir / (base_name + ext)
            if path.exists():
                return path

        for ext in nmr_extensions:
            matches = list(spectra_dir.rglob(f"{base_name}{ext}"))
            if matches:
                return matches[0]

    return None


def ppm_to_index(ppm_axis: np.ndarray, ppm_value: float) -> int:
    """Convert PPM value to nearest index."""
    return int(np.argmin(np.abs(ppm_axis - ppm_value)))


def extract_patch(
    spectrum: np.ndarray,
    y_idx: int,
    x_idx: int,
    patch_size: int = 21,
) -> Optional[np.ndarray]:
    """
    Extract a patch centered at (y_idx, x_idx).

    Returns None if patch would be mostly out of bounds.
    """
    half = patch_size // 2
    ny, nx = spectrum.shape

    # Check if center is too close to edge
    if y_idx < 2 or y_idx >= ny - 2 or x_idx < 2 or x_idx >= nx - 2:
        return None

    # Extract with zero-padding
    y_start = max(0, y_idx - half)
    y_end = min(ny, y_idx + half + 1)
    x_start = max(0, x_idx - half)
    x_end = min(nx, x_idx + half + 1)

    patch = np.zeros((patch_size, patch_size), dtype=np.float32)

    py_start = half - (y_idx - y_start)
    py_end = py_start + (y_end - y_start)
    px_start = half - (x_idx - x_start)
    px_end = px_start + (x_end - x_start)

    patch[py_start:py_end, px_start:px_end] = spectrum[y_start:y_end, x_start:x_end]

    return patch


def find_hard_negatives(
    spectrum: np.ndarray,
    ppm_x: np.ndarray,
    ppm_y: np.ndarray,
    peak_positions: List[Tuple[float, float]],
    tolerance_f1: float = 0.3,
    tolerance_f2: float = 0.03,
    max_negatives: int = 500,
) -> List[Tuple[int, int]]:
    """
    Find hard negative locations (local maxima that aren't real peaks).

    Args:
        spectrum: 2D spectrum array
        ppm_x, ppm_y: PPM axes
        peak_positions: List of (pos_f1, pos_f2) for real peaks
        tolerance_f1, tolerance_f2: Exclusion radius around real peaks
        max_negatives: Maximum number of negatives to return

    Returns:
        List of (y_idx, x_idx) tuples for hard negatives
    """
    ny, nx = spectrum.shape

    # Find all local maxima
    local_max = maximum_filter(spectrum, size=3)
    is_max = (spectrum == local_max) & (spectrum > 0)

    # Get noise threshold (exclude very weak candidates)
    corners = [
        spectrum[:20, :20],
        spectrum[:20, -20:],
        spectrum[-20:, :20],
        spectrum[-20:, -20:],
    ]
    noise_level = np.std(np.concatenate([c.flatten() for c in corners]))
    is_max = is_max & (spectrum > 3 * noise_level)

    y_indices, x_indices = np.where(is_max)

    # Filter out locations near real peaks
    negatives = []
    for y_idx, x_idx in zip(y_indices, x_indices):
        pos_f1 = ppm_y[y_idx]
        pos_f2 = ppm_x[x_idx]

        is_near_peak = False
        for peak_f1, peak_f2 in peak_positions:
            if abs(pos_f1 - peak_f1) < tolerance_f1 and abs(pos_f2 - peak_f2) < tolerance_f2:
                is_near_peak = True
                break

        if not is_near_peak:
            negatives.append((y_idx, x_idx))

    # Subsample if too many
    if len(negatives) > max_negatives:
        random.shuffle(negatives)
        negatives = negatives[:max_negatives]

    return negatives


def extract_training_data(
    training_data: Dict,
    spectra_dirs: List[Path],
    patch_size: int = 21,
    min_r2: float = 0.85,
    max_negatives_per_spectrum: int = 500,
) -> Tuple[Dict[str, np.ndarray], Dict[str, np.ndarray], Dict[str, np.ndarray]]:
    """
    Extract patches from training data with spectrum-level grouping.

    Returns:
        Tuple of (train_data, val_data, test_data) dicts
        Each dict has 'patches' and 'labels' arrays
    """
    samples = training_data.get('samples', [])
    logger.info(f"Processing {len(samples)} samples from training data")

    # Group by spectrum
    spectrum_patches = {}  # spectrum_name -> (patches, labels)

    for i, sample in enumerate(samples):
        spectrum_name = sample.get('spectrum_name', '')
        if not spectrum_name:
            continue

        # Find spectrum file
        spectrum_path = find_spectrum_file(spectrum_name, spectra_dirs)
        if spectrum_path is None:
            logger.debug(f"Spectrum not found: {spectrum_name}")
            continue

        try:
            spectrum, ppm_x, ppm_y = load_spectrum(spectrum_path)
        except Exception as e:
            logger.warning(f"Error loading {spectrum_name}: {e}")
            continue

        peaks = sample.get('peaks', [])
        good_peaks = [
            p for p in peaks
            if p.get('r_squared', 0) >= min_r2
            and p.get('pos_f1', 0) > 0
            and p.get('pos_f2', 0) > 0
        ]

        if len(good_peaks) < 5:
            continue

        # Extract positive patches (real peaks)
        pos_patches = []
        peak_positions = []
        for peak in good_peaks:
            pos_f1 = peak['pos_f1']
            pos_f2 = peak['pos_f2']
            peak_positions.append((pos_f1, pos_f2))

            y_idx = ppm_to_index(ppm_y, pos_f1)
            x_idx = ppm_to_index(ppm_x, pos_f2)

            patch = extract_patch(spectrum, y_idx, x_idx, patch_size)
            if patch is not None:
                pos_patches.append(patch)

        # Find and extract hard negatives
        neg_locations = find_hard_negatives(
            spectrum, ppm_x, ppm_y, peak_positions,
            max_negatives=max_negatives_per_spectrum,
        )

        neg_patches = []
        for y_idx, x_idx in neg_locations:
            patch = extract_patch(spectrum, y_idx, x_idx, patch_size)
            if patch is not None:
                neg_patches.append(patch)

        if len(pos_patches) > 0 and len(neg_patches) > 0:
            patches = np.array(pos_patches + neg_patches)
            labels = np.array([1.0] * len(pos_patches) + [0.0] * len(neg_patches))
            spectrum_patches[spectrum_name] = (patches, labels)

            logger.info(f"[{i+1}/{len(samples)}] {spectrum_name}: "
                       f"{len(pos_patches)} peaks, {len(neg_patches)} negatives")

    # Spectrum-level split
    spectrum_names = list(spectrum_patches.keys())
    random.shuffle(spectrum_names)

    n_spectra = len(spectrum_names)
    n_test = max(1, int(0.15 * n_spectra))
    n_val = max(1, int(0.15 * n_spectra))

    test_spectra = spectrum_names[:n_test]
    val_spectra = spectrum_names[n_test:n_test + n_val]
    train_spectra = spectrum_names[n_test + n_val:]

    logger.info(f"Split: {len(train_spectra)} train, {len(val_spectra)} val, {len(test_spectra)} test spectra")

    def combine_spectra(names):
        all_patches = []
        all_labels = []
        for name in names:
            patches, labels = spectrum_patches[name]
            all_patches.append(patches)
            all_labels.append(labels)
        if all_patches:
            return {
                'patches': np.concatenate(all_patches),
                'labels': np.concatenate(all_labels),
                'spectra': names,
            }
        return {'patches': np.array([]), 'labels': np.array([]), 'spectra': []}

    return (
        combine_spectra(train_spectra),
        combine_spectra(val_spectra),
        combine_spectra(test_spectra),
    )


def train_model(
    model: nn.Module,
    train_loader: DataLoader,
    val_loader: DataLoader,
    device: torch.device,
    epochs: int = 50,
    learning_rate: float = 1e-3,
    weight_decay: float = 1e-4,
    pos_weight: float = 1.0,
    patience: int = 10,
) -> Tuple[nn.Module, Dict]:
    """
    Train the CNN model.

    Args:
        model: The CNN model
        train_loader: Training data loader
        val_loader: Validation data loader
        device: torch device
        epochs: Maximum epochs
        learning_rate: Initial learning rate
        weight_decay: L2 regularization
        pos_weight: Weight for positive class (for imbalance)
        patience: Early stopping patience

    Returns:
        Tuple of (best_model, training_history)
    """
    model = model.to(device)

    # Loss with class weighting
    criterion = nn.BCEWithLogitsLoss(pos_weight=torch.tensor([pos_weight]).to(device))

    # Optimizer
    optimizer = optim.AdamW(model.parameters(), lr=learning_rate, weight_decay=weight_decay)

    # Learning rate scheduler
    scheduler = optim.lr_scheduler.CosineAnnealingLR(optimizer, T_max=epochs)

    history = {
        'train_loss': [],
        'val_loss': [],
        'train_acc': [],
        'val_acc': [],
    }

    best_val_loss = float('inf')
    best_model_state = None
    epochs_without_improvement = 0

    for epoch in range(epochs):
        # Training
        model.train()
        train_loss = 0.0
        train_correct = 0
        train_total = 0

        for patches, labels in train_loader:
            patches = patches.to(device)
            labels = labels.to(device)

            optimizer.zero_grad()
            outputs = model(patches).squeeze(-1)
            loss = criterion(outputs, labels)
            loss.backward()
            optimizer.step()

            train_loss += loss.item() * patches.size(0)
            preds = (torch.sigmoid(outputs) >= 0.5).float()
            train_correct += (preds == labels).sum().item()
            train_total += labels.size(0)

        train_loss /= train_total
        train_acc = train_correct / train_total

        # Validation
        model.eval()
        val_loss = 0.0
        val_correct = 0
        val_total = 0

        with torch.no_grad():
            for patches, labels in val_loader:
                patches = patches.to(device)
                labels = labels.to(device)

                outputs = model(patches).squeeze(-1)
                loss = criterion(outputs, labels)

                val_loss += loss.item() * patches.size(0)
                preds = (torch.sigmoid(outputs) >= 0.5).float()
                val_correct += (preds == labels).sum().item()
                val_total += labels.size(0)

        val_loss /= val_total
        val_acc = val_correct / val_total

        scheduler.step()

        # Record history
        history['train_loss'].append(train_loss)
        history['val_loss'].append(val_loss)
        history['train_acc'].append(train_acc)
        history['val_acc'].append(val_acc)

        # Early stopping
        if val_loss < best_val_loss:
            best_val_loss = val_loss
            best_model_state = model.state_dict().copy()
            epochs_without_improvement = 0
            marker = " *"
        else:
            epochs_without_improvement += 1
            marker = ""

        logger.info(
            f"Epoch {epoch+1}/{epochs}: "
            f"train_loss={train_loss:.4f}, train_acc={train_acc:.4f}, "
            f"val_loss={val_loss:.4f}, val_acc={val_acc:.4f}{marker}"
        )

        if epochs_without_improvement >= patience:
            logger.info(f"Early stopping after {epoch+1} epochs")
            break

    # Load best model
    if best_model_state is not None:
        model.load_state_dict(best_model_state)

    return model, history


def main():
    parser = argparse.ArgumentParser(description='Train CNN peak classifier')
    parser.add_argument('--training-dir', type=Path, required=True,
                       help='Directory containing training_data.json')
    parser.add_argument('--spectra-dirs', type=Path, nargs='+', required=True,
                       help='Directories to search for spectrum files')
    parser.add_argument('--output', type=Path, required=True,
                       help='Output path for trained model (.pt)')
    parser.add_argument('--min-r2', type=float, default=0.85,
                       help='Minimum R² for positive samples (default: 0.85)')
    parser.add_argument('--patch-size', type=int, default=21,
                       help='Patch size in pixels (default: 21)')
    parser.add_argument('--epochs', type=int, default=50,
                       help='Maximum training epochs (default: 50)')
    parser.add_argument('--batch-size', type=int, default=64,
                       help='Batch size (default: 64)')
    parser.add_argument('--learning-rate', type=float, default=1e-3,
                       help='Learning rate (default: 1e-3)')
    parser.add_argument('--patience', type=int, default=10,
                       help='Early stopping patience (default: 10)')
    parser.add_argument('--no-augment', action='store_true',
                       help='Disable data augmentation')
    parser.add_argument('--seed', type=int, default=42,
                       help='Random seed (default: 42)')

    args = parser.parse_args()

    # Set random seeds
    random.seed(args.seed)
    np.random.seed(args.seed)
    torch.manual_seed(args.seed)

    # Determine device
    if torch.cuda.is_available():
        device = torch.device('cuda')
    elif hasattr(torch.backends, 'mps') and torch.backends.mps.is_available():
        device = torch.device('mps')
    else:
        device = torch.device('cpu')
    logger.info(f"Using device: {device}")

    # Load training data
    training_data_path = args.training_dir / 'training_data.json'
    if not training_data_path.exists():
        logger.error(f"Training data not found: {training_data_path}")
        sys.exit(1)

    logger.info(f"Loading training data from {training_data_path}")
    with open(training_data_path) as f:
        training_data = json.load(f)

    # Extract patches
    logger.info("Extracting patches from spectra...")
    t_start = time.time()
    train_data, val_data, test_data = extract_training_data(
        training_data,
        args.spectra_dirs,
        patch_size=args.patch_size,
        min_r2=args.min_r2,
    )
    t_extract = time.time() - t_start

    logger.info(f"Extraction complete in {t_extract:.1f}s")
    logger.info(f"Train: {len(train_data['labels'])} samples "
               f"({sum(train_data['labels'])} peaks, {len(train_data['labels']) - sum(train_data['labels']):.0f} negatives)")
    logger.info(f"Val: {len(val_data['labels'])} samples")
    logger.info(f"Test: {len(test_data['labels'])} samples")

    if len(train_data['labels']) < 100:
        logger.error("Not enough training samples!")
        sys.exit(1)

    # Create datasets
    train_dataset = PeakPatchDataset(
        train_data['patches'],
        train_data['labels'],
        augment=not args.no_augment,
    )
    val_dataset = PeakPatchDataset(
        val_data['patches'],
        val_data['labels'],
        augment=False,
    )

    # Create data loaders
    train_loader = DataLoader(
        train_dataset,
        batch_size=args.batch_size,
        shuffle=True,
        num_workers=0,  # Avoid multiprocessing issues
        pin_memory=device.type == 'cuda',
    )
    val_loader = DataLoader(
        val_dataset,
        batch_size=args.batch_size,
        shuffle=False,
        num_workers=0,
        pin_memory=device.type == 'cuda',
    )

    # Calculate class weight
    n_pos = sum(train_data['labels'])
    n_neg = len(train_data['labels']) - n_pos
    pos_weight = n_neg / n_pos if n_pos > 0 else 1.0
    logger.info(f"Class balance: {n_pos:.0f} positive, {n_neg:.0f} negative (pos_weight={pos_weight:.2f})")

    # Create model
    config = CNNConfig(patch_size=args.patch_size)
    model = create_peak_cnn(config)
    n_params = count_parameters(model)
    logger.info(f"Created PeakCNN with {n_params:,} parameters")

    # Train
    logger.info("Starting training...")
    t_start = time.time()
    model, history = train_model(
        model,
        train_loader,
        val_loader,
        device,
        epochs=args.epochs,
        learning_rate=args.learning_rate,
        pos_weight=pos_weight,
        patience=args.patience,
    )
    t_train = time.time() - t_start
    logger.info(f"Training complete in {t_train:.1f}s")

    # Save model
    args.output.parent.mkdir(parents=True, exist_ok=True)

    checkpoint = {
        'model_state_dict': model.state_dict(),
        'config': {
            'patch_size': config.patch_size,
            'n_channels': config.n_channels,
            'dropout': config.dropout,
            'global_mean': 0.0,
            'global_std': 1.0,
        },
        'metadata': {
            'train_spectra': train_data['spectra'],
            'val_spectra': val_data['spectra'],
            'test_spectra': test_data['spectra'],
            'n_train_samples': len(train_data['labels']),
            'n_val_samples': len(val_data['labels']),
            'n_test_samples': len(test_data['labels']),
            'best_val_loss': min(history['val_loss']),
            'best_val_acc': max(history['val_acc']),
            'epochs_trained': len(history['train_loss']),
            'pos_weight': pos_weight,
        },
        'history': history,
        'version': '1.0.0',
    }

    torch.save(checkpoint, args.output)
    logger.info(f"Saved model to {args.output}")

    # Print summary
    print("\n" + "=" * 60)
    print("TRAINING SUMMARY")
    print("=" * 60)
    print(f"Model parameters: {n_params:,}")
    print(f"Training samples: {len(train_data['labels']):,}")
    print(f"Validation samples: {len(val_data['labels']):,}")
    print(f"Test samples: {len(test_data['labels']):,}")
    print(f"Best validation loss: {min(history['val_loss']):.4f}")
    print(f"Best validation accuracy: {max(history['val_acc']):.4f}")
    print(f"Epochs trained: {len(history['train_loss'])}")
    print(f"Training time: {t_train:.1f}s")
    print("=" * 60)


if __name__ == '__main__':
    main()

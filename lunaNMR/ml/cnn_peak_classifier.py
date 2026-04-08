# ABOUTME: CNN-based peak classifier for NMR spectra using PyTorch.
# ABOUTME: Classifies 2D spectrum patches as peak or noise.

"""
CNN Peak Classifier for NMR Spectra

This module provides a small CNN that classifies 2D spectrum patches as peaks or noise.
The CNN learns discriminative features directly from raw data instead of hand-crafted features.

Architecture: 3 conv layers (~25K params) → fast inference (~0.6ms per candidate on CPU)
"""

import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union

import numpy as np

logger = logging.getLogger(__name__)

# Lazy import PyTorch to avoid import errors if not installed
_torch = None
_nn = None


def _import_torch():
    """Lazy import of PyTorch."""
    global _torch, _nn
    if _torch is None:
        try:
            import torch
            import torch.nn as nn
            _torch = torch
            _nn = nn
        except ImportError:
            raise ImportError(
                "PyTorch is required for CNN peak classification. "
                "Install with: pip install torch torchvision"
            )
    return _torch, _nn


@dataclass
class CNNConfig:
    """Configuration for CNN peak classifier."""
    patch_size: int = 21  # Patch size in pixels (21x21)
    n_channels: int = 1   # Single channel (intensity)
    dropout: float = 0.3  # Dropout rate

    # Normalization parameters (learned from training data)
    global_mean: float = 0.0
    global_std: float = 1.0


def create_peak_cnn(config: Optional[CNNConfig] = None):
    """
    Create the PeakCNN model.

    Architecture:
        Input: (batch, 1, 21, 21)
        Conv2D(1→16, 3x3) → BatchNorm → ReLU → MaxPool(2x2)
        Conv2D(16→32, 3x3) → BatchNorm → ReLU → MaxPool(2x2)
        Conv2D(32→64, 3x3) → BatchNorm → ReLU → AdaptiveAvgPool(1x1)
        Flatten → Linear(64→32) → ReLU → Dropout → Linear(32→1) → Sigmoid

    Returns:
        nn.Module: The CNN model
    """
    torch, nn = _import_torch()

    if config is None:
        config = CNNConfig()

    class PeakCNN(nn.Module):
        """Small CNN for peak classification from 2D patches."""

        def __init__(self):
            super().__init__()

            # Convolutional layers
            self.conv1 = nn.Sequential(
                nn.Conv2d(1, 16, kernel_size=3, padding=1),
                nn.BatchNorm2d(16),
                nn.ReLU(inplace=True),
                nn.MaxPool2d(2, 2)  # 21→10
            )

            self.conv2 = nn.Sequential(
                nn.Conv2d(16, 32, kernel_size=3, padding=1),
                nn.BatchNorm2d(32),
                nn.ReLU(inplace=True),
                nn.MaxPool2d(2, 2)  # 10→5
            )

            self.conv3 = nn.Sequential(
                nn.Conv2d(32, 64, kernel_size=3, padding=1),
                nn.BatchNorm2d(64),
                nn.ReLU(inplace=True),
                nn.AdaptiveAvgPool2d(1)  # 5→1
            )

            # Classifier
            self.classifier = nn.Sequential(
                nn.Flatten(),
                nn.Linear(64, 32),
                nn.ReLU(inplace=True),
                nn.Dropout(config.dropout),
                nn.Linear(32, 1),
            )

        def forward(self, x):
            """Forward pass."""
            x = self.conv1(x)
            x = self.conv2(x)
            x = self.conv3(x)
            x = self.classifier(x)
            return x

        def get_embedding(self, x):
            """Get 64-dim embedding (for hybrid models)."""
            x = self.conv1(x)
            x = self.conv2(x)
            x = self.conv3(x)
            return x.view(x.size(0), -1)  # (batch, 64)

    return PeakCNN()


class CNNPeakClassifier:
    """
    Wrapper class for CNN peak classification.

    Handles:
    - Model loading/saving
    - Patch extraction and normalization
    - Batch inference
    """

    def __init__(
        self,
        model_path: Optional[Union[str, Path]] = None,
        device: Optional[str] = None,
    ):
        """
        Initialize the classifier.

        Args:
            model_path: Path to trained model (.pt file)
            device: 'cpu', 'cuda', or 'mps'. Auto-detected if None.
        """
        torch, nn = _import_torch()

        # Determine device
        if device is None:
            if torch.cuda.is_available():
                device = 'cuda'
            elif hasattr(torch.backends, 'mps') and torch.backends.mps.is_available():
                device = 'mps'
            else:
                device = 'cpu'

        self.device = torch.device(device)
        self.model = None
        self.config = CNNConfig()

        if model_path is not None:
            self.load(model_path)

    def load(self, model_path: Union[str, Path]) -> None:
        """
        Load trained model from file.

        Args:
            model_path: Path to .pt file containing model state and config
        """
        torch, nn = _import_torch()

        model_path = Path(model_path)
        if not model_path.exists():
            raise FileNotFoundError(f"Model not found: {model_path}")

        checkpoint = torch.load(model_path, map_location=self.device)

        # Load config
        if 'config' in checkpoint:
            config_dict = checkpoint['config']
            self.config = CNNConfig(**config_dict)

        # Create and load model
        self.model = create_peak_cnn(self.config)
        self.model.load_state_dict(checkpoint['model_state_dict'])
        self.model.to(self.device)
        self.model.eval()

        logger.info(f"Loaded CNN model from {model_path} (device: {self.device})")

    def save(self, model_path: Union[str, Path], metadata: Optional[Dict] = None) -> None:
        """
        Save model to file.

        Args:
            model_path: Output path for .pt file
            metadata: Optional metadata to include (e.g., training info)
        """
        torch, _ = _import_torch()

        if self.model is None:
            raise RuntimeError("No model to save")

        checkpoint = {
            'model_state_dict': self.model.state_dict(),
            'config': {
                'patch_size': self.config.patch_size,
                'n_channels': self.config.n_channels,
                'dropout': self.config.dropout,
                'global_mean': self.config.global_mean,
                'global_std': self.config.global_std,
            },
            'version': '1.0.0',
        }

        if metadata:
            checkpoint['metadata'] = metadata

        torch.save(checkpoint, model_path)
        logger.info(f"Saved CNN model to {model_path}")

    def extract_patches(
        self,
        spectrum: np.ndarray,
        candidates: List[Dict],
        normalize: bool = True,
    ) -> np.ndarray:
        """
        Extract patches at candidate locations.

        Args:
            spectrum: 2D spectrum array (ny, nx)
            candidates: List of dicts with 'y_idx', 'x_idx' keys
            normalize: Whether to z-score normalize each patch

        Returns:
            patches: Array of shape (n_candidates, 1, patch_size, patch_size)
        """
        patch_size = self.config.patch_size
        half = patch_size // 2
        ny, nx = spectrum.shape

        patches = []
        for cand in candidates:
            y_idx = cand['y_idx']
            x_idx = cand['x_idx']

            # Extract patch with zero-padding at edges
            y_start = max(0, y_idx - half)
            y_end = min(ny, y_idx + half + 1)
            x_start = max(0, x_idx - half)
            x_end = min(nx, x_idx + half + 1)

            # Create zero-padded patch
            patch = np.zeros((patch_size, patch_size), dtype=np.float32)

            # Calculate where to place the extracted region
            py_start = half - (y_idx - y_start)
            py_end = py_start + (y_end - y_start)
            px_start = half - (x_idx - x_start)
            px_end = px_start + (x_end - x_start)

            patch[py_start:py_end, px_start:px_end] = spectrum[y_start:y_end, x_start:x_end]

            # Normalize
            if normalize:
                mean = patch.mean()
                std = patch.std()
                if std > 1e-10:
                    patch = (patch - mean) / std
                else:
                    patch = patch - mean

            patches.append(patch)

        # Shape: (n_candidates, 1, patch_size, patch_size)
        return np.array(patches, dtype=np.float32)[:, np.newaxis, :, :]

    def classify(
        self,
        patches: np.ndarray,
        threshold: float = 0.5,
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Classify patches as peak or noise.

        Args:
            patches: Array of shape (n, 1, patch_size, patch_size)
            threshold: Classification threshold (default 0.5)

        Returns:
            is_peak: Boolean array of shape (n,)
            confidence: Probability array of shape (n,)
        """
        torch, _ = _import_torch()

        if self.model is None:
            raise RuntimeError("Model not loaded. Call load() first.")

        # Convert to tensor (use torch.tensor for numpy 2.x compatibility)
        patches_tensor = torch.tensor(patches, dtype=torch.float32).to(self.device)

        # Inference
        with torch.no_grad():
            logits = self.model(patches_tensor)
            probs = torch.sigmoid(logits).squeeze(-1)

        # Convert to numpy via list to avoid numpy 2.x compatibility issues
        confidence = np.array(probs.cpu().tolist(), dtype=np.float32)
        is_peak = confidence >= threshold

        return is_peak, confidence

    def classify_spectrum(
        self,
        spectrum: np.ndarray,
        candidates: List[Dict],
        threshold: float = 0.5,
        batch_size: int = 256,
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Classify candidates in a spectrum.

        Args:
            spectrum: 2D spectrum array
            candidates: List of candidate dicts with 'y_idx', 'x_idx'
            threshold: Classification threshold
            batch_size: Batch size for inference

        Returns:
            is_peak: Boolean array
            confidence: Probability array
        """
        torch, _ = _import_torch()

        if len(candidates) == 0:
            return np.array([], dtype=bool), np.array([], dtype=np.float32)

        # Extract all patches
        patches = self.extract_patches(spectrum, candidates)

        # Classify in batches
        n = len(candidates)
        all_probs = []

        for i in range(0, n, batch_size):
            batch = patches[i:i + batch_size]
            _, probs = self.classify(batch, threshold=0.0)  # Get raw probs
            all_probs.append(probs)

        confidence = np.concatenate(all_probs)
        is_peak = confidence >= threshold

        return is_peak, confidence


def count_parameters(model) -> int:
    """Count trainable parameters in model."""
    return sum(p.numel() for p in model.parameters() if p.requires_grad)


# Test block
if __name__ == "__main__":
    # Test model creation and parameter count
    model = create_peak_cnn()
    n_params = count_parameters(model)
    print(f"PeakCNN parameters: {n_params:,}")

    # Test forward pass
    import torch
    x = torch.randn(4, 1, 21, 21)
    out = model(x)
    print(f"Input shape: {x.shape}")
    print(f"Output shape: {out.shape}")

    # Test embedding
    emb = model.get_embedding(x)
    print(f"Embedding shape: {emb.shape}")

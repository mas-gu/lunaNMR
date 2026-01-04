# ABOUTME: Random Forest classifier wrapper for peak vs noise classification.
# ABOUTME: Provides training, prediction, and model persistence for ML-based peak detection.

"""
Peak Classifier for ML-based peak detection.

Wraps scikit-learn RandomForestClassifier with:
- Feature scaling
- Confidence-based classification
- Model persistence (joblib)
- Feature importance analysis
"""

import numpy as np
from typing import Tuple, Optional, List, Dict
from pathlib import Path
import logging

# Try to import sklearn, gracefully handle if not available
try:
    from sklearn.ensemble import RandomForestClassifier
    from sklearn.preprocessing import StandardScaler
    from sklearn.model_selection import cross_val_score
    import joblib
    SKLEARN_AVAILABLE = True
except ImportError:
    SKLEARN_AVAILABLE = False
    logging.warning("scikit-learn not available. Peak classifier will not work.")

from .peak_feature_extractor import PeakFeatures, PeakFeatureExtractor


class PeakClassifier:
    """
    Random Forest classifier for peak vs noise classification.

    Features:
    - StandardScaler for feature normalization
    - Balanced class weights (handles imbalanced data)
    - Probability output for confidence thresholding
    - Feature importance analysis
    """

    def __init__(
        self,
        n_estimators: int = 200,
        max_depth: int = 15,
        min_samples_leaf: int = 3,
        random_state: int = 42,
    ):
        """
        Initialize classifier.

        Args:
            n_estimators: Number of trees in forest
            max_depth: Maximum depth of trees
            min_samples_leaf: Minimum samples required in leaf
            random_state: Random seed for reproducibility
        """
        if not SKLEARN_AVAILABLE:
            raise ImportError("scikit-learn is required for PeakClassifier")

        self.scaler = StandardScaler()
        self.model = RandomForestClassifier(
            n_estimators=n_estimators,
            max_depth=max_depth,
            min_samples_leaf=min_samples_leaf,
            class_weight='balanced',
            n_jobs=-1,
            random_state=random_state,
        )
        self.is_fitted = False
        self.feature_names = PeakFeatures.feature_names()

    def fit(
        self,
        X_train: np.ndarray,
        y_train: np.ndarray,
        X_val: Optional[np.ndarray] = None,
        y_val: Optional[np.ndarray] = None,
    ) -> Dict:
        """
        Train the classifier.

        Args:
            X_train: Training features, shape (n_samples, n_features)
            y_train: Training labels, shape (n_samples,), 1=peak, 0=noise
            X_val: Optional validation features
            y_val: Optional validation labels

        Returns:
            Dict with training statistics
        """
        # Scale features
        X_train_scaled = self.scaler.fit_transform(X_train)

        # Train model
        self.model.fit(X_train_scaled, y_train)
        self.is_fitted = True

        # Compute statistics
        stats = {
            'n_train_samples': len(y_train),
            'n_peaks': int(np.sum(y_train == 1)),
            'n_noise': int(np.sum(y_train == 0)),
            'train_accuracy': float(self.model.score(X_train_scaled, y_train)),
        }

        # Validation metrics if provided
        if X_val is not None and y_val is not None:
            X_val_scaled = self.scaler.transform(X_val)
            stats['val_accuracy'] = float(self.model.score(X_val_scaled, y_val))

            # Compute precision/recall on validation
            y_pred = self.model.predict(X_val_scaled)
            tp = np.sum((y_pred == 1) & (y_val == 1))
            fp = np.sum((y_pred == 1) & (y_val == 0))
            fn = np.sum((y_pred == 0) & (y_val == 1))

            precision = tp / (tp + fp) if (tp + fp) > 0 else 0.0
            recall = tp / (tp + fn) if (tp + fn) > 0 else 0.0
            f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0.0

            stats['val_precision'] = float(precision)
            stats['val_recall'] = float(recall)
            stats['val_f1'] = float(f1)

        return stats

    def cross_validate(
        self,
        X: np.ndarray,
        y: np.ndarray,
        cv: int = 5,
    ) -> Dict:
        """
        Perform cross-validation.

        Args:
            X: Features
            y: Labels
            cv: Number of folds

        Returns:
            Dict with CV statistics
        """
        X_scaled = self.scaler.fit_transform(X)
        scores = cross_val_score(self.model, X_scaled, y, cv=cv, scoring='accuracy')

        return {
            'cv_mean_accuracy': float(np.mean(scores)),
            'cv_std_accuracy': float(np.std(scores)),
            'cv_scores': scores.tolist(),
        }

    def predict_proba(self, X: np.ndarray) -> np.ndarray:
        """
        Predict probability of being a peak.

        Args:
            X: Features, shape (n_samples, n_features)

        Returns:
            Probabilities, shape (n_samples,)
        """
        if not self.is_fitted:
            raise RuntimeError("Model not fitted. Call fit() first.")

        X_scaled = self.scaler.transform(X)
        probs = self.model.predict_proba(X_scaled)

        # Return probability of class 1 (peak)
        return probs[:, 1]

    def classify(
        self,
        X: np.ndarray,
        threshold: float = 0.5,
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Classify candidates as peak or noise.

        Args:
            X: Features, shape (n_samples, n_features)
            threshold: Classification threshold (default 0.5)

        Returns:
            Tuple of (is_peak, confidence) arrays
        """
        probs = self.predict_proba(X)
        is_peak = probs >= threshold
        return is_peak, probs

    def get_feature_importance(self) -> Dict[str, float]:
        """
        Get feature importance from trained model.

        Returns:
            Dict mapping feature name to importance
        """
        if not self.is_fitted:
            raise RuntimeError("Model not fitted. Call fit() first.")

        importances = self.model.feature_importances_
        return dict(zip(self.feature_names, importances.tolist()))

    def save(self, path: Path) -> None:
        """
        Save model to file.

        Args:
            path: Path to save model (.joblib)
        """
        if not self.is_fitted:
            raise RuntimeError("Model not fitted. Call fit() first.")

        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)

        model_data = {
            'scaler': self.scaler,
            'model': self.model,
            'feature_names': self.feature_names,
            'version': '1.0.0',
        }
        joblib.dump(model_data, path)

    @classmethod
    def load(cls, path: Path) -> 'PeakClassifier':
        """
        Load model from file.

        Args:
            path: Path to saved model (.joblib)

        Returns:
            Loaded PeakClassifier instance
        """
        if not SKLEARN_AVAILABLE:
            raise ImportError("scikit-learn is required to load PeakClassifier")

        model_data = joblib.load(path)

        classifier = cls.__new__(cls)
        classifier.scaler = model_data['scaler']
        classifier.model = model_data['model']
        classifier.feature_names = model_data.get('feature_names', PeakFeatures.feature_names())
        classifier.is_fitted = True

        return classifier


class MLPeakDetector:
    """
    High-level interface for ML-based peak detection.

    Combines feature extraction and classification.
    """

    def __init__(
        self,
        classifier: Optional[PeakClassifier] = None,
        model_path: Optional[Path] = None,
        noise_level: float = 1.0,
        confidence_threshold: float = 0.4,
    ):
        """
        Initialize ML peak detector.

        Args:
            classifier: Trained PeakClassifier (or load from model_path)
            model_path: Path to load classifier from
            noise_level: Noise level for feature extraction
            confidence_threshold: Classification threshold (aggressive default)
        """
        if classifier is not None:
            self.classifier = classifier
        elif model_path is not None:
            self.classifier = PeakClassifier.load(model_path)
        else:
            raise ValueError("Must provide classifier or model_path")

        self.feature_extractor = PeakFeatureExtractor(noise_level=noise_level)
        self.confidence_threshold = confidence_threshold

    def detect_peaks(
        self,
        spectrum: np.ndarray,
        candidates: List[Dict],
    ) -> List[Dict]:
        """
        Classify candidate peaks as real or noise.

        Args:
            spectrum: 2D spectrum data
            candidates: List of candidate dicts with 'y_idx' and 'x_idx'

        Returns:
            List of accepted peaks (classified as real)
        """
        if not candidates:
            return []

        # Extract features for all candidates
        features = self.feature_extractor.extract_features_batch(spectrum, candidates)

        # Classify
        is_peak, confidence = self.classifier.classify(features, self.confidence_threshold)

        # Return accepted peaks
        accepted = []
        for i, candidate in enumerate(candidates):
            if is_peak[i]:
                candidate = candidate.copy()
                candidate['ml_confidence'] = float(confidence[i])
                accepted.append(candidate)

        return accepted

    def update_noise_level(self, noise_level: float) -> None:
        """Update noise level for feature extraction."""
        self.feature_extractor.noise_level = max(noise_level, 1e-10)

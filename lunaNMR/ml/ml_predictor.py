# ABOUTME: Random Forest ML predictor for NMR parameter optimization
# ABOUTME: Wraps scikit-learn with confidence estimation and graceful fallback

from pathlib import Path
from typing import List, Optional, Tuple
import numpy as np
import logging

from .feature_extractor import SpectrumFeatures, OptimalParameters

# Try to import scikit-learn, gracefully handle if not available
try:
    from sklearn.ensemble import RandomForestRegressor
    import joblib
    SKLEARN_AVAILABLE = True
except ImportError:
    SKLEARN_AVAILABLE = False
    RandomForestRegressor = None
    joblib = None


logger = logging.getLogger(__name__)


class MLPredictor:
    """
    Random Forest predictor for NMR processing parameters.

    Predicts optimal linewidths, integration windows, and overlap thresholds
    based on spectrum features. Includes confidence estimation based on
    tree agreement.
    """

    # Model configuration
    N_ESTIMATORS = 100
    MAX_DEPTH = 10
    MIN_SAMPLES_LEAF = 5
    RANDOM_STATE = 42

    def __init__(self):
        """Initialize the ML predictor."""
        self.model: Optional["RandomForestRegressor"] = None
        self._is_trained = False

        if not SKLEARN_AVAILABLE:
            logger.warning(
                "scikit-learn not installed. ML predictions will not be available. "
                "Install with: pip install scikit-learn"
            )

    def is_available(self) -> bool:
        """Check if model is available for predictions."""
        return SKLEARN_AVAILABLE and self._is_trained and self.model is not None

    def train(
        self,
        features_list: List[SpectrumFeatures],
        targets_list: List[OptimalParameters]
    ) -> bool:
        """
        Train the Random Forest model on collected data.

        Parameters
        ----------
        features_list : list of SpectrumFeatures
            Input features for each training sample
        targets_list : list of OptimalParameters
            Target parameters for each training sample

        Returns
        -------
        bool
            True if training succeeded, False otherwise
        """
        if not SKLEARN_AVAILABLE:
            logger.error("Cannot train: scikit-learn not installed")
            return False

        if len(features_list) < 10:
            logger.warning(f"Insufficient training data: {len(features_list)} samples (need >= 10)")
            return False

        if len(features_list) != len(targets_list):
            logger.error("Features and targets lists must have same length")
            return False

        try:
            # Convert to numpy arrays
            X = np.array([f.to_array() for f in features_list])
            y = np.array([t.to_array() for t in targets_list])

            # Create and train model
            self.model = RandomForestRegressor(
                n_estimators=self.N_ESTIMATORS,
                max_depth=self.MAX_DEPTH,
                min_samples_leaf=self.MIN_SAMPLES_LEAF,
                n_jobs=-1,  # Use all cores
                random_state=self.RANDOM_STATE,
            )

            self.model.fit(X, y)
            self._is_trained = True

            logger.info(f"ML model trained on {len(features_list)} samples")
            return True

        except Exception as e:
            logger.error(f"Training failed: {e}")
            return False

    def predict(self, features: SpectrumFeatures) -> Optional[OptimalParameters]:
        """
        Predict optimal parameters for a spectrum.

        Parameters
        ----------
        features : SpectrumFeatures
            Features extracted from the spectrum

        Returns
        -------
        OptimalParameters or None
            Predicted parameters with confidence, or None if unavailable
        """
        if not self.is_available():
            return None

        try:
            X = features.to_array().reshape(1, -1)

            # Get prediction and confidence
            prediction, confidence = self._predict_with_confidence(X)

            return OptimalParameters.from_array(
                prediction,
                confidence=confidence,
                source="ml"
            )

        except Exception as e:
            logger.error(f"Prediction failed: {e}")
            return None

    def _predict_with_confidence(self, X: np.ndarray) -> Tuple[np.ndarray, float]:
        """
        Make prediction with confidence estimation.

        Confidence is based on agreement between trees in the forest.
        Higher agreement = higher confidence.

        Parameters
        ----------
        X : np.ndarray
            Input features (1, n_features)

        Returns
        -------
        tuple
            (prediction array, confidence float)
        """
        # Get predictions from all trees
        tree_predictions = np.array([
            tree.predict(X) for tree in self.model.estimators_
        ])  # Shape: (n_trees, 1, n_targets)

        # Mean prediction across trees
        prediction = np.mean(tree_predictions, axis=0).flatten()

        # Calculate confidence from coefficient of variation
        std = np.std(tree_predictions, axis=0).flatten()
        cv = std / (np.abs(prediction) + 1e-6)  # Avoid division by zero

        # Confidence is inverse of mean CV, scaled to 0-1
        mean_cv = np.mean(cv)
        confidence = 1.0 / (1.0 + mean_cv)

        # Clip to valid range
        confidence = float(np.clip(confidence, 0.0, 1.0))

        return prediction, confidence

    def save(self, path: Path) -> bool:
        """
        Save trained model to file.

        Parameters
        ----------
        path : Path
            Output file path (.joblib extension recommended)

        Returns
        -------
        bool
            True if save succeeded
        """
        if not self.is_available():
            logger.error("Cannot save: no trained model")
            return False

        if not SKLEARN_AVAILABLE:
            logger.error("Cannot save: joblib not available")
            return False

        try:
            path = Path(path)
            path.parent.mkdir(parents=True, exist_ok=True)
            joblib.dump(self.model, path)
            logger.info(f"Model saved to {path}")
            return True

        except Exception as e:
            logger.error(f"Failed to save model: {e}")
            return False

    def load(self, path: Path) -> bool:
        """
        Load trained model from file.

        Parameters
        ----------
        path : Path
            Input file path (.joblib)

        Returns
        -------
        bool
            True if load succeeded
        """
        if not SKLEARN_AVAILABLE:
            logger.error("Cannot load: scikit-learn not installed")
            return False

        try:
            path = Path(path)
            if not path.exists():
                logger.error(f"Model file not found: {path}")
                return False

            self.model = joblib.load(path)
            self._is_trained = True

            logger.info(f"Model loaded from {path}")
            return True

        except Exception as e:
            logger.error(f"Failed to load model: {e}")
            return False

    def get_feature_importance(self) -> Optional[dict]:
        """
        Get feature importance from trained model.

        Returns
        -------
        dict or None
            Dictionary mapping feature names to importance scores
        """
        if not self.is_available():
            return None

        try:
            importances = self.model.feature_importances_
            names = SpectrumFeatures.field_names()

            # Average importance across all targets
            if importances.ndim > 1:
                importances = np.mean(importances, axis=0)

            return dict(zip(names, importances.tolist()))

        except Exception as e:
            logger.error(f"Failed to get feature importance: {e}")
            return None

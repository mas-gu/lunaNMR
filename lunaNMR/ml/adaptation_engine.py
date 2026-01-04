# ABOUTME: Adaptation engine for local ML model retraining
# ABOUTME: Uses collected training data to adapt models to user's spectra

from typing import Dict, Optional, Any
from datetime import datetime
import logging

from .feature_extractor import SpectrumFeatures, OptimalParameters
from .ml_predictor import MLPredictor, SKLEARN_AVAILABLE
from .training_collector import TrainingCollector


logger = logging.getLogger(__name__)


class AdaptationEngine:
    """
    Adapts pre-trained ML models using locally collected training data.

    Uses a weighted ensemble approach: blends predictions from the base
    (pre-trained) model with a locally trained model. The blend weight
    increases as more local samples are collected.
    """

    # Default configuration
    DEFAULT_MIN_SAMPLES = 10  # Minimum samples to train local model
    DEFAULT_MAX_LOCAL_WEIGHT = 0.7  # Maximum weight for local model
    DEFAULT_SAMPLES_FOR_MAX_WEIGHT = 100  # Samples needed for max weight

    def __init__(
        self,
        min_samples_for_adaptation: int = DEFAULT_MIN_SAMPLES,
        max_local_weight: float = DEFAULT_MAX_LOCAL_WEIGHT,
        samples_for_max_weight: int = DEFAULT_SAMPLES_FOR_MAX_WEIGHT,
    ):
        """
        Initialize the adaptation engine.

        Parameters
        ----------
        min_samples_for_adaptation : int
            Minimum samples required to train local model
        max_local_weight : float
            Maximum weight for local model in blending (0-1)
        samples_for_max_weight : int
            Number of samples needed to reach max weight
        """
        self.min_samples_for_adaptation = min_samples_for_adaptation
        self.max_local_weight = max_local_weight
        self.samples_for_max_weight = samples_for_max_weight

        # Local models by nucleus
        self._local_models: Dict[str, Optional[MLPredictor]] = {
            "15N": None,
            "13C": None,
        }

        # Adaptation history
        self._adaptation_history: Dict[str, Optional[Dict]] = {
            "15N": None,
            "13C": None,
        }

    def should_adapt(self, nucleus: str, collector: TrainingCollector) -> bool:
        """
        Check if adaptation should be triggered.

        Parameters
        ----------
        nucleus : str
            Nucleus type
        collector : TrainingCollector
            Training data collector

        Returns
        -------
        bool
            True if adaptation should proceed
        """
        sample_count = collector.get_sample_count(nucleus)
        return sample_count >= self.min_samples_for_adaptation

    def adapt(
        self,
        nucleus: str,
        collector: TrainingCollector,
    ) -> Optional[MLPredictor]:
        """
        Train a local model from collected data.

        Parameters
        ----------
        nucleus : str
            Nucleus type ("15N" or "13C")
        collector : TrainingCollector
            Training data collector

        Returns
        -------
        MLPredictor or None
            Trained local model, or None if training failed
        """
        if not SKLEARN_AVAILABLE:
            logger.warning("Cannot adapt: scikit-learn not installed")
            return None

        sample_count = collector.get_sample_count(nucleus)
        if sample_count < self.min_samples_for_adaptation:
            logger.info(
                f"Not enough samples for {nucleus} adaptation: "
                f"{sample_count} < {self.min_samples_for_adaptation}"
            )
            return None

        # Get training data
        features_list, targets_list = collector.get_training_data(nucleus)

        if not features_list:
            return None

        # Train local model
        local_model = MLPredictor()

        # Use smaller forest for local model
        local_model.N_ESTIMATORS = 50
        local_model.MAX_DEPTH = 8

        if local_model.train(features_list, targets_list):
            self._local_models[nucleus] = local_model
            self.record_adaptation(nucleus, sample_count)
            logger.info(f"Trained local {nucleus} model on {sample_count} samples")
            return local_model
        else:
            logger.error(f"Failed to train local {nucleus} model")
            return None

    def get_local_model(self, nucleus: str) -> Optional[MLPredictor]:
        """Get the locally trained model for a nucleus type."""
        return self._local_models.get(nucleus)

    def get_blend_weight(
        self,
        nucleus: str,
        collector: TrainingCollector,
    ) -> float:
        """
        Calculate blend weight for local model based on sample count.

        Weight increases linearly from 0 to max_local_weight as samples
        increase from min_samples to samples_for_max_weight.

        Parameters
        ----------
        nucleus : str
            Nucleus type
        collector : TrainingCollector
            Training data collector

        Returns
        -------
        float
            Blend weight for local model (0 to max_local_weight)
        """
        sample_count = collector.get_sample_count(nucleus)

        if sample_count < self.min_samples_for_adaptation:
            return 0.0

        # Linear interpolation
        effective_samples = sample_count - self.min_samples_for_adaptation
        max_effective = self.samples_for_max_weight - self.min_samples_for_adaptation

        if max_effective <= 0:
            return self.max_local_weight

        weight = min(1.0, effective_samples / max_effective) * self.max_local_weight
        return weight

    def blend_predictions(
        self,
        base_model: Optional[MLPredictor],
        local_model: Optional[MLPredictor],
        features: SpectrumFeatures,
        local_weight: float,
    ) -> Optional[OptimalParameters]:
        """
        Blend predictions from base and local models.

        Parameters
        ----------
        base_model : MLPredictor or None
            Pre-trained base model
        local_model : MLPredictor or None
            Locally trained model
        features : SpectrumFeatures
            Input features
        local_weight : float
            Weight for local model (0-1)

        Returns
        -------
        OptimalParameters or None
            Blended prediction
        """
        base_pred = None
        local_pred = None

        if base_model is not None and base_model.is_available():
            base_pred = base_model.predict(features)

        if local_model is not None and local_model.is_available():
            local_pred = local_model.predict(features)

        # Handle cases where one or both models unavailable
        if base_pred is None and local_pred is None:
            return None

        if base_pred is None:
            return local_pred

        if local_pred is None:
            return base_pred

        # Blend predictions
        base_weight = 1.0 - local_weight

        blended = OptimalParameters(
            lw_f1_median=base_weight * base_pred.lw_f1_median + local_weight * local_pred.lw_f1_median,
            lw_f2_median=base_weight * base_pred.lw_f2_median + local_weight * local_pred.lw_f2_median,
            rad_f1=base_weight * base_pred.rad_f1 + local_weight * local_pred.rad_f1,
            rad_f2=base_weight * base_pred.rad_f2 + local_weight * local_pred.rad_f2,
            overlap_threshold_f1=base_weight * base_pred.overlap_threshold_f1 + local_weight * local_pred.overlap_threshold_f1,
            overlap_threshold_f2=base_weight * base_pred.overlap_threshold_f2 + local_weight * local_pred.overlap_threshold_f2,
            achievable_r2=base_weight * base_pred.achievable_r2 + local_weight * local_pred.achievable_r2,
            confidence=(base_pred.confidence + local_pred.confidence) / 2,
            source=f"blended_{int(local_weight*100)}pct_local",
        )

        return blended

    def record_adaptation(self, nucleus: str, sample_count: int) -> None:
        """
        Record that adaptation was performed.

        Parameters
        ----------
        nucleus : str
            Nucleus type
        sample_count : int
            Number of samples used
        """
        self._adaptation_history[nucleus] = {
            "timestamp": datetime.now().isoformat(),
            "sample_count": sample_count,
        }

    def get_last_adaptation(self, nucleus: str) -> Optional[Dict[str, Any]]:
        """
        Get information about the last adaptation.

        Parameters
        ----------
        nucleus : str
            Nucleus type

        Returns
        -------
        dict or None
            Adaptation info with timestamp and sample_count
        """
        return self._adaptation_history.get(nucleus)

    def get_status(self) -> Dict[str, Any]:
        """
        Get overall adaptation status.

        Returns
        -------
        dict
            Status information for all nucleus types
        """
        return {
            "15N": {
                "has_local_model": self._local_models.get("15N") is not None,
                "last_adaptation": self._adaptation_history.get("15N"),
            },
            "13C": {
                "has_local_model": self._local_models.get("13C") is not None,
                "last_adaptation": self._adaptation_history.get("13C"),
            },
            "config": {
                "min_samples": self.min_samples_for_adaptation,
                "max_local_weight": self.max_local_weight,
                "samples_for_max_weight": self.samples_for_max_weight,
            },
        }

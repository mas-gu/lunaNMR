# ABOUTME: Model manager coordinating ML and statistical predictors
# ABOUTME: Handles model loading, prediction orchestration, and persistence

from pathlib import Path
from typing import Dict, Optional, Any
import logging

from .feature_extractor import SpectrumFeatures, OptimalParameters
from .ml_predictor import MLPredictor
from .stats_predictor import StatisticalPredictor
from .decision_engine import DecisionEngine


logger = logging.getLogger(__name__)


# Default configuration
DEFAULT_CONFIG = {
    "enabled": False,
    "use_predictions": False,
    "collect_training_data": False,
    "use_statistical_fallback": True,
    "ml_threshold": 0.6,
    "stats_threshold": 0.4,
    "divergence_threshold": 0.3,
}


class ModelManager:
    """
    Coordinate ML and statistical predictors for parameter optimization.

    Manages model loading, prediction selection, and data persistence.
    Acts as the main interface for the ML learning center.
    """

    def __init__(self, config: Optional[Dict[str, Any]] = None):
        """
        Initialize the model manager.

        Parameters
        ----------
        config : dict, optional
            Configuration overrides
        """
        # Merge with defaults
        self.config = {**DEFAULT_CONFIG}
        if config:
            self.config.update(config)

        # Initialize components
        self._ml_predictors: Dict[str, MLPredictor] = {
            "15N": MLPredictor(),
            "13C": MLPredictor(),
        }
        self._stats_predictor = StatisticalPredictor()
        self._decision_engine = DecisionEngine(
            ml_threshold=self.config["ml_threshold"],
            stats_threshold=self.config["stats_threshold"],
            divergence_threshold=self.config["divergence_threshold"],
        )

    def predict(self, features: SpectrumFeatures) -> OptimalParameters:
        """
        Get optimal parameters for a spectrum.

        Uses ML model if available and confident, falls back to
        statistical predictions or defaults.

        Parameters
        ----------
        features : SpectrumFeatures
            Features extracted from the spectrum

        Returns
        -------
        OptimalParameters
            Predicted optimal parameters
        """
        if not self.config["enabled"] or not self.config["use_predictions"]:
            return self._get_default_params(features.nucleus_type)

        # Get ML prediction
        ml_pred = None
        nucleus = features.nucleus_type
        if nucleus in self._ml_predictors:
            ml_predictor = self._ml_predictors[nucleus]
            if ml_predictor.is_available():
                ml_pred = ml_predictor.predict(features)

        # Get statistical prediction
        stats_pred = None
        if self.config["use_statistical_fallback"]:
            stats_pred = self._stats_predictor.predict(features)

        # Use decision engine to select best prediction
        result = self._decision_engine.decide(ml_pred, stats_pred)

        logger.debug(
            f"Prediction for {nucleus}: source={result.source}, "
            f"confidence={result.confidence:.2f}"
        )

        return result

    def update_statistics(self, nucleus: str, params: OptimalParameters) -> None:
        """
        Add successful parameters to statistical history.

        Parameters
        ----------
        nucleus : str
            Nucleus type ("15N" or "13C")
        params : OptimalParameters
            Parameters from successful processing
        """
        if not self.config["collect_training_data"]:
            return

        self._stats_predictor.update(nucleus, params)

    def load_ml_model(self, nucleus: str, path: Path) -> bool:
        """
        Load a pretrained ML model for a nucleus type.

        Parameters
        ----------
        nucleus : str
            Nucleus type ("15N" or "13C")
        path : Path
            Path to .joblib model file

        Returns
        -------
        bool
            True if load succeeded
        """
        if nucleus not in self._ml_predictors:
            self._ml_predictors[nucleus] = MLPredictor()

        return self._ml_predictors[nucleus].load(path)

    def save_ml_model(self, nucleus: str, path: Path) -> bool:
        """
        Save ML model for a nucleus type.

        Parameters
        ----------
        nucleus : str
            Nucleus type
        path : Path
            Output path

        Returns
        -------
        bool
            True if save succeeded
        """
        if nucleus not in self._ml_predictors:
            return False

        return self._ml_predictors[nucleus].save(path)

    def load_statistics(self, path: Path) -> bool:
        """
        Load statistical history from file.

        Parameters
        ----------
        path : Path
            Path to JSON statistics file

        Returns
        -------
        bool
            True if load succeeded
        """
        return self._stats_predictor.load(path)

    def save_statistics(self, path: Path) -> bool:
        """
        Save statistical history to file.

        Parameters
        ----------
        path : Path
            Output path

        Returns
        -------
        bool
            True if save succeeded
        """
        return self._stats_predictor.save(path)

    def get_status(self) -> Dict[str, Any]:
        """
        Get current status of the model manager.

        Returns
        -------
        dict
            Status information including model availability and sample counts
        """
        ml_available = any(
            p.is_available() for p in self._ml_predictors.values()
        )

        return {
            "enabled": self.config["enabled"],
            "ml_available": ml_available,
            "ml_15N_available": self._ml_predictors["15N"].is_available(),
            "ml_13C_available": self._ml_predictors["13C"].is_available(),
            "stats_samples_15N": self._stats_predictor.get_sample_count("15N"),
            "stats_samples_13C": self._stats_predictor.get_sample_count("13C"),
            "use_predictions": self.config["use_predictions"],
            "collect_training_data": self.config["collect_training_data"],
        }

    def get_statistics_summary(self, nucleus: str) -> Optional[Dict]:
        """
        Get summary statistics for a nucleus type.

        Parameters
        ----------
        nucleus : str
            Nucleus type

        Returns
        -------
        dict or None
            Summary statistics if available
        """
        return self._stats_predictor.get_statistics_summary(nucleus)

    def _get_default_params(self, nucleus: str) -> OptimalParameters:
        """Get default parameters for a nucleus type."""
        # Use nucleus-specific defaults
        if nucleus == "15N":
            return OptimalParameters(
                lw_f1_median=0.4,
                lw_f2_median=0.02,
                rad_f1=0.4,
                rad_f2=0.04,
                overlap_threshold_f1=0.4,
                overlap_threshold_f2=0.04,
                achievable_r2=0.85,
                confidence=0.0,
                source="default",
            )
        elif nucleus == "13C":
            return OptimalParameters(
                lw_f1_median=0.15,
                lw_f2_median=0.02,
                rad_f1=0.15,
                rad_f2=0.04,
                overlap_threshold_f1=0.15,
                overlap_threshold_f2=0.04,
                achievable_r2=0.85,
                confidence=0.0,
                source="default",
            )
        else:
            # Generic default
            return OptimalParameters(
                lw_f1_median=0.4,
                lw_f2_median=0.02,
                rad_f1=0.4,
                rad_f2=0.04,
                overlap_threshold_f1=0.4,
                overlap_threshold_f2=0.04,
                achievable_r2=0.85,
                confidence=0.0,
                source="default",
            )

# ABOUTME: Decision engine for choosing between ML and statistical predictions
# ABOUTME: Implements confidence-based selection with divergence detection

from typing import Optional
import logging
import numpy as np

from .feature_extractor import OptimalParameters


logger = logging.getLogger(__name__)


class DecisionEngine:
    """
    Decide between ML and statistical predictions based on confidence.

    Uses confidence thresholds to select the best prediction source.
    Detects when predictions diverge significantly and logs warnings.
    """

    # Default thresholds
    DEFAULT_ML_THRESHOLD = 0.6
    DEFAULT_STATS_THRESHOLD = 0.4
    DEFAULT_DIVERGENCE_THRESHOLD = 0.3  # 30% difference

    def __init__(
        self,
        ml_threshold: float = DEFAULT_ML_THRESHOLD,
        stats_threshold: float = DEFAULT_STATS_THRESHOLD,
        divergence_threshold: float = DEFAULT_DIVERGENCE_THRESHOLD,
    ):
        """
        Initialize the decision engine.

        Parameters
        ----------
        ml_threshold : float
            Minimum confidence to prefer ML prediction
        stats_threshold : float
            Minimum confidence to prefer stats prediction
        divergence_threshold : float
            Relative difference threshold for divergence warning
        """
        self.ml_threshold = ml_threshold
        self.stats_threshold = stats_threshold
        self.divergence_threshold = divergence_threshold

    def decide(
        self,
        ml_pred: Optional[OptimalParameters],
        stats_pred: Optional[OptimalParameters],
    ) -> OptimalParameters:
        """
        Choose between ML and statistical predictions.

        Decision logic:
        1. If ML confidence >= ml_threshold: use ML
        2. Elif stats confidence >= stats_threshold: use stats
        3. Elif ML confidence > stats confidence: use ML (low conf)
        4. Else: use stats (low conf) or default

        Parameters
        ----------
        ml_pred : OptimalParameters or None
            Prediction from ML model
        stats_pred : OptimalParameters or None
            Prediction from statistical model

        Returns
        -------
        OptimalParameters
            Selected prediction with source field updated
        """
        # Handle None cases
        if ml_pred is None and stats_pred is None:
            return self._get_default_params()

        if ml_pred is None:
            return self._select(stats_pred, "stats")

        if stats_pred is None:
            return self._select(ml_pred, "ml")

        # Check for divergence and warn
        if self.predictions_diverge(ml_pred, stats_pred):
            logger.warning(
                "ML and stats predictions diverge significantly. "
                f"ML: lw_f1={ml_pred.lw_f1_median:.3f}, "
                f"Stats: lw_f1={stats_pred.lw_f1_median:.3f}"
            )

        # Decision based on confidence thresholds
        ml_conf = ml_pred.confidence
        stats_conf = stats_pred.confidence

        # Case 1: ML confident
        if ml_conf >= self.ml_threshold:
            return self._select(ml_pred, "ml")

        # Case 2: ML not confident, but stats is
        if stats_conf >= self.stats_threshold:
            return self._select(stats_pred, "stats")

        # Case 3: Neither confident, prefer higher confidence
        if ml_conf > stats_conf:
            return self._select(ml_pred, "ml_low_conf")

        # Case 4: Stats has equal or higher confidence
        return self._select(stats_pred, "stats_low_conf")

    def predictions_diverge(
        self,
        ml_pred: OptimalParameters,
        stats_pred: OptimalParameters,
    ) -> bool:
        """
        Check if predictions diverge significantly.

        Parameters
        ----------
        ml_pred : OptimalParameters
            ML prediction
        stats_pred : OptimalParameters
            Stats prediction

        Returns
        -------
        bool
            True if predictions differ by more than threshold
        """
        # Compare key parameters
        params_to_check = [
            (ml_pred.lw_f1_median, stats_pred.lw_f1_median),
            (ml_pred.lw_f2_median, stats_pred.lw_f2_median),
            (ml_pred.rad_f1, stats_pred.rad_f1),
            (ml_pred.rad_f2, stats_pred.rad_f2),
        ]

        for ml_val, stats_val in params_to_check:
            # Calculate relative difference
            avg = (abs(ml_val) + abs(stats_val)) / 2
            if avg > 0:
                rel_diff = abs(ml_val - stats_val) / avg
                if rel_diff > self.divergence_threshold:
                    return True

        return False

    def _select(
        self,
        pred: OptimalParameters,
        source: str,
    ) -> OptimalParameters:
        """Create a copy of prediction with updated source."""
        return OptimalParameters(
            lw_f1_median=pred.lw_f1_median,
            lw_f2_median=pred.lw_f2_median,
            rad_f1=pred.rad_f1,
            rad_f2=pred.rad_f2,
            overlap_threshold_f1=pred.overlap_threshold_f1,
            overlap_threshold_f2=pred.overlap_threshold_f2,
            achievable_r2=pred.achievable_r2,
            confidence=pred.confidence,
            source=source,
        )

    def _get_default_params(self) -> OptimalParameters:
        """Return default parameters when no predictions available."""
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

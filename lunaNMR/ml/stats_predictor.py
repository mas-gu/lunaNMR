# ABOUTME: Statistical fallback predictor using histogram-based estimates
# ABOUTME: Provides predictions when ML model is unavailable or low confidence

from pathlib import Path
from typing import Dict, List, Optional
from dataclasses import dataclass, field
from datetime import datetime
import json
import logging
import numpy as np

from .feature_extractor import SpectrumFeatures, OptimalParameters


logger = logging.getLogger(__name__)


@dataclass
class NucleusStatistics:
    """Statistics for a single nucleus type."""
    lw_f1_values: List[float] = field(default_factory=list)
    lw_f2_values: List[float] = field(default_factory=list)
    rad_f1_values: List[float] = field(default_factory=list)
    rad_f2_values: List[float] = field(default_factory=list)
    overlap_f1_values: List[float] = field(default_factory=list)
    overlap_f2_values: List[float] = field(default_factory=list)
    r2_values: List[float] = field(default_factory=list)
    last_updated: str = ""

    def to_dict(self) -> dict:
        """Convert to dictionary for JSON serialization."""
        return {
            "lw_f1": self.lw_f1_values,
            "lw_f2": self.lw_f2_values,
            "rad_f1": self.rad_f1_values,
            "rad_f2": self.rad_f2_values,
            "overlap_f1": self.overlap_f1_values,
            "overlap_f2": self.overlap_f2_values,
            "r2": self.r2_values,
            "last_updated": self.last_updated,
            "sample_count": len(self.lw_f1_values),
        }

    @classmethod
    def from_dict(cls, data: dict) -> "NucleusStatistics":
        """Create from dictionary."""
        return cls(
            lw_f1_values=data.get("lw_f1", []),
            lw_f2_values=data.get("lw_f2", []),
            rad_f1_values=data.get("rad_f1", []),
            rad_f2_values=data.get("rad_f2", []),
            overlap_f1_values=data.get("overlap_f1", []),
            overlap_f2_values=data.get("overlap_f2", []),
            r2_values=data.get("r2", []),
            last_updated=data.get("last_updated", ""),
        )


class StatisticalPredictor:
    """
    Statistical fallback predictor using histogram-based estimates.

    Maintains running statistics per nucleus type and predicts optimal
    parameters using median values. Provides confidence based on sample
    count and variance.
    """

    # Minimum samples for reasonable prediction
    MIN_SAMPLES = 5

    # Target sample count for max confidence
    MAX_CONFIDENCE_SAMPLES = 50

    def __init__(self):
        """Initialize the statistical predictor."""
        self.statistics: Dict[str, NucleusStatistics] = {
            "15N": NucleusStatistics(),
            "13C": NucleusStatistics(),
        }

    def update(self, nucleus: str, params: OptimalParameters) -> None:
        """
        Add a successful parameter set to the statistics.

        Parameters
        ----------
        nucleus : str
            Nucleus type ("15N" or "13C")
        params : OptimalParameters
            Parameters from successful processing
        """
        if nucleus not in self.statistics:
            self.statistics[nucleus] = NucleusStatistics()

        stats = self.statistics[nucleus]
        stats.lw_f1_values.append(params.lw_f1_median)
        stats.lw_f2_values.append(params.lw_f2_median)
        stats.rad_f1_values.append(params.rad_f1)
        stats.rad_f2_values.append(params.rad_f2)
        stats.overlap_f1_values.append(params.overlap_threshold_f1)
        stats.overlap_f2_values.append(params.overlap_threshold_f2)
        stats.r2_values.append(params.achievable_r2)
        stats.last_updated = datetime.now().isoformat()

    def predict(self, features: SpectrumFeatures) -> OptimalParameters:
        """
        Predict parameters using statistical estimates.

        Parameters
        ----------
        features : SpectrumFeatures
            Features from the spectrum (nucleus_type used for selection)

        Returns
        -------
        OptimalParameters
            Predicted parameters with confidence
        """
        nucleus = features.nucleus_type

        if nucleus not in self.statistics:
            return self._get_default_params()

        stats = self.statistics[nucleus]
        n_samples = len(stats.lw_f1_values)

        if n_samples < self.MIN_SAMPLES:
            return self._get_default_params()

        # Calculate medians (robust to outliers)
        lw_f1_median = float(np.median(stats.lw_f1_values))
        lw_f2_median = float(np.median(stats.lw_f2_values))
        rad_f1 = float(np.median(stats.rad_f1_values))
        rad_f2 = float(np.median(stats.rad_f2_values))
        overlap_f1 = float(np.median(stats.overlap_f1_values))
        overlap_f2 = float(np.median(stats.overlap_f2_values))
        achievable_r2 = float(np.median(stats.r2_values))

        # Calculate confidence based on sample count
        confidence = min(1.0, n_samples / self.MAX_CONFIDENCE_SAMPLES)

        return OptimalParameters(
            lw_f1_median=lw_f1_median,
            lw_f2_median=lw_f2_median,
            rad_f1=rad_f1,
            rad_f2=rad_f2,
            overlap_threshold_f1=overlap_f1,
            overlap_threshold_f2=overlap_f2,
            achievable_r2=achievable_r2,
            confidence=confidence,
            source="stats",
        )

    def get_sample_count(self, nucleus: str) -> int:
        """Get number of samples collected for a nucleus type."""
        if nucleus not in self.statistics:
            return 0
        return len(self.statistics[nucleus].lw_f1_values)

    def save(self, path: Path) -> bool:
        """
        Save statistics to JSON file.

        Parameters
        ----------
        path : Path
            Output file path

        Returns
        -------
        bool
            True if save succeeded
        """
        try:
            path = Path(path)
            path.parent.mkdir(parents=True, exist_ok=True)

            data = {
                nucleus: stats.to_dict()
                for nucleus, stats in self.statistics.items()
            }

            with open(path, 'w') as f:
                json.dump(data, f, indent=2)

            logger.info(f"Statistics saved to {path}")
            return True

        except Exception as e:
            logger.error(f"Failed to save statistics: {e}")
            return False

    def load(self, path: Path) -> bool:
        """
        Load statistics from JSON file.

        Parameters
        ----------
        path : Path
            Input file path

        Returns
        -------
        bool
            True if load succeeded
        """
        try:
            path = Path(path)
            if not path.exists():
                logger.warning(f"Statistics file not found: {path}")
                return False

            with open(path, 'r') as f:
                data = json.load(f)

            for nucleus, stats_dict in data.items():
                self.statistics[nucleus] = NucleusStatistics.from_dict(stats_dict)

            logger.info(f"Statistics loaded from {path}")
            return True

        except Exception as e:
            logger.error(f"Failed to load statistics: {e}")
            return False

    def _get_default_params(self) -> OptimalParameters:
        """Return default parameters when insufficient data."""
        return OptimalParameters(
            lw_f1_median=0.4,
            lw_f2_median=0.02,
            rad_f1=0.4,
            rad_f2=0.04,
            overlap_threshold_f1=0.4,
            overlap_threshold_f2=0.04,
            achievable_r2=0.85,
            confidence=0.0,
            source="stats_default",
        )

    def get_statistics_summary(self, nucleus: str) -> Optional[dict]:
        """
        Get summary statistics for a nucleus type.

        Returns
        -------
        dict or None
            Summary with median, std, min, max for each parameter
        """
        if nucleus not in self.statistics:
            return None

        stats = self.statistics[nucleus]
        if len(stats.lw_f1_values) < self.MIN_SAMPLES:
            return None

        def summarize(values: List[float]) -> dict:
            arr = np.array(values)
            return {
                "median": float(np.median(arr)),
                "std": float(np.std(arr)),
                "min": float(np.min(arr)),
                "max": float(np.max(arr)),
                "n": len(values),
            }

        return {
            "lw_f1": summarize(stats.lw_f1_values),
            "lw_f2": summarize(stats.lw_f2_values),
            "rad_f1": summarize(stats.rad_f1_values),
            "rad_f2": summarize(stats.rad_f2_values),
            "overlap_f1": summarize(stats.overlap_f1_values),
            "overlap_f2": summarize(stats.overlap_f2_values),
            "r2": summarize(stats.r2_values),
            "sample_count": len(stats.lw_f1_values),
            "last_updated": stats.last_updated,
        }

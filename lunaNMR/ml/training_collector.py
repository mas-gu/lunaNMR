# ABOUTME: Training data collector for ML model adaptation
# ABOUTME: Collects successful fit results for local model retraining

from pathlib import Path
from typing import Dict, List, Optional, Tuple
from dataclasses import dataclass, field, asdict
from datetime import datetime
import json
import logging

from .feature_extractor import SpectrumFeatures, OptimalParameters


logger = logging.getLogger(__name__)


@dataclass
class TrainingSample:
    """A single training sample with features and targets."""
    features: Dict  # Serialized SpectrumFeatures
    targets: Dict   # Serialized OptimalParameters
    median_r2: float
    timestamp: str
    nucleus: str


class TrainingCollector:
    """
    Collects training data from successful spectrum processing.

    Stores spectrum features and optimal parameters for model adaptation.
    Filters by quality threshold and manages sample limits.
    """

    # Default storage location
    DEFAULT_STORAGE_NAME = "local_training_data.json"

    # Default limits
    DEFAULT_MIN_R2 = 0.80
    DEFAULT_MAX_SAMPLES = 500  # Per nucleus type

    def __init__(
        self,
        storage_path: Optional[Path] = None,
        min_r2: float = DEFAULT_MIN_R2,
        max_samples_per_nucleus: int = DEFAULT_MAX_SAMPLES,
    ):
        """
        Initialize the training collector.

        Parameters
        ----------
        storage_path : Path, optional
            Path to store training data. If None, uses default location.
        min_r2 : float
            Minimum R² to accept a sample
        max_samples_per_nucleus : int
            Maximum samples to keep per nucleus type (FIFO)
        """
        self.storage_path = storage_path
        self.min_r2 = min_r2
        self.max_samples_per_nucleus = max_samples_per_nucleus

        # Storage by nucleus type
        self._samples: Dict[str, List[TrainingSample]] = {
            "15N": [],
            "13C": [],
        }

    def collect(
        self,
        features: SpectrumFeatures,
        params: OptimalParameters,
        median_r2: float,
    ) -> bool:
        """
        Collect a training sample from successful processing.

        Parameters
        ----------
        features : SpectrumFeatures
            Spectrum features
        params : OptimalParameters
            Optimal parameters from successful processing
        median_r2 : float
            Median R² achieved across peaks

        Returns
        -------
        bool
            True if sample was accepted, False if rejected
        """
        # Quality filter
        if median_r2 < self.min_r2:
            logger.debug(f"Rejected sample: R² {median_r2:.3f} < {self.min_r2}")
            return False

        nucleus = features.nucleus_type
        if nucleus not in self._samples:
            self._samples[nucleus] = []

        # Create sample
        sample = TrainingSample(
            features=self._serialize_features(features),
            targets=self._serialize_params(params),
            median_r2=median_r2,
            timestamp=datetime.now().isoformat(),
            nucleus=nucleus,
        )

        # Add to collection
        self._samples[nucleus].append(sample)

        # Enforce max samples (FIFO - drop oldest)
        if len(self._samples[nucleus]) > self.max_samples_per_nucleus:
            self._samples[nucleus] = self._samples[nucleus][-self.max_samples_per_nucleus:]

        logger.debug(f"Collected {nucleus} sample: R² {median_r2:.3f}")
        return True

    def get_sample_count(self, nucleus: str) -> int:
        """Get number of samples for a nucleus type."""
        if nucleus not in self._samples:
            return 0
        return len(self._samples[nucleus])

    def get_training_data(
        self,
        nucleus: str,
    ) -> Tuple[List[SpectrumFeatures], List[OptimalParameters]]:
        """
        Get training data for model training.

        Parameters
        ----------
        nucleus : str
            Nucleus type ("15N" or "13C")

        Returns
        -------
        tuple
            (features_list, targets_list) for training
        """
        if nucleus not in self._samples:
            return [], []

        features_list = []
        targets_list = []

        for sample in self._samples[nucleus]:
            features = self._deserialize_features(sample.features)
            targets = self._deserialize_params(sample.targets)
            features_list.append(features)
            targets_list.append(targets)

        return features_list, targets_list

    def save(self) -> bool:
        """
        Save collected samples to storage.

        Returns
        -------
        bool
            True if save succeeded
        """
        if self.storage_path is None:
            logger.warning("No storage path configured")
            return False

        try:
            path = Path(self.storage_path)
            path.parent.mkdir(parents=True, exist_ok=True)

            # Convert to serializable format
            data = {
                "version": "1.0.0",
                "updated": datetime.now().isoformat(),
                "samples": {}
            }

            for nucleus, samples in self._samples.items():
                data["samples"][nucleus] = [
                    {
                        "features": s.features,
                        "targets": s.targets,
                        "median_r2": s.median_r2,
                        "timestamp": s.timestamp,
                        "nucleus": s.nucleus,
                    }
                    for s in samples
                ]

            with open(path, 'w') as f:
                json.dump(data, f, indent=2)

            logger.info(f"Saved {sum(len(s) for s in self._samples.values())} training samples to {path}")
            return True

        except Exception as e:
            logger.error(f"Failed to save training data: {e}")
            return False

    def load(self) -> bool:
        """
        Load samples from storage.

        Returns
        -------
        bool
            True if load succeeded
        """
        if self.storage_path is None:
            return False

        try:
            path = Path(self.storage_path)
            if not path.exists():
                return False

            with open(path, 'r') as f:
                data = json.load(f)

            # Clear existing
            self._samples = {"15N": [], "13C": []}

            # Load samples
            for nucleus, samples in data.get("samples", {}).items():
                if nucleus not in self._samples:
                    self._samples[nucleus] = []

                for s in samples:
                    sample = TrainingSample(
                        features=s["features"],
                        targets=s["targets"],
                        median_r2=s["median_r2"],
                        timestamp=s["timestamp"],
                        nucleus=s["nucleus"],
                    )
                    self._samples[nucleus].append(sample)

            logger.info(f"Loaded {sum(len(s) for s in self._samples.values())} training samples from {path}")
            return True

        except Exception as e:
            logger.error(f"Failed to load training data: {e}")
            return False

    def _serialize_features(self, features: SpectrumFeatures) -> Dict:
        """Serialize SpectrumFeatures to dict."""
        return {
            "nucleus_type": features.nucleus_type,
            "snr_estimate": features.snr_estimate,
            "noise_level": features.noise_level,
            "dynamic_range": features.dynamic_range,
            "peak_count": features.peak_count,
            "peak_density": features.peak_density,
            "shift_range_f1_min": features.shift_range_f1_min,
            "shift_range_f1_max": features.shift_range_f1_max,
            "shift_range_f2_min": features.shift_range_f2_min,
            "shift_range_f2_max": features.shift_range_f2_max,
            "field_strength_mhz": features.field_strength_mhz,
        }

    def _deserialize_features(self, data: Dict) -> SpectrumFeatures:
        """Deserialize dict to SpectrumFeatures."""
        return SpectrumFeatures(
            nucleus_type=data["nucleus_type"],
            snr_estimate=data["snr_estimate"],
            noise_level=data["noise_level"],
            dynamic_range=data["dynamic_range"],
            peak_count=data["peak_count"],
            peak_density=data["peak_density"],
            shift_range_f1_min=data["shift_range_f1_min"],
            shift_range_f1_max=data["shift_range_f1_max"],
            shift_range_f2_min=data["shift_range_f2_min"],
            shift_range_f2_max=data["shift_range_f2_max"],
            field_strength_mhz=data["field_strength_mhz"],
        )

    def _serialize_params(self, params: OptimalParameters) -> Dict:
        """Serialize OptimalParameters to dict."""
        return {
            "lw_f1_median": params.lw_f1_median,
            "lw_f2_median": params.lw_f2_median,
            "rad_f1": params.rad_f1,
            "rad_f2": params.rad_f2,
            "overlap_threshold_f1": params.overlap_threshold_f1,
            "overlap_threshold_f2": params.overlap_threshold_f2,
            "achievable_r2": params.achievable_r2,
        }

    def _deserialize_params(self, data: Dict) -> OptimalParameters:
        """Deserialize dict to OptimalParameters."""
        return OptimalParameters(
            lw_f1_median=data["lw_f1_median"],
            lw_f2_median=data["lw_f2_median"],
            rad_f1=data["rad_f1"],
            rad_f2=data["rad_f2"],
            overlap_threshold_f1=data["overlap_threshold_f1"],
            overlap_threshold_f2=data["overlap_threshold_f2"],
            achievable_r2=data["achievable_r2"],
        )

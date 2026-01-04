# ABOUTME: Feature extraction for ML parameter prediction
# ABOUTME: Extracts spectrum characteristics and optimal parameters from fit results

from dataclasses import dataclass, field
from typing import List, Dict, Any, Optional, Tuple
import numpy as np


@dataclass
class SpectrumFeatures:
    """Features extracted from a spectrum for ML prediction."""

    # Categorical
    nucleus_type: str  # "15N" or "13C"

    # Signal quality
    snr_estimate: float
    noise_level: float
    dynamic_range: float

    # Spectrum characteristics
    peak_count: int
    peak_density: float  # peaks per ppm² area

    # Chemical shift ranges
    shift_range_f1_min: float
    shift_range_f1_max: float
    shift_range_f2_min: float
    shift_range_f2_max: float

    # Optional metadata
    field_strength_mhz: float = 600.0

    def to_array(self) -> np.ndarray:
        """Convert features to numpy array for ML model input."""
        # Encode nucleus type as binary (15N=1, 13C=0)
        nucleus_encoded = 1.0 if self.nucleus_type == "15N" else 0.0

        return np.array([
            nucleus_encoded,
            self.snr_estimate,
            self.noise_level,
            self.dynamic_range,
            float(self.peak_count),
            self.peak_density,
            self.shift_range_f1_min,
            self.shift_range_f1_max,
            self.shift_range_f2_min,
            self.shift_range_f2_max,
            self.field_strength_mhz,
        ], dtype=np.float64)

    @staticmethod
    def field_names() -> List[str]:
        """Return list of feature field names for ML model."""
        return [
            "nucleus_encoded",
            "snr_estimate",
            "noise_level",
            "dynamic_range",
            "peak_count",
            "peak_density",
            "shift_range_f1_min",
            "shift_range_f1_max",
            "shift_range_f2_min",
            "shift_range_f2_max",
            "field_strength_mhz",
        ]


@dataclass
class OptimalParameters:
    """Optimal parameters learned from successful fits."""

    # Linewidths (ppm)
    lw_f1_median: float
    lw_f2_median: float

    # Integration window radii (ppm)
    rad_f1: float
    rad_f2: float

    # Overlap detection thresholds (ppm)
    overlap_threshold_f1: float
    overlap_threshold_f2: float

    # Quality expectation
    achievable_r2: float

    # Metadata (not used in training)
    confidence: float = 0.0
    source: str = "default"

    def to_array(self) -> np.ndarray:
        """Convert parameters to numpy array for ML training targets."""
        return np.array([
            self.lw_f1_median,
            self.lw_f2_median,
            self.rad_f1,
            self.rad_f2,
            self.overlap_threshold_f1,
            self.overlap_threshold_f2,
            self.achievable_r2,
        ], dtype=np.float64)

    @staticmethod
    def field_names() -> List[str]:
        """Return list of target parameter names."""
        return [
            "lw_f1_median",
            "lw_f2_median",
            "rad_f1",
            "rad_f2",
            "overlap_threshold_f1",
            "overlap_threshold_f2",
            "achievable_r2",
        ]

    @classmethod
    def from_array(cls, arr: np.ndarray, confidence: float = 0.0,
                   source: str = "default") -> "OptimalParameters":
        """Create OptimalParameters from numpy array."""
        return cls(
            lw_f1_median=float(arr[0]),
            lw_f2_median=float(arr[1]),
            rad_f1=float(arr[2]),
            rad_f2=float(arr[3]),
            overlap_threshold_f1=float(arr[4]),
            overlap_threshold_f2=float(arr[5]),
            achievable_r2=float(arr[6]),
            confidence=confidence,
            source=source,
        )


class FeatureExtractor:
    """Extract features from spectra and fit results for ML."""

    # Nucleus detection ranges (ppm)
    NUCLEUS_RANGES = {
        "15N": (90.0, 150.0),   # 15N typically 100-140 ppm
        "13C": (0.0, 230.0),    # 13C wide range
        "1H": (0.0, 15.0),      # 1H proton range
    }

    def __init__(self):
        """Initialize the feature extractor."""
        pass

    def extract_from_spectrum(
        self,
        spectrum_data: np.ndarray,
        ppm_f1: np.ndarray,
        ppm_f2: np.ndarray,
        peak_count: int,
        field_strength_mhz: float = 600.0
    ) -> SpectrumFeatures:
        """
        Extract features from raw spectrum data.

        Parameters
        ----------
        spectrum_data : np.ndarray
            2D spectrum intensity array
        ppm_f1 : np.ndarray
            F1 dimension ppm axis (15N or 13C)
        ppm_f2 : np.ndarray
            F2 dimension ppm axis (1H)
        peak_count : int
            Number of detected peaks
        field_strength_mhz : float
            Spectrometer field strength in MHz

        Returns
        -------
        SpectrumFeatures
            Extracted features for ML prediction
        """
        # Detect nucleus type from F1 range
        nucleus_type = self._detect_nucleus_type(ppm_f1)

        # Calculate signal quality metrics
        noise_level = self._estimate_noise(spectrum_data)
        signal_level = np.percentile(np.abs(spectrum_data), 95)
        snr_estimate = signal_level / noise_level if noise_level > 0 else 0.0
        dynamic_range = float(np.max(spectrum_data) - np.min(spectrum_data))

        # Calculate spectrum area and peak density
        f1_range = float(np.max(ppm_f1) - np.min(ppm_f1))
        f2_range = float(np.max(ppm_f2) - np.min(ppm_f2))
        spectrum_area = f1_range * f2_range
        peak_density = peak_count / spectrum_area if spectrum_area > 0 else 0.0

        return SpectrumFeatures(
            nucleus_type=nucleus_type,
            snr_estimate=float(snr_estimate),
            noise_level=float(noise_level),
            dynamic_range=dynamic_range,
            peak_count=peak_count,
            peak_density=peak_density,
            shift_range_f1_min=float(np.min(ppm_f1)),
            shift_range_f1_max=float(np.max(ppm_f1)),
            shift_range_f2_min=float(np.min(ppm_f2)),
            shift_range_f2_max=float(np.max(ppm_f2)),
            field_strength_mhz=field_strength_mhz,
        )

    def extract_optimal_params(
        self,
        fit_results: List[Dict[str, Any]],
        min_r2: float = 0.80
    ) -> OptimalParameters:
        """
        Extract optimal parameters from successful fit results.

        Parameters
        ----------
        fit_results : list
            List of fit result dictionaries from processing
        min_r2 : float
            Minimum R² threshold for including a fit

        Returns
        -------
        OptimalParameters
            Optimal parameters derived from successful fits
        """
        # Filter to successful, high-quality fits
        good_fits = [
            f for f in fit_results
            if f.get("success", True) and f.get("r_squared", 0) >= min_r2
        ]

        if not good_fits:
            # Return defaults if no good fits
            return self._get_default_params()

        # Extract linewidths from good fits
        lw_f1_values = []
        lw_f2_values = []
        r2_values = []

        for fit in good_fits:
            if "lw_f1" in fit:
                lw_f1_values.append(fit["lw_f1"])
            if "lw_f2" in fit:
                lw_f2_values.append(fit["lw_f2"])
            if "r_squared" in fit:
                r2_values.append(fit["r_squared"])

        # Calculate medians (robust to outliers)
        lw_f1_median = float(np.median(lw_f1_values)) if lw_f1_values else 0.4
        lw_f2_median = float(np.median(lw_f2_values)) if lw_f2_values else 0.02
        achievable_r2 = float(np.median(r2_values)) if r2_values else 0.85

        # Integration windows typically 1-2x linewidth
        rad_f1 = lw_f1_median * 1.2
        rad_f2 = lw_f2_median * 1.2

        # Overlap thresholds similar to integration windows
        overlap_threshold_f1 = lw_f1_median
        overlap_threshold_f2 = lw_f2_median

        return OptimalParameters(
            lw_f1_median=lw_f1_median,
            lw_f2_median=lw_f2_median,
            rad_f1=rad_f1,
            rad_f2=rad_f2,
            overlap_threshold_f1=overlap_threshold_f1,
            overlap_threshold_f2=overlap_threshold_f2,
            achievable_r2=achievable_r2,
            confidence=min(1.0, len(good_fits) / 50.0),  # Confidence from sample size
            source="extracted",
        )

    def _detect_nucleus_type(self, ppm_f1: np.ndarray) -> str:
        """
        Detect nucleus type from F1 chemical shift range.

        Uses range heuristics:
        - 15N: 90-150 ppm (typically 100-140)
        - 13C: 0-230 ppm (wide range, >50 ppm span)
        - 1H: 0-15 ppm (narrow range)
        """
        ppm_min = float(np.min(ppm_f1))
        ppm_max = float(np.max(ppm_f1))
        ppm_center = (ppm_min + ppm_max) / 2
        ppm_span = ppm_max - ppm_min

        # 15N detection: center around 100-140, span < 100
        if 90 <= ppm_center <= 150 and ppm_span < 100:
            return "15N"

        # 13C detection: wide range (>50 ppm span) or center in 0-200
        if ppm_span > 50 or (0 <= ppm_center <= 200 and ppm_span > 20):
            return "13C"

        # Default to 1H for narrow ranges
        return "1H"

    def _estimate_noise(self, spectrum_data: np.ndarray) -> float:
        """
        Estimate noise level from spectrum corners.

        Uses median absolute deviation (MAD) of corner regions
        for robust noise estimation.
        """
        height, width = spectrum_data.shape
        corner_size = max(10, min(height // 10, width // 10))

        # Extract corner regions
        corners = [
            spectrum_data[:corner_size, :corner_size],
            spectrum_data[:corner_size, -corner_size:],
            spectrum_data[-corner_size:, :corner_size],
            spectrum_data[-corner_size:, -corner_size:],
        ]

        # Combine corners and compute MAD
        corner_data = np.concatenate([c.flatten() for c in corners])
        median = np.median(corner_data)
        mad = np.median(np.abs(corner_data - median))

        # Scale MAD to standard deviation equivalent
        noise_estimate = mad * 1.4826

        return float(noise_estimate) if noise_estimate > 0 else 1.0

    def _get_default_params(self) -> OptimalParameters:
        """Return default parameters when no good fits available."""
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

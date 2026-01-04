# ABOUTME: Comprehensive ML training data collector with ~155 parameters
# ABOUTME: Captures spectrum features, fitting results, adaptive optimization, and improvements

import json
import logging
import numpy as np
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Any
from dataclasses import dataclass, field, asdict
from scipy import stats as scipy_stats

from .storage import get_training_data_path, get_session_log_dir, get_metadata_path

logger = logging.getLogger(__name__)


@dataclass
class ComprehensiveTrainingSample:
    """Complete training sample with all ~143 parameters."""

    # === Processing Metadata ===
    timestamp: str = ""
    spectrum_name: str = ""
    software_version: str = "1.0.0"
    processing_mode: str = "sequential"  # sequential or parallel
    is_reference_spectrum: bool = False
    series_index: int = 0

    # === Spectrum Features ===
    nucleus_type: str = "15N"
    field_strength_mhz: float = 600.0
    snr_estimate: float = 0.0
    noise_level: float = 0.0
    dynamic_range: float = 0.0
    peak_count: int = 0
    peak_density: float = 0.0
    shift_range_f1_min: float = 0.0
    shift_range_f1_max: float = 0.0
    shift_range_f2_min: float = 0.0
    shift_range_f2_max: float = 0.0

    # === Spectrum Quality ===
    # Volume statistics (fitted integral for normalized Voigt)
    volume_mean: float = 0.0
    volume_std: float = 0.0
    volume_skewness: float = 0.0
    volume_kurtosis: float = 0.0
    volume_cv: float = 0.0
    volume_dynamic_range: float = 0.0
    # Height statistics (peak maximum at center, different from volume!)
    height_mean: float = 0.0
    height_std: float = 0.0
    height_skewness: float = 0.0
    height_kurtosis: float = 0.0
    height_cv: float = 0.0
    height_dynamic_range: float = 0.0
    # Detected intensity statistics (from peak detection before fitting)
    detected_intensity_mean: float = 0.0
    detected_intensity_std: float = 0.0
    detected_intensity_skewness: float = 0.0
    detected_intensity_kurtosis: float = 0.0
    detected_intensity_cv: float = 0.0
    detected_intensity_dynamic_range: float = 0.0
    # Baseline
    baseline_level: float = 0.0
    baseline_std: float = 0.0

    # === Peak Distribution & Clustering ===
    n_isolated_peaks: int = 0
    n_clustered_peaks: int = 0
    n_clusters: int = 0
    cluster_sizes: List[int] = field(default_factory=list)
    max_cluster_size: int = 0
    mean_cluster_size: float = 0.0
    fraction_clustered: float = 0.0
    mean_peak_separation_f1: float = 0.0
    mean_peak_separation_f2: float = 0.0
    min_peak_separation_f1: float = 0.0
    min_peak_separation_f2: float = 0.0

    # === Protein Characteristics ===
    f1_dispersion: float = 0.0
    f2_dispersion: float = 0.0
    is_idp_like: bool = False
    estimated_protein_size: str = "medium"
    folding_indicator: float = 0.0
    dispersion_ratio: float = 0.0  # f1_dispersion / f2_dispersion

    # === Pre-fitting Clustering Estimates (computed from peak positions only) ===
    n_close_pairs_f1: int = 0           # Pairs within 2*typical_lw in F1
    n_close_pairs_f2: int = 0           # Pairs within 2*typical_lw in F2
    n_close_pairs_2d: int = 0           # Pairs within elliptical distance
    fraction_potentially_overlapping: float = 0.0  # Fraction of peaks near others

    # === Peak List Intensity Stats (from peak list before fitting) ===
    peaklist_intensity_mean: float = 0.0
    peaklist_intensity_std: float = 0.0
    peaklist_intensity_cv: float = 0.0
    peaklist_intensity_dynamic_range: float = 0.0

    # === Crowding Metrics ===
    max_local_density: float = 0.0      # Peaks per ppm² in densest region
    crowding_hotspot_fraction: float = 0.0  # Fraction of peaks in densest 10%

    # === Peak Detection Parameters ===
    detection_search_window_x: float = 0.01   # 1H search window (ppm)
    detection_search_window_y: float = 0.04   # 15N search window (ppm)
    detection_noise_threshold: float = 3.0    # S/N threshold multiplier
    detection_noise_level: float = 0.0        # Estimated noise level

    # === Reference vs Detected Statistics ===
    reference_peak_count: int = 0              # Original peak list size
    detected_peak_count: int = 0               # Peaks successfully matched
    reference_retained_count: int = 0          # Peaks using reference position
    dummy_peaks_added: int = 0                 # Auto-added nearby peaks
    detection_rate: float = 0.0                # Percentage matched

    # Position accuracy statistics
    mean_distance_from_reference: float = 0.0  # Mean shift from reference (ppm)
    max_distance_from_reference: float = 0.0   # Max shift from reference (ppm)
    std_distance_from_reference: float = 0.0   # Std of shifts

    # === PASS1 Learned Statistics ===
    pass1_lw_f1_median: float = 0.0
    pass1_lw_f2_median: float = 0.0
    pass1_lw_f1_std: float = 0.0
    pass1_lw_f2_std: float = 0.0
    pass1_lw_f1_min: float = 0.0
    pass1_lw_f1_max: float = 0.0
    pass1_lw_f2_min: float = 0.0
    pass1_lw_f2_max: float = 0.0
    pass1_n_good_peaks: int = 0
    pass1_mean_r_squared: float = 0.0

    # === Adaptive Optimization Results ===
    adaptive_success: bool = False
    adaptive_fallback_reason: str = ""
    adaptive_radF1: float = 0.0
    adaptive_radF2: float = 0.0
    adaptive_overlap_threshold_x: float = 0.0
    adaptive_overlap_threshold_y: float = 0.0
    adaptive_multiplier_f1: float = 1.5
    adaptive_multiplier_f2: float = 1.5
    adaptive_train_score: float = 0.0
    adaptive_validation_score: float = 0.0
    adaptive_generalization_gap: float = 0.0
    adaptive_n_peaks_used: int = 0
    adaptive_n_train: int = 0
    adaptive_n_val: int = 0
    adaptive_search_history: List[Dict] = field(default_factory=list)

    # === PASS1 -> PASS1-bis Improvement ===
    pass1bis_enabled: bool = False
    pass1bis_n_peaks_refit: int = 0
    pass1bis_mean_r2_before: float = 0.0
    pass1bis_mean_r2_after: float = 0.0
    pass1bis_r2_improvement: float = 0.0
    pass1bis_n_improved: int = 0
    pass1bis_n_degraded: int = 0
    pass1bis_n_unchanged: int = 0

    # === Re-clustering Results ===
    recluster_performed: bool = False
    recluster_n_originally_clustered: int = 0
    recluster_n_now_isolated: int = 0
    recluster_n_still_clustered: int = 0
    recluster_new_cluster_count: int = 0
    recluster_cluster_size_change: float = 0.0

    # === PASS2 Results ===
    pass2_mean_r_squared: float = 0.0
    pass2_n_peaks_fitted: int = 0
    pass2_n_clusters_fitted: int = 0
    pass2_total_time: float = 0.0

    # === Overall Improvement ===
    overall_r2_pass1: float = 0.0
    overall_r2_pass1bis: float = 0.0
    overall_r2_pass2: float = 0.0
    improvement_pass1_to_pass1bis: float = 0.0
    improvement_pass1bis_to_pass2: float = 0.0
    total_improvement: float = 0.0

    # === Timing ===
    total_processing_time: float = 0.0
    pass1_time: float = 0.0
    adaptive_time: float = 0.0
    pass1bis_time: float = 0.0
    pass2_time: float = 0.0

    # === Per-Peak Data (stored separately) ===
    peaks: List[Dict] = field(default_factory=list)

    def to_dict(self) -> Dict:
        """Convert to dictionary for JSON serialization."""
        return asdict(self)

    @classmethod
    def from_dict(cls, data: Dict) -> 'ComprehensiveTrainingSample':
        """Create from dictionary."""
        return cls(**data)


@dataclass
class PeakTrainingData:
    """Per-peak training data."""

    # Identification
    peak_index: int = 0
    peak_id: str = ""
    assignment: str = ""  # Peak assignment (e.g., "G45N-H", "A123CA")

    # Detection results (reference list → detected peaks)
    was_detected: bool = True              # True if matched to detected peak
    distance_from_reference: float = 0.0   # Euclidean distance (ppm) - legacy
    distance_from_reference_x: float = 0.0 # X/F2 dimension distance (1H, ppm)
    distance_from_reference_y: float = 0.0 # Y/F1 dimension distance (15N, ppm)
    distance_from_reference_elliptical: float = 0.0  # Normalized elliptical distance
    detection_quality: str = "Matched"     # "Matched", "Reference", "Auto-added"

    # Fitted parameters
    pos_f1: float = 0.0
    pos_f2: float = 0.0
    lw_lor_f1: float = 0.0
    lw_gau_f1: float = 0.0
    lw_lor_f2: float = 0.0
    lw_gau_f2: float = 0.0
    lw_total_f1: float = 0.0
    lw_total_f2: float = 0.0
    lg_ratio_f1: float = 0.0
    lg_ratio_f2: float = 0.0
    volume: float = 0.0     # Fitted integral (normalized Voigt parameter)
    height: float = 0.0     # Peak maximum at center (different from volume!)
    detected_intensity: float = 0.0

    # Refinement metrics
    position_shift_f1: float = 0.0
    position_shift_f2: float = 0.0
    linewidth_change_f1: float = 0.0
    linewidth_change_f2: float = 0.0
    initial_pos_f1: float = 0.0
    initial_pos_f2: float = 0.0
    initial_lw_f1: float = 0.0
    initial_lw_f2: float = 0.0

    # Context flags
    is_isolated: bool = True
    is_clustered: bool = False
    cluster_id: int = -1
    cluster_size: int = 1
    tooclose_flag: bool = False
    heavy_overlap_flag: bool = False
    was_refit_in_pass1bis: bool = False

    # Uncertainties
    intensity_uncertainty: float = 0.0

    # Fit quality
    r_squared: float = 0.0
    chi2: float = 0.0
    iterations: int = 0
    convergence_type: str = ""  # "formal", "pragmatic_r2", "chi2_reduction", or ""
    success: bool = False

    # Fitting strategy
    fix_positions_applied: bool = False
    fix_linewidths_applied: bool = False
    fitting_mode: str = "1D"  # "1D" or "2D"
    n_peaks_in_fit: int = 1
    used_lg_penalty: bool = False
    used_intensity_ratio: bool = False
    pass_number: float = 1  # 1, 1.5, or 2

    # Cluster-level (for 2D fits)
    cluster_mean_lg_f1: float = 0.0
    cluster_mean_lg_f2: float = 0.0
    cluster_r_squared: float = 0.0
    cluster_chi2: float = 0.0
    cluster_iterations: int = 0

    def to_dict(self) -> Dict:
        """Convert to dictionary for JSON serialization."""
        return asdict(self)


class ComprehensiveTrainingCollector:
    """
    Comprehensive training data collector for ML model training.

    Collects ~143 parameters from spectrum processing including:
    - Spectrum features and quality metrics
    - Peak distribution and clustering
    - PASS1 learned statistics
    - Adaptive optimization results
    - PASS1-bis improvements
    - Per-peak fitting results
    """

    VERSION = "2.0.0"
    DEFAULT_MIN_R2 = 0.70  # Lower threshold for diversity

    def __init__(
        self,
        storage_dir: Optional[Path] = None,
        min_r2: float = DEFAULT_MIN_R2,
        auto_save: bool = True,
    ):
        """
        Initialize the comprehensive collector.

        Parameters
        ----------
        storage_dir : Path, optional
            Custom storage directory. If None, uses platform default.
        min_r2 : float
            Minimum median R² to accept a sample (default 0.70)
        auto_save : bool
            Whether to auto-save after each collection
        """
        self.storage_dir = storage_dir
        self.min_r2 = min_r2
        self.auto_save = auto_save

        # Current session samples
        self.session_samples: List[ComprehensiveTrainingSample] = []
        self.session_start = datetime.now()

        # Statistics
        self.n_collected = 0
        self.n_rejected = 0

    def collect_spectrum_processing(
        self,
        spectrum_data: np.ndarray,
        peak_list: List[Dict],
        fit_results: List[Dict],
        pass1_statistics: Optional[Dict] = None,
        adaptive_results: Optional[Dict] = None,
        pass1_results: Optional[List[Dict]] = None,
        pass1bis_results: Optional[List[Dict]] = None,
        pass2_results: Optional[List[Dict]] = None,
        cluster_info: Optional[Dict] = None,
        timing_info: Optional[Dict] = None,
        spectrum_name: str = "",
        nucleus_type: str = "15N",
        field_strength_mhz: float = 600.0,
        ppm_ranges: Optional[Dict] = None,
        processing_mode: str = "sequential",
        is_reference: bool = False,
        series_index: int = 0,
        detection_params: Optional[Dict] = None,
        detection_statistics: Optional[Dict] = None,
    ) -> bool:
        """
        Collect comprehensive training data from spectrum processing.

        Parameters
        ----------
        spectrum_data : np.ndarray
            2D spectrum data
        peak_list : list
            List of peak dictionaries
        fit_results : list
            Final fit results for all peaks
        pass1_statistics : dict, optional
            Statistics from PASS1 (lw_f1_median, lw_f2_median, etc.)
        adaptive_results : dict, optional
            Results from adaptive optimization
        pass1_results : list, optional
            Fit results from PASS1 (before adaptive)
        pass1bis_results : list, optional
            Fit results from PASS1-bis (after adaptive)
        pass2_results : list, optional
            Fit results from PASS2
        cluster_info : dict, optional
            Clustering information
        timing_info : dict, optional
            Processing timing information
        spectrum_name : str
            Spectrum file name
        nucleus_type : str
            Nucleus type ("15N" or "13C")
        field_strength_mhz : float
            Spectrometer frequency
        ppm_ranges : dict, optional
            PPM ranges {'f1': (min, max), 'f2': (min, max)}
        processing_mode : str
            "sequential" or "parallel"
        is_reference : bool
            Whether this is a series reference spectrum
        series_index : int
            Position in series
        detection_params : dict, optional
            Peak detection parameters: search_window_x, search_window_y,
            noise_threshold, noise_level
        detection_statistics : dict, optional
            Detection results: total_peaks, detected_peaks, reference_retained,
            dummy_peaks_added, detection_rate

        Returns
        -------
        bool
            True if sample was accepted
        """
        # Calculate median R² from fit results (handle both r_squared and avg_r_squared naming)
        r_squared_values = [r.get('r_squared', r.get('avg_r_squared', 0)) for r in fit_results if r.get('success', False)]
        if not r_squared_values:
            logger.debug("Rejected: no successful fits")
            self.n_rejected += 1
            return False

        median_r2 = np.median(r_squared_values)
        if median_r2 < self.min_r2:
            logger.debug(f"Rejected: median R² {median_r2:.3f} < {self.min_r2}")
            self.n_rejected += 1
            return False

        # Create sample
        sample = ComprehensiveTrainingSample(
            timestamp=datetime.now().isoformat(),
            spectrum_name=spectrum_name,
            processing_mode=processing_mode,
            is_reference_spectrum=is_reference,
            series_index=series_index,
            nucleus_type=nucleus_type,
            field_strength_mhz=field_strength_mhz,
        )

        # Extract spectrum features
        self._extract_spectrum_features(sample, spectrum_data, peak_list, ppm_ranges)

        # Extract spectrum quality (use fit_results for intensity values)
        self._extract_spectrum_quality(sample, spectrum_data, peak_list, fit_results)

        # Extract clustering info (pass fit_results for overlap_group_size extraction)
        self._extract_clustering_info(sample, cluster_info, peak_list, fit_results)

        # Extract protein characteristics
        self._extract_protein_characteristics(sample)

        # Extract pre-fitting clustering estimates (from peak positions only)
        self._extract_prefitting_clustering(sample, peak_list)

        # Extract peak list intensity stats (if available)
        self._extract_peaklist_intensity_stats(sample, peak_list)

        # Extract crowding metrics
        self._extract_crowding_metrics(sample, peak_list)

        # Extract detection info (reference list → detected peaks transformation)
        self._extract_detection_info(sample, fit_results, detection_params, detection_statistics)

        # Extract PASS1 statistics
        if pass1_statistics:
            self._extract_pass1_statistics(sample, pass1_statistics, pass1_results)

        # Extract adaptive optimization results
        if adaptive_results:
            self._extract_adaptive_results(sample, adaptive_results)

        # Extract PASS1-bis improvement
        if pass1_results and pass1bis_results:
            self._extract_pass1bis_improvement(sample, pass1_results, pass1bis_results, pass1_statistics)

        # Extract re-clustering info
        if cluster_info and adaptive_results:
            self._extract_recluster_info(sample, cluster_info)

        # Extract PASS2 results (use fit_results as pass2_results if not provided)
        effective_pass2 = pass2_results if pass2_results else fit_results
        if effective_pass2:
            self._extract_pass2_results(sample, effective_pass2)

        # Extract overall improvement
        self._extract_overall_improvement(sample, pass1_results, pass1bis_results, effective_pass2)

        # Extract timing
        if timing_info:
            self._extract_timing(sample, timing_info)

        # Extract per-peak data (pass PASS1 stats for initial linewidths)
        self._extract_peak_data(sample, fit_results, peak_list, cluster_info, pass1_statistics)

        # Add to session
        self.session_samples.append(sample)
        self.n_collected += 1

        logger.info(f"Collected sample: {spectrum_name}, {len(fit_results)} peaks, R²={median_r2:.3f}")

        # Auto-save
        if self.auto_save:
            self.save()

        return True

    def _extract_spectrum_features(
        self,
        sample: ComprehensiveTrainingSample,
        spectrum_data: np.ndarray,
        peak_list: List[Dict],
        ppm_ranges: Optional[Dict],
    ):
        """Extract spectrum-level features."""
        # Noise estimation (corner method)
        ny, nx = spectrum_data.shape
        size_y = max(1, ny // 10)
        size_x = max(1, nx // 10)
        corners = [
            spectrum_data[:size_y, :size_x],
            spectrum_data[:size_y, -size_x:],
            spectrum_data[-size_y:, :size_x],
            spectrum_data[-size_y:, -size_x:],
        ]
        noise_level = np.median([np.std(c) for c in corners])
        sample.noise_level = float(noise_level)

        # Dynamic range
        sample.dynamic_range = float(np.max(spectrum_data) - np.min(spectrum_data))

        # SNR estimate
        if noise_level > 0:
            sample.snr_estimate = float(np.max(spectrum_data) / noise_level)

        # Peak count and density
        sample.peak_count = len(peak_list)

        # Calculate shift ranges from actual peak positions (for IDP detection)
        # Position_Y = F1 (15N/13C), Position_X = F2 (1H)
        f1_positions = []
        f2_positions = []
        for peak in peak_list:
            f1_pos = peak.get('Position_Y', peak.get('y_ppm', peak.get('pos_y')))
            f2_pos = peak.get('Position_X', peak.get('x_ppm', peak.get('pos_x')))
            if f1_pos is not None:
                f1_positions.append(float(f1_pos))
            if f2_pos is not None:
                f2_positions.append(float(f2_pos))

        if f1_positions:
            sample.shift_range_f1_min = float(np.min(f1_positions))
            sample.shift_range_f1_max = float(np.max(f1_positions))
        if f2_positions:
            sample.shift_range_f2_min = float(np.min(f2_positions))
            sample.shift_range_f2_max = float(np.max(f2_positions))

        # Peak density uses spectrum ppm_ranges (not peak positions)
        if ppm_ranges:
            f1_range = ppm_ranges.get('f1', (0, 1))
            f2_range = ppm_ranges.get('f2', (0, 1))
            area = abs(f1_range[1] - f1_range[0]) * abs(f2_range[1] - f2_range[0])
            if area > 0:
                sample.peak_density = float(len(peak_list) / area)

    def _extract_spectrum_quality(
        self,
        sample: ComprehensiveTrainingSample,
        spectrum_data: np.ndarray,
        peak_list: List[Dict],
        fit_results: Optional[List[Dict]] = None,
    ):
        """Extract spectrum quality metrics.

        Intensities are extracted from fit_results (preferred) or peak_list.
        fit_results contains: intensity (fitted), volume, height, detected_intensity
        """
        # Get intensities from fit_results (preferred) or peak_list
        volumes = []
        heights = []
        detected_intensities = []

        if fit_results:
            for result in fit_results:
                if not result.get('success', False):
                    continue
                # Volume (fitted integral for normalized Voigt)
                volume = result.get('volume', result.get('intensity', result.get('Volume', 0)))
                if volume > 0:
                    volumes.append(volume)
                # Height (peak maximum at center, different from volume!)
                height = result.get('height', result.get('amplitude', result.get('Height', 0)))
                if height > 0:
                    heights.append(height)
                # Detected intensity (from peak detection before fitting)
                det_int = result.get('detected_intensity', 0)
                if det_int and det_int > 0:
                    detected_intensities.append(det_int)

        # Fallback to peak_list if no fit_results
        if not volumes:
            for peak in peak_list:
                vol = peak.get('Volume', peak.get('volume', peak.get('Intensity', peak.get('intensity', 0))))
                if vol > 0:
                    volumes.append(vol)

        # Helper to calculate statistics for a value array
        def calc_stats(values):
            if not values:
                return 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
            arr = np.array(values)
            mean_val = float(np.mean(arr))
            std_val = float(np.std(arr))
            cv = float(std_val / mean_val) if mean_val > 0 else 0.0
            skew = float(scipy_stats.skew(arr)) if len(arr) > 2 else 0.0
            kurt = float(scipy_stats.kurtosis(arr)) if len(arr) > 2 else 0.0
            min_val = np.min(arr)
            dyn_range = float(np.max(arr) / min_val) if min_val > 0 else 0.0
            return mean_val, std_val, skew, kurt, cv, dyn_range

        # Calculate volume statistics
        if volumes:
            stats = calc_stats(volumes)
            sample.volume_mean = stats[0]
            sample.volume_std = stats[1]
            sample.volume_skewness = stats[2]
            sample.volume_kurtosis = stats[3]
            sample.volume_cv = stats[4]
            sample.volume_dynamic_range = stats[5]

        # Calculate height statistics (peak maximum, different from volume)
        if heights:
            stats = calc_stats(heights)
            sample.height_mean = stats[0]
            sample.height_std = stats[1]
            sample.height_skewness = stats[2]
            sample.height_kurtosis = stats[3]
            sample.height_cv = stats[4]
            sample.height_dynamic_range = stats[5]

        # Calculate detected intensity statistics (pre-fitting peak detection)
        if detected_intensities:
            stats = calc_stats(detected_intensities)
            sample.detected_intensity_mean = stats[0]
            sample.detected_intensity_std = stats[1]
            sample.detected_intensity_skewness = stats[2]
            sample.detected_intensity_kurtosis = stats[3]
            sample.detected_intensity_cv = stats[4]
            sample.detected_intensity_dynamic_range = stats[5]

        # Baseline estimation
        ny, nx = spectrum_data.shape
        size_y = max(1, ny // 10)
        size_x = max(1, nx // 10)
        corners = np.concatenate([
            spectrum_data[:size_y, :size_x].ravel(),
            spectrum_data[:size_y, -size_x:].ravel(),
            spectrum_data[-size_y:, :size_x].ravel(),
            spectrum_data[-size_y:, -size_x:].ravel(),
        ])
        sample.baseline_level = float(np.median(corners))
        sample.baseline_std = float(np.std(corners))

    def _extract_clustering_info(
        self,
        sample: ComprehensiveTrainingSample,
        cluster_info: Optional[Dict],
        peak_list: List[Dict],
        fit_results: Optional[List[Dict]] = None,
    ):
        """Extract clustering information.

        Computes cluster stats from fit results using overlap_group_size,
        since cluster_info may not be populated.
        """
        # Try to extract from fit_results using overlap_group_size
        if fit_results:
            isolated_count = 0
            clustered_count = 0
            cluster_sizes_seen = []

            for result in fit_results:
                overlap_size = int(result.get('overlap_group_size', 1))
                if overlap_size <= 1:
                    isolated_count += 1
                else:
                    clustered_count += 1
                    cluster_sizes_seen.append(overlap_size)

            sample.n_isolated_peaks = isolated_count
            sample.n_clustered_peaks = clustered_count

            # Estimate number of unique clusters from sizes
            # (this is approximate since we don't have true cluster IDs)
            if cluster_sizes_seen:
                # Get unique cluster sizes and estimate cluster count
                unique_sizes = set(cluster_sizes_seen)
                sample.cluster_sizes = sorted(unique_sizes, reverse=True)
                sample.max_cluster_size = max(cluster_sizes_seen)
                sample.mean_cluster_size = float(np.mean(cluster_sizes_seen))
                # Rough estimate of cluster count
                sample.n_clusters = len([s for s in cluster_sizes_seen if cluster_sizes_seen.count(s) == s])
                if sample.n_clusters == 0:
                    sample.n_clusters = len(unique_sizes)

        # Fall back to cluster_info if available
        elif cluster_info:
            isolated = cluster_info.get('isolated_clusters', [])
            multi = cluster_info.get('multi_peak_clusters', [])

            sample.n_isolated_peaks = len(isolated)
            sample.n_clustered_peaks = sum(len(c) for c in multi)
            sample.n_clusters = len(multi)

            if multi:
                cluster_sizes = [len(c) for c in multi]
                sample.cluster_sizes = cluster_sizes
                sample.max_cluster_size = max(cluster_sizes)
                sample.mean_cluster_size = float(np.mean(cluster_sizes))
        else:
            sample.n_isolated_peaks = len(peak_list)

        total_peaks = sample.n_isolated_peaks + sample.n_clustered_peaks
        if total_peaks > 0:
            sample.fraction_clustered = float(sample.n_clustered_peaks / total_peaks)

        # Peak separations (lunaNMR uses Position_Y for F1, Position_X for F2)
        if len(peak_list) > 1:
            positions_f1 = [p.get('Position_Y', p.get('y_ppm', p.get('pos_y', 0))) for p in peak_list]
            positions_f2 = [p.get('Position_X', p.get('x_ppm', p.get('pos_x', 0))) for p in peak_list]

            # Filter out zeros (missing data)
            valid_f1 = [pos for pos in positions_f1 if pos != 0]
            valid_f2 = [pos for pos in positions_f2 if pos != 0]

            if len(valid_f1) > 1:
                separations_f1 = []
                for i in range(len(valid_f1)):
                    for j in range(i + 1, len(valid_f1)):
                        separations_f1.append(abs(valid_f1[i] - valid_f1[j]))
                if separations_f1:
                    sample.mean_peak_separation_f1 = float(np.mean(separations_f1))
                    sample.min_peak_separation_f1 = float(np.min(separations_f1))

            if len(valid_f2) > 1:
                separations_f2 = []
                for i in range(len(valid_f2)):
                    for j in range(i + 1, len(valid_f2)):
                        separations_f2.append(abs(valid_f2[i] - valid_f2[j]))
                if separations_f2:
                    sample.mean_peak_separation_f2 = float(np.mean(separations_f2))
                    sample.min_peak_separation_f2 = float(np.min(separations_f2))

    def _extract_protein_characteristics(self, sample: ComprehensiveTrainingSample):
        """Extract derived protein characteristics."""
        # Dispersions
        sample.f1_dispersion = sample.shift_range_f1_max - sample.shift_range_f1_min
        sample.f2_dispersion = sample.shift_range_f2_max - sample.shift_range_f2_min

        # Dispersion ratio (spectral shape indicator)
        if sample.f2_dispersion > 0:
            sample.dispersion_ratio = sample.f1_dispersion / sample.f2_dispersion
        else:
            sample.dispersion_ratio = 0.0

        # IDP detection (narrow 1H dispersion)
        sample.is_idp_like = sample.f2_dispersion < 1.5
        sample.folding_indicator = sample.f2_dispersion

        # Protein size estimate from peak count
        if sample.peak_count < 80:
            sample.estimated_protein_size = "small"
        elif sample.peak_count < 200:
            sample.estimated_protein_size = "medium"
        else:
            sample.estimated_protein_size = "large"

    def _extract_prefitting_clustering(
        self,
        sample: ComprehensiveTrainingSample,
        peak_list: List[Dict],
    ):
        """Extract pre-fitting clustering estimates from peak positions.

        Uses typical linewidth thresholds to estimate overlap before fitting.
        For 15N-HSQC: radF1=0.4 ppm, radF2=0.04 ppm (standard values).
        """
        if len(peak_list) < 2:
            return

        # Get peak positions
        positions = []
        for peak in peak_list:
            f1 = peak.get('Position_Y', peak.get('y_ppm', peak.get('pos_y')))
            f2 = peak.get('Position_X', peak.get('x_ppm', peak.get('pos_x')))
            if f1 is not None and f2 is not None:
                positions.append((float(f1), float(f2)))

        if len(positions) < 2:
            return

        # Typical radii for overlap detection (before we know actual linewidths)
        # These are conservative estimates based on typical 15N-HSQC values
        if sample.nucleus_type == '15N':
            threshold_f1 = 0.4  # ppm, typical 15N linewidth radius
            threshold_f2 = 0.04  # ppm, typical 1H linewidth radius
        else:  # 13C
            threshold_f1 = 0.8  # ppm, typical 13C linewidth radius
            threshold_f2 = 0.04  # ppm, typical 1H linewidth radius

        n_close_f1 = 0
        n_close_f2 = 0
        n_close_2d = 0
        peaks_with_neighbors = set()

        n_peaks = len(positions)
        for i in range(n_peaks):
            for j in range(i + 1, n_peaks):
                f1_i, f2_i = positions[i]
                f1_j, f2_j = positions[j]

                dist_f1 = abs(f1_i - f1_j)
                dist_f2 = abs(f2_i - f2_j)

                # Check F1 proximity
                if dist_f1 < 2 * threshold_f1:
                    n_close_f1 += 1

                # Check F2 proximity
                if dist_f2 < 2 * threshold_f2:
                    n_close_f2 += 1

                # Check 2D elliptical distance
                elliptical_dist = np.sqrt(
                    (dist_f1 / threshold_f1) ** 2 + (dist_f2 / threshold_f2) ** 2
                )
                if elliptical_dist < 2.0:  # Within 2 radii
                    n_close_2d += 1
                    peaks_with_neighbors.add(i)
                    peaks_with_neighbors.add(j)

        sample.n_close_pairs_f1 = n_close_f1
        sample.n_close_pairs_f2 = n_close_f2
        sample.n_close_pairs_2d = n_close_2d
        sample.fraction_potentially_overlapping = len(peaks_with_neighbors) / n_peaks

    def _extract_peaklist_intensity_stats(
        self,
        sample: ComprehensiveTrainingSample,
        peak_list: List[Dict],
    ):
        """Extract intensity statistics from peak list (before fitting).

        Peak lists may contain reference intensities/heights from detection
        or from a reference spectrum.
        """
        intensities = []
        for peak in peak_list:
            # Try various intensity field names
            intensity = peak.get('Intensity', peak.get('intensity',
                         peak.get('Height', peak.get('height',
                         peak.get('Volume', peak.get('volume'))))))
            if intensity is not None and intensity > 0:
                intensities.append(float(intensity))

        if len(intensities) < 2:
            return

        arr = np.array(intensities)
        sample.peaklist_intensity_mean = float(np.mean(arr))
        sample.peaklist_intensity_std = float(np.std(arr))

        if sample.peaklist_intensity_mean > 0:
            sample.peaklist_intensity_cv = float(
                sample.peaklist_intensity_std / sample.peaklist_intensity_mean
            )

        min_val = np.min(arr)
        if min_val > 0:
            sample.peaklist_intensity_dynamic_range = float(np.max(arr) / min_val)

    def _extract_crowding_metrics(
        self,
        sample: ComprehensiveTrainingSample,
        peak_list: List[Dict],
    ):
        """Extract local crowding metrics from peak positions.

        Identifies regions with high peak density.
        """
        if len(peak_list) < 3:
            return

        # Get peak positions
        positions = []
        for peak in peak_list:
            f1 = peak.get('Position_Y', peak.get('y_ppm', peak.get('pos_y')))
            f2 = peak.get('Position_X', peak.get('x_ppm', peak.get('pos_x')))
            if f1 is not None and f2 is not None:
                positions.append((float(f1), float(f2)))

        if len(positions) < 3:
            return

        f1_coords = np.array([p[0] for p in positions])
        f2_coords = np.array([p[1] for p in positions])

        # Define window size based on nucleus type
        if sample.nucleus_type == '15N':
            window_f1 = 2.0  # ppm
            window_f2 = 0.2  # ppm
        else:  # 13C
            window_f1 = 4.0  # ppm
            window_f2 = 0.2  # ppm

        window_area = window_f1 * window_f2

        # Sliding window to find max density
        # Use peak positions as window centers
        max_density = 0.0
        densities = []

        for f1_center, f2_center in positions:
            # Count peaks in window centered at this peak
            in_window = np.sum(
                (np.abs(f1_coords - f1_center) <= window_f1 / 2) &
                (np.abs(f2_coords - f2_center) <= window_f2 / 2)
            )
            density = in_window / window_area
            densities.append(density)
            if density > max_density:
                max_density = density

        sample.max_local_density = float(max_density)

        # Fraction in densest 10% of density values
        if densities:
            threshold = np.percentile(densities, 90)
            n_in_hotspot = np.sum(np.array(densities) >= threshold)
            sample.crowding_hotspot_fraction = float(n_in_hotspot / len(positions))

    def _extract_detection_info(
        self,
        sample: ComprehensiveTrainingSample,
        fit_results: List[Dict],
        detection_params: Optional[Dict],
        detection_statistics: Optional[Dict],
    ):
        """Extract peak detection parameters and reference vs detected statistics.

        Parameters
        ----------
        sample : ComprehensiveTrainingSample
            Sample to populate
        fit_results : list
            Fit results containing per-peak detection flags
        detection_params : dict, optional
            Detection parameters: search_window_x, search_window_y,
            noise_threshold, noise_level
        detection_statistics : dict, optional
            Detection results: total_peaks, detected_peaks, reference_retained,
            dummy_peaks_added, detection_rate
        """
        # Extract detection parameters
        if detection_params:
            sample.detection_search_window_x = float(detection_params.get('search_window_x', 0.01))
            sample.detection_search_window_y = float(detection_params.get('search_window_y', 0.04))
            sample.detection_noise_threshold = float(detection_params.get('noise_threshold', 3.0))
            sample.detection_noise_level = float(detection_params.get('noise_level', 0))

        # Extract detection statistics
        if detection_statistics:
            sample.reference_peak_count = int(detection_statistics.get('total_peaks', 0))
            sample.detected_peak_count = int(detection_statistics.get('detected_peaks', 0))
            sample.reference_retained_count = int(detection_statistics.get('reference_retained', 0))
            sample.dummy_peaks_added = int(detection_statistics.get('dummy_peaks_added', 0))
            sample.detection_rate = float(detection_statistics.get('detection_rate', 0))

        # Calculate distance from reference statistics from fit_results
        distances = []
        for result in fit_results:
            if result.get('success', False):
                dist = result.get('distance_from_reference', 0)
                if dist is not None and dist > 0:
                    distances.append(float(dist))

        if distances:
            sample.mean_distance_from_reference = float(np.mean(distances))
            sample.max_distance_from_reference = float(np.max(distances))
            sample.std_distance_from_reference = float(np.std(distances))

    def _extract_pass1_statistics(
        self,
        sample: ComprehensiveTrainingSample,
        stats: Dict,
        pass1_results: Optional[List[Dict]],
    ):
        """Extract PASS1 learned statistics."""
        sample.pass1_lw_f1_median = float(stats.get('lw_f1_median', 0))
        sample.pass1_lw_f2_median = float(stats.get('lw_f2_median', 0))
        sample.pass1_lw_f1_std = float(stats.get('lw_f1_std', 0))
        sample.pass1_lw_f2_std = float(stats.get('lw_f2_std', 0))
        sample.pass1_lw_f1_min = float(stats.get('lw_f1_min', 0))
        sample.pass1_lw_f1_max = float(stats.get('lw_f1_max', 0))
        sample.pass1_lw_f2_min = float(stats.get('lw_f2_min', 0))
        sample.pass1_lw_f2_max = float(stats.get('lw_f2_max', 0))
        sample.pass1_n_good_peaks = int(stats.get('n_good_peaks', 0))

        # Get mean R² from pass1_results if available, otherwise from stats
        if pass1_results:
            r2_values = [r.get('r_squared', r.get('avg_r_squared', 0)) for r in pass1_results if r.get('success', False)]
            if r2_values:
                sample.pass1_mean_r_squared = float(np.mean(r2_values))
        elif 'mean_r_squared' in stats:
            sample.pass1_mean_r_squared = float(stats.get('mean_r_squared', 0))

    def _extract_adaptive_results(
        self,
        sample: ComprehensiveTrainingSample,
        results: Dict,
    ):
        """Extract adaptive optimization results."""
        sample.adaptive_success = bool(results.get('success', False))
        sample.adaptive_fallback_reason = str(results.get('fallback_reason', ''))
        sample.adaptive_radF1 = float(results.get('radF1', 0))
        sample.adaptive_radF2 = float(results.get('radF2', 0))
        sample.adaptive_overlap_threshold_x = float(results.get('overlap_threshold_x', 0))
        sample.adaptive_overlap_threshold_y = float(results.get('overlap_threshold_y', 0))
        sample.adaptive_multiplier_f1 = float(results.get('multiplier_f1', 1.5))
        sample.adaptive_multiplier_f2 = float(results.get('multiplier_f2', 1.5))
        sample.adaptive_train_score = float(results.get('train_score', 0) or 0)
        sample.adaptive_validation_score = float(results.get('validation_score', 0) or 0)
        sample.adaptive_generalization_gap = float(results.get('generalization_gap', 0) or 0)
        sample.adaptive_n_peaks_used = int(results.get('n_peaks_used', 0))
        # Note: n_train and n_val are not stored separately in adaptive optimizer
        # The optimizer uses all n_peaks_used peaks with internal train/val split
        sample.adaptive_n_train = 0  # Not available from optimizer
        sample.adaptive_n_val = 0    # Not available from optimizer
        sample.adaptive_search_history = results.get('search_history', [])

    def _extract_pass1bis_improvement(
        self,
        sample: ComprehensiveTrainingSample,
        pass1_results: List[Dict],
        pass1bis_results: List[Dict],
        pass1_statistics: Optional[Dict],
    ):
        """Extract PASS1 to PASS1-bis improvement metrics."""
        sample.pass1bis_enabled = True
        sample.pass1bis_n_peaks_refit = len(pass1bis_results)

        # R² comparison
        pass1_r2 = [r.get('r_squared', r.get('avg_r_squared', 0)) for r in pass1_results if r.get('success', False)]
        pass1bis_r2 = [r.get('r_squared', r.get('avg_r_squared', 0)) for r in pass1bis_results if r.get('success', False)]

        if pass1_r2:
            sample.pass1bis_mean_r2_before = float(np.mean(pass1_r2))
        if pass1bis_r2:
            sample.pass1bis_mean_r2_after = float(np.mean(pass1bis_r2))

        sample.pass1bis_r2_improvement = sample.pass1bis_mean_r2_after - sample.pass1bis_mean_r2_before

        # Count improved/degraded/unchanged
        n_improved = 0
        n_degraded = 0
        n_unchanged = 0

        # Match by peak position if possible
        for i, r1 in enumerate(pass1_results):
            if i < len(pass1bis_results):
                r2 = pass1bis_results[i]
                r1_val = r1.get('r_squared', r1.get('avg_r_squared', 0))
                r2_val = r2.get('r_squared', r2.get('avg_r_squared', 0))
                diff = r2_val - r1_val
                if diff > 0.01:
                    n_improved += 1
                elif diff < -0.01:
                    n_degraded += 1
                else:
                    n_unchanged += 1

        sample.pass1bis_n_improved = n_improved
        sample.pass1bis_n_degraded = n_degraded
        sample.pass1bis_n_unchanged = n_unchanged

    def _extract_recluster_info(
        self,
        sample: ComprehensiveTrainingSample,
        cluster_info: Dict,
    ):
        """Extract re-clustering results."""
        if 'recluster' in cluster_info:
            recluster = cluster_info['recluster']
            sample.recluster_performed = True
            sample.recluster_n_originally_clustered = int(recluster.get('n_originally_clustered', 0))
            sample.recluster_n_now_isolated = int(recluster.get('n_now_isolated', 0))
            sample.recluster_n_still_clustered = int(recluster.get('n_still_clustered', 0))
            sample.recluster_new_cluster_count = int(recluster.get('new_cluster_count', 0))
            sample.recluster_cluster_size_change = float(recluster.get('cluster_size_change', 0))

    def _extract_pass2_results(
        self,
        sample: ComprehensiveTrainingSample,
        pass2_results: List[Dict],
    ):
        """Extract PASS2 results."""
        r2_values = [r.get('r_squared', r.get('avg_r_squared', 0)) for r in pass2_results if r.get('success', False)]
        if r2_values:
            sample.pass2_mean_r_squared = float(np.mean(r2_values))
        sample.pass2_n_peaks_fitted = len([r for r in pass2_results if r.get('success', False)])

        # Count clusters - approximate by grouping peaks with same overlap_group_size
        # Each cluster has N peaks with overlap_group_size=N, so n_clusters ≈ n_clustered_peaks / avg_size
        clustered_peaks = [r for r in pass2_results if r.get('overlap_group_size', 1) > 1 and r.get('success', False)]
        if clustered_peaks:
            total_size = sum(r.get('overlap_group_size', 1) for r in clustered_peaks)
            avg_size = total_size / len(clustered_peaks)
            sample.pass2_n_clusters_fitted = int(len(clustered_peaks) / avg_size) if avg_size > 0 else 0
        else:
            sample.pass2_n_clusters_fitted = 0

    def _extract_overall_improvement(
        self,
        sample: ComprehensiveTrainingSample,
        pass1_results: Optional[List[Dict]],
        pass1bis_results: Optional[List[Dict]],
        pass2_results: Optional[List[Dict]],
    ):
        """Extract overall improvement metrics."""
        if pass1_results:
            r2 = [r.get('r_squared', r.get('avg_r_squared', 0)) for r in pass1_results if r.get('success', False)]
            if r2:
                sample.overall_r2_pass1 = float(np.mean(r2))
        else:
            # Fallback to pass1_mean_r_squared if pass1_results not available
            if sample.pass1_mean_r_squared > 0:
                sample.overall_r2_pass1 = sample.pass1_mean_r_squared

        if pass1bis_results:
            r2 = [r.get('r_squared', r.get('avg_r_squared', 0)) for r in pass1bis_results if r.get('success', False)]
            if r2:
                sample.overall_r2_pass1bis = float(np.mean(r2))

        if pass2_results:
            r2 = [r.get('r_squared', r.get('avg_r_squared', 0)) for r in pass2_results if r.get('success', False)]
            if r2:
                sample.overall_r2_pass2 = float(np.mean(r2))

        sample.improvement_pass1_to_pass1bis = sample.overall_r2_pass1bis - sample.overall_r2_pass1
        sample.improvement_pass1bis_to_pass2 = sample.overall_r2_pass2 - sample.overall_r2_pass1bis
        sample.total_improvement = sample.overall_r2_pass2 - sample.overall_r2_pass1

    def _extract_timing(
        self,
        sample: ComprehensiveTrainingSample,
        timing_info: Dict,
    ):
        """Extract timing information."""
        sample.total_processing_time = float(timing_info.get('total', 0))
        sample.pass1_time = float(timing_info.get('pass1', 0))
        sample.adaptive_time = float(timing_info.get('adaptive', 0))
        sample.pass1bis_time = float(timing_info.get('pass1bis', 0))
        sample.pass2_time = float(timing_info.get('pass2', 0))

    def _extract_peak_data(
        self,
        sample: ComprehensiveTrainingSample,
        fit_results: List[Dict],
        peak_list: List[Dict],
        cluster_info: Optional[Dict],
        pass1_statistics: Optional[Dict] = None,
    ):
        """Extract per-peak training data."""
        peaks_data = []

        for i, result in enumerate(fit_results):
            peak_data = PeakTrainingData(
                peak_index=i,
                peak_id=str(result.get('peak_id', i)),
                assignment=str(result.get('assignment', '')),
            )

            # Fitted parameters - handle multiple naming conventions
            # pos_f1/pos_f2 (PS2D), center_y/center_x (converted), peak_y/peak_x (legacy)
            peak_data.pos_f1 = float(result.get('pos_f1', result.get('center_y', result.get('peak_y', 0))))
            peak_data.pos_f2 = float(result.get('pos_f2', result.get('center_x', result.get('peak_x', 0))))

            # Linewidth parameters - handle naming conventions
            # PS2D 2D: lw_lor_f1/lw_gau_f1 or gamma_y/sigma_y at top level
            # 1D fits: F1 params in y_fit nested dict, F2 params in x_fit or top-level gamma/sigma
            y_fit = result.get('y_fit', {})
            x_fit = result.get('x_fit', {})

            # F1 (Y-dimension) linewidths - check y_fit for 1D fits
            peak_data.lw_lor_f1 = float(result.get('lw_lor_f1', result.get('gamma_y', y_fit.get('gamma', 0))))
            peak_data.lw_gau_f1 = float(result.get('lw_gau_f1', result.get('sigma_y', y_fit.get('sigma', 0))))

            # F2 (X-dimension) linewidths - top-level gamma/sigma are F2 values for 1D fits
            peak_data.lw_lor_f2 = float(result.get('lw_lor_f2', result.get('gamma_x', result.get('gamma', x_fit.get('gamma', 0)))))
            peak_data.lw_gau_f2 = float(result.get('lw_gau_f2', result.get('sigma_x', result.get('sigma', x_fit.get('sigma', 0)))))

            # Total linewidths
            peak_data.lw_total_f1 = peak_data.lw_lor_f1 + peak_data.lw_gau_f1
            peak_data.lw_total_f2 = peak_data.lw_lor_f2 + peak_data.lw_gau_f2

            # L/G ratios
            if peak_data.lw_total_f1 > 0:
                peak_data.lg_ratio_f1 = peak_data.lw_lor_f1 / peak_data.lw_total_f1
            if peak_data.lw_total_f2 > 0:
                peak_data.lg_ratio_f2 = peak_data.lw_lor_f2 / peak_data.lw_total_f2

            # Intensity metrics
            # volume = fitted integral (normalized Voigt parameter)
            # height = peak maximum at center (different from volume!)
            peak_data.volume = float(result.get('volume', result.get('intensity', result.get('Volume', 0))))
            peak_data.height = float(result.get('height', result.get('amplitude', result.get('Height', 0))))
            peak_data.detected_intensity = float(result.get('detected_intensity', 0))

            # Initial positions from peak list (lunaNMR: Position_Y=F1, Position_X=F2)
            if i < len(peak_list):
                initial = peak_list[i]
                peak_data.initial_pos_f1 = float(initial.get('Position_Y', initial.get('y_ppm', initial.get('pos_y', 0))))
                peak_data.initial_pos_f2 = float(initial.get('Position_X', initial.get('x_ppm', initial.get('pos_x', 0))))
                peak_data.position_shift_f1 = peak_data.pos_f1 - peak_data.initial_pos_f1
                peak_data.position_shift_f2 = peak_data.pos_f2 - peak_data.initial_pos_f2

            # Initial linewidths from PASS1 learned statistics (used as initial guess for PASS2)
            if pass1_statistics:
                peak_data.initial_lw_f1 = float(pass1_statistics.get('lw_f1_median', 0))
                peak_data.initial_lw_f2 = float(pass1_statistics.get('lw_f2_median', 0))
                # Calculate linewidth change from initial guess
                peak_data.linewidth_change_f1 = peak_data.lw_total_f1 - peak_data.initial_lw_f1
                peak_data.linewidth_change_f2 = peak_data.lw_total_f2 - peak_data.initial_lw_f2

            # Fit quality - handle r_squared vs avg_r_squared naming
            peak_data.r_squared = float(result.get('r_squared', result.get('avg_r_squared', 0)))
            peak_data.chi2 = float(result.get('chi2', 0))
            peak_data.iterations = int(result.get('iterations', 0))
            peak_data.success = bool(result.get('success', False))

            # Determine convergence type (which criterion was met)
            if result.get('formal_convergence', False):
                peak_data.convergence_type = "formal"
            elif result.get('pragmatic_acceptance', False):
                peak_data.convergence_type = "pragmatic_r2"
            elif result.get('chi2_reduction_success', False):
                peak_data.convergence_type = "chi2_reduction"
            else:
                peak_data.convergence_type = ""

            # Fitting mode - lunaNMR uses 'method' key with values like:
            # '2d_simultaneous_multi_peak', '2d_simultaneous', 'voigt_2d_fitting', etc.
            method = result.get('method', '')
            if '2d_simultaneous' in method or '2d' in method.lower():
                peak_data.fitting_mode = '2D'
            else:
                peak_data.fitting_mode = '1D'

            # Cluster info - lunaNMR uses 'overlap_group_size' for cluster membership
            overlap_group_size = int(result.get('overlap_group_size', 1))
            peak_data.n_peaks_in_fit = overlap_group_size
            peak_data.cluster_size = overlap_group_size
            peak_data.is_clustered = overlap_group_size > 1
            peak_data.is_isolated = overlap_group_size == 1

            # Get cluster_id from cluster_idx in fit result (added in single_spectrum_processor/parallel_voigt_processor)
            # -1 means isolated peak (legacy), otherwise use actual cluster index
            cluster_idx = result.get('cluster_idx', -1 if peak_data.is_isolated else 0)
            peak_data.cluster_id = int(cluster_idx)

            # Cluster-level statistics (for 2D multi-peak fits)
            # These are calculated from all_peaks in the result
            if peak_data.is_clustered:
                # R² for the cluster fit
                peak_data.cluster_r_squared = float(result.get('r_squared', result.get('avg_r_squared', 0)))

                # Chi² and iterations from 2D fitter (if available)
                peak_data.cluster_chi2 = float(result.get('chi2', 0))
                peak_data.cluster_iterations = int(result.get('iterations', 0))

                # Calculate mean L/G ratios from all_peaks in the cluster
                all_peaks = result.get('all_peaks', [])
                if all_peaks:
                    lg_f1_values = []
                    lg_f2_values = []
                    for peak in all_peaks:
                        lw_lor_f1 = peak.get('lw_lor_f1', 0)
                        lw_gau_f1 = peak.get('lw_gau_f1', 0)
                        lw_lor_f2 = peak.get('lw_lor_f2', 0)
                        lw_gau_f2 = peak.get('lw_gau_f2', 0)
                        total_f1 = lw_lor_f1 + lw_gau_f1
                        total_f2 = lw_lor_f2 + lw_gau_f2
                        if total_f1 > 0:
                            lg_f1_values.append(lw_lor_f1 / total_f1)
                        if total_f2 > 0:
                            lg_f2_values.append(lw_lor_f2 / total_f2)
                    if lg_f1_values:
                        peak_data.cluster_mean_lg_f1 = float(np.mean(lg_f1_values))
                    if lg_f2_values:
                        peak_data.cluster_mean_lg_f2 = float(np.mean(lg_f2_values))

            # Compute tooclose_flag and heavy_overlap_flag from peak positions
            # Only for clustered peaks - isolated peaks by definition have no significant overlap
            # Uses same logic as PS2D algorithm in core_integrator.py
            # Thresholds: heavy_overlap_threshold=1.5, tooclose_threshold=0.8
            # Radii: radF1=0.4 ppm (15N), radF2=0.04 ppm (1H)
            peak_data.tooclose_flag = False
            peak_data.heavy_overlap_flag = False

            # Only compute overlap flags for clustered peaks
            if peak_data.is_clustered and peak_data.pos_f1 > 0 and peak_data.pos_f2 > 0:
                rad_f1 = 0.4  # Standard 15N radius
                rad_f2 = 0.04  # Standard 1H radius
                heavy_overlap_threshold = 1.5
                tooclose_threshold = 0.8

                # Compare to all other peaks in results
                for j, other_result in enumerate(fit_results):
                    if i == j:
                        continue
                    # Handle multiple naming conventions: pos_f1/pos_f2, center_y/center_x, peak_y/peak_x
                    # Also check nested y_fit/x_fit for 1D fits
                    other_y_fit = other_result.get('y_fit', {})
                    other_x_fit = other_result.get('x_fit', {})
                    other_f1 = float(other_result.get('pos_f1', other_result.get('center_y', other_result.get('peak_y', other_y_fit.get('center', 0)))))
                    other_f2 = float(other_result.get('pos_f2', other_result.get('center_x', other_result.get('peak_x', other_x_fit.get('center', 0)))))
                    if other_f1 <= 0 or other_f2 <= 0:
                        continue

                    # Calculate elliptical distance
                    dist_f1 = abs(peak_data.pos_f1 - other_f1) / rad_f1
                    dist_f2 = abs(peak_data.pos_f2 - other_f2) / rad_f2
                    elliptical_distance = np.sqrt(dist_f1**2 + dist_f2**2)

                    if elliptical_distance < heavy_overlap_threshold:
                        peak_data.heavy_overlap_flag = True
                    if elliptical_distance < tooclose_threshold:
                        peak_data.tooclose_flag = True

                    # Early exit if both flags set
                    if peak_data.tooclose_flag and peak_data.heavy_overlap_flag:
                        break

            # Detection info (reference list → detected peaks)
            peak_data.was_detected = result.get('detected', True)
            peak_data.distance_from_reference = float(result.get('distance_from_reference', 0))
            peak_data.distance_from_reference_x = float(result.get('distance_from_reference_x', 0))
            peak_data.distance_from_reference_y = float(result.get('distance_from_reference_y', 0))
            peak_data.distance_from_reference_elliptical = float(result.get('distance_from_reference_elliptical', 0))
            peak_data.detection_quality = str(result.get('detection_quality', 'Matched'))

            peaks_data.append(peak_data.to_dict())

        sample.peaks = peaks_data

    def save(self) -> bool:
        """Save collected samples to storage."""
        if not self.session_samples:
            return True

        try:
            data_path = get_training_data_path(self.storage_dir)

            # Load existing data
            existing_samples = []
            if data_path.exists():
                with open(data_path, 'r') as f:
                    existing_data = json.load(f)
                    existing_samples = existing_data.get('samples', [])

            # Append new samples
            for sample in self.session_samples:
                existing_samples.append(sample.to_dict())

            # Save
            data = {
                'version': self.VERSION,
                'updated': datetime.now().isoformat(),
                'n_samples': len(existing_samples),
                'samples': existing_samples,
            }

            with open(data_path, 'w') as f:
                json.dump(data, f, indent=2)

            # Save session log (before clearing session_samples)
            self._save_session_log()

            # Update metadata
            self._update_metadata(len(existing_samples))

            # Clear session samples after successful save to prevent re-adding on next auto-save
            saved_count = len(self.session_samples)
            self.session_samples = []

            logger.info(f"Saved {saved_count} samples to {data_path}")
            return True

        except Exception as e:
            logger.error(f"Failed to save training data: {e}")
            return False

    def _save_session_log(self):
        """Save session-specific log."""
        log_dir = get_session_log_dir(self.storage_dir)
        log_file = log_dir / f"{self.session_start.strftime('%Y-%m-%d_%H%M%S')}.json"

        session_data = {
            'session_start': self.session_start.isoformat(),
            'session_end': datetime.now().isoformat(),
            'n_collected': self.n_collected,
            'n_rejected': self.n_rejected,
            'samples': [s.to_dict() for s in self.session_samples],
        }

        with open(log_file, 'w') as f:
            json.dump(session_data, f, indent=2)

    def _update_metadata(self, total_samples: int):
        """Update training metadata file."""
        metadata_path = get_metadata_path(self.storage_dir)

        metadata = {
            'version': self.VERSION,
            'last_updated': datetime.now().isoformat(),
            'total_samples': total_samples,
            'min_r2_threshold': self.min_r2,
            'collection_stats': {
                'session_collected': self.n_collected,
                'session_rejected': self.n_rejected,
            },
        }

        with open(metadata_path, 'w') as f:
            json.dump(metadata, f, indent=2)

    def get_statistics(self) -> Dict:
        """Get collection statistics."""
        return {
            'session_start': self.session_start.isoformat(),
            'n_collected': self.n_collected,
            'n_rejected': self.n_rejected,
            'n_session_samples': len(self.session_samples),
        }

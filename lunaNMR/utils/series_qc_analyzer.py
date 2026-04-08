# ABOUTME: Analyzes series integration results for quality control.
# ABOUTME: Computes per-peak and per-spectrum statistics, flags anomalies.

"""
Series QC Analyzer module for LunaNMR.

Provides quality control analysis for series integration results stored in
training_data.json files. Computes per-peak statistics (CV, means), tracks
cluster stability, detects anomalies, and generates QC reports.
"""

from dataclasses import dataclass, field
from enum import Enum
from typing import List, Dict, Tuple, Optional
from pathlib import Path
import json
import numpy as np


class SeriesType(Enum):
    """Type of series affecting which QC metrics are flagged."""
    REPLICATE = "replicate"    # Volumes should be stable
    RELAXATION = "relaxation"  # T1/T2 decay expected
    TITRATION = "titration"    # Intensity changes with ligand
    KINETIC = "kinetic"        # Volume/intensity changes expected


@dataclass
class PeakQCStats:
    """QC statistics for a single peak across all spectra in series."""
    peak_id: str
    assignment: str  # SPARKY assignment (e.g., "A.2.ASP.H")

    # Time series (per-spectrum values, indexed by spectrum order)
    volumes: List[Optional[float]] = field(default_factory=list)
    heights: List[Optional[float]] = field(default_factory=list)
    r2_values: List[Optional[float]] = field(default_factory=list)
    lw_f1_values: List[Optional[float]] = field(default_factory=list)  # 15N or 13C
    lw_f2_values: List[Optional[float]] = field(default_factory=list)  # 1H
    cluster_sizes: List[int] = field(default_factory=list)
    pos_f1_values: List[Optional[float]] = field(default_factory=list)
    pos_f2_values: List[Optional[float]] = field(default_factory=list)

    # Aggregates (computed after loading)
    volume_cv: float = 0.0
    height_cv: float = 0.0
    r2_mean: float = 0.0
    r2_min: float = 0.0
    lw_f1_cv: float = 0.0
    lw_f2_cv: float = 0.0
    pos_f1_drift: float = 0.0  # max - min
    pos_f2_drift: float = 0.0

    # Cluster stability
    is_cluster_stable: bool = True
    cluster_transitions: List[Tuple[int, int, int]] = field(default_factory=list)

    # LW jump detection (sudden changes between consecutive spectra)
    lw_f1_jumps: List[Tuple[int, float]] = field(default_factory=list)  # (spectrum_idx, rel_change)
    lw_f2_jumps: List[Tuple[int, float]] = field(default_factory=list)

    # Outlier detection
    outlier_spectra: List[int] = field(default_factory=list)

    # Flags
    flags: List[str] = field(default_factory=list)


@dataclass
class SeriesQCReport:
    """Complete QC report for a series integration run."""
    source_path: str
    n_peaks: int
    n_spectra: int
    spectrum_names: List[str] = field(default_factory=list)

    # Processing context (affects expected CV values)
    processing_mode: str = "unknown"  # "parallel" or "sequential"
    nucleus_type: str = "unknown"     # "15N", "13C", etc.
    series_type: SeriesType = SeriesType.REPLICATE
    fix_positions_used: bool = False
    fix_linewidths_used: bool = False

    # Aggregates
    volume_cv_median: float = 0.0
    volume_cv_mean: float = 0.0
    r2_mean: float = 0.0
    r2_median: float = 0.0
    n_cluster_stable: int = 0
    n_lw_stable: int = 0
    fit_success_rate: float = 0.0

    # LW stability aggregates
    lw_f1_cv_median: float = 0.0
    lw_f2_cv_median: float = 0.0
    n_lw_f1_stable: int = 0  # Peaks with lw_f1_cv < threshold
    n_lw_f2_stable: int = 0
    n_lw_jumps: int = 0  # Total peaks with LW jumps

    # Position drift aggregates
    pos_f1_drift_median: float = 0.0
    pos_f2_drift_median: float = 0.0
    n_pos_f1_stable: int = 0  # Peaks with drift < threshold
    n_pos_f2_stable: int = 0

    # Per-spectrum quality (identify bad spectra)
    per_spectrum_r2_mean: List[float] = field(default_factory=list)
    per_spectrum_n_failed: List[int] = field(default_factory=list)
    outlier_spectra: List[int] = field(default_factory=list)

    # Per-peak
    peak_stats: Dict[str, PeakQCStats] = field(default_factory=dict)
    flagged_peaks: List[str] = field(default_factory=list)

    # Special categories
    dummy_peak_ids: List[str] = field(default_factory=list)
    failed_peak_ids: List[str] = field(default_factory=list)


class SeriesQCAnalyzer:
    """Analyzes training_data.json for quality control metrics."""

    # Default thresholds
    VOLUME_CV_WARN = 0.3
    HEIGHT_CV_WARN = 0.3
    R2_POOR = 0.7
    LW_CV_WARN = 0.4  # 40% CV threshold for linewidth stability
    LW_JUMP_THRESHOLD = 0.75  # 75% change between consecutive spectra
    POS_F1_DRIFT_WARN = 0.05  # ppm (for 15N/13C)
    POS_F2_DRIFT_WARN = 0.02  # ppm (for 1H)

    def __init__(self, series_type: SeriesType = SeriesType.REPLICATE):
        self.series_type = series_type
        self.data = None
        self.report = None
        self.spectrum_names: List[str] = []
        self.peak_ids: List[str] = []
        self.processing_mode: str = "unknown"
        self.nucleus_type: str = "unknown"
        self.peak_stats: Dict[str, PeakQCStats] = {}
        self.dummy_peak_ids: List[str] = []
        self.failed_peak_ids: List[str] = []
        self._source_path: str = ""

    def load_training_data(self, path: str) -> bool:
        """Load and parse training_data.json.

        Args:
            path: Path to training_data.json file

        Returns:
            True if successfully loaded, False otherwise
        """
        filepath = Path(path)
        if not filepath.exists():
            return False

        try:
            with open(filepath, 'r') as f:
                data = json.load(f)
        except json.JSONDecodeError:
            return False

        # Validate required fields
        if "samples" not in data or "n_samples" not in data:
            return False

        if not isinstance(data["samples"], list) or len(data["samples"]) == 0:
            return False

        self.data = data
        self._source_path = str(filepath)

        # Extract spectrum names
        self.spectrum_names = [
            sample.get("spectrum_name", f"spectrum_{i}")
            for i, sample in enumerate(data["samples"])
        ]

        # Extract unique peak IDs across all spectra
        peak_id_set = set()
        for sample in data["samples"]:
            for peak in sample.get("peaks", []):
                peak_id_set.add(peak.get("peak_id", ""))
        self.peak_ids = sorted(list(peak_id_set))

        # Extract processing context from first sample
        first_sample = data["samples"][0]
        self.processing_mode = first_sample.get("processing_mode", "unknown")
        self.nucleus_type = first_sample.get("nucleus_type", "unknown")

        return True

    def compute_per_peak_stats(self) -> None:
        """Compute CV, means, and flags for each peak."""
        if self.data is None:
            return

        samples = self.data["samples"]
        n_spectra = len(samples)

        # Build a lookup: peak_id -> list of peak dicts (one per spectrum)
        peak_data_by_id: Dict[str, List[Optional[dict]]] = {
            pid: [None] * n_spectra for pid in self.peak_ids
        }

        for spec_idx, sample in enumerate(samples):
            for peak in sample.get("peaks", []):
                pid = peak.get("peak_id", "")
                if pid in peak_data_by_id:
                    peak_data_by_id[pid][spec_idx] = peak

        # Track dummy peaks (peaks that were never detected) and failed peaks
        dummy_peak_ids = set()
        failed_peak_ids = set()

        # Create PeakQCStats for each peak
        for pid in self.peak_ids:
            peak_dicts = peak_data_by_id[pid]

            # Get assignment from first non-None entry
            assignment = ""
            for pd in peak_dicts:
                if pd is not None:
                    assignment = pd.get("assignment", "")
                    break

            stats = PeakQCStats(peak_id=pid, assignment=assignment)

            # Extract time series
            for pd in peak_dicts:
                if pd is None:
                    stats.volumes.append(None)
                    stats.heights.append(None)
                    stats.r2_values.append(None)
                    stats.lw_f1_values.append(None)
                    stats.lw_f2_values.append(None)
                    stats.pos_f1_values.append(None)
                    stats.pos_f2_values.append(None)
                    stats.cluster_sizes.append(0)
                else:
                    stats.volumes.append(pd.get("volume"))
                    stats.heights.append(pd.get("height"))
                    stats.r2_values.append(pd.get("r_squared"))
                    stats.lw_f1_values.append(pd.get("lw_total_f1"))
                    stats.lw_f2_values.append(pd.get("lw_total_f2"))
                    stats.pos_f1_values.append(pd.get("pos_f1"))
                    stats.pos_f2_values.append(pd.get("pos_f2"))
                    stats.cluster_sizes.append(pd.get("cluster_size", 1))

                    # Check if dummy peak
                    if not pd.get("was_detected", True):
                        dummy_peak_ids.add(pid)

                    # Check if fit failed
                    if not pd.get("success", True):
                        failed_peak_ids.add(pid)

            # Compute aggregates
            stats.volume_cv = self._compute_cv(stats.volumes)
            stats.height_cv = self._compute_cv(stats.heights)
            stats.lw_f1_cv = self._compute_cv(stats.lw_f1_values)
            stats.lw_f2_cv = self._compute_cv(stats.lw_f2_values)

            stats.r2_mean = self._compute_mean(stats.r2_values)
            stats.r2_min = self._compute_min(stats.r2_values)

            stats.pos_f1_drift = self._compute_drift(stats.pos_f1_values)
            stats.pos_f2_drift = self._compute_drift(stats.pos_f2_values)

            self.peak_stats[pid] = stats

        self.dummy_peak_ids = sorted(list(dummy_peak_ids))
        self.failed_peak_ids = sorted(list(failed_peak_ids))

    def _compute_cv(self, values: List[Optional[float]]) -> float:
        """Compute coefficient of variation (std/mean) for non-None values."""
        valid = [v for v in values if v is not None]
        if len(valid) < 2:
            return 0.0
        arr = np.array(valid)
        mean = np.mean(arr)
        if mean == 0:
            return 0.0
        return float(np.std(arr, ddof=1) / mean)

    def _compute_mean(self, values: List[Optional[float]]) -> float:
        """Compute mean for non-None values."""
        valid = [v for v in values if v is not None]
        if len(valid) == 0:
            return 0.0
        return float(np.mean(valid))

    def _compute_min(self, values: List[Optional[float]]) -> float:
        """Compute min for non-None values."""
        valid = [v for v in values if v is not None]
        if len(valid) == 0:
            return 0.0
        return float(np.min(valid))

    def _compute_drift(self, values: List[Optional[float]]) -> float:
        """Compute drift (max - min) for non-None values."""
        valid = [v for v in values if v is not None]
        if len(valid) < 2:
            return 0.0
        return float(np.max(valid) - np.min(valid))

    def analyze_cluster_stability(self) -> None:
        """Track cluster size changes across spectra."""
        for pid, stats in self.peak_stats.items():
            cluster_sizes = stats.cluster_sizes
            if len(cluster_sizes) < 2:
                stats.is_cluster_stable = True
                continue

            # Check for transitions
            transitions = []
            for i in range(1, len(cluster_sizes)):
                prev_size = cluster_sizes[i - 1]
                curr_size = cluster_sizes[i]
                if prev_size != curr_size:
                    transitions.append((i, prev_size, curr_size))

            stats.cluster_transitions = transitions
            stats.is_cluster_stable = len(transitions) == 0

    def detect_lw_jumps(self) -> None:
        """Detect sudden linewidth changes between consecutive spectra."""
        for pid, stats in self.peak_stats.items():
            # F1 jumps
            stats.lw_f1_jumps = self._find_jumps(stats.lw_f1_values)
            # F2 jumps
            stats.lw_f2_jumps = self._find_jumps(stats.lw_f2_values)

    def _find_jumps(self, values: List[Optional[float]]) -> List[Tuple[int, float]]:
        """Find indices where relative change exceeds threshold."""
        jumps = []
        prev_val = None
        prev_idx = -1

        for i, val in enumerate(values):
            if val is None or val == 0:
                continue
            if prev_val is not None and prev_val != 0:
                rel_change = abs(val - prev_val) / prev_val
                if rel_change > self.LW_JUMP_THRESHOLD:
                    jumps.append((i, rel_change))
            prev_val = val
            prev_idx = i

        return jumps

    def identify_outlier_spectra(self) -> None:
        """Find spectra with many problems (low R², many failures)."""
        # TODO: Implement if needed
        pass

    def generate_report(self) -> SeriesQCReport:
        """Generate complete QC report."""
        if self.data is None:
            return None

        n_spectra = len(self.data["samples"])
        n_peaks = len(self.peak_ids)

        report = SeriesQCReport(
            source_path=self._source_path,
            n_peaks=n_peaks,
            n_spectra=n_spectra,
            spectrum_names=self.spectrum_names.copy(),
            processing_mode=self.processing_mode,
            nucleus_type=self.nucleus_type,
            series_type=self.series_type,
        )

        # Compute aggregate statistics from peak_stats
        volume_cvs = [s.volume_cv for s in self.peak_stats.values()]
        r2_means = [s.r2_mean for s in self.peak_stats.values()]
        lw_f1_cvs = [s.lw_f1_cv for s in self.peak_stats.values() if s.lw_f1_cv > 0]
        lw_f2_cvs = [s.lw_f2_cv for s in self.peak_stats.values() if s.lw_f2_cv > 0]
        pos_f1_drifts = [s.pos_f1_drift for s in self.peak_stats.values() if s.pos_f1_drift > 0]
        pos_f2_drifts = [s.pos_f2_drift for s in self.peak_stats.values() if s.pos_f2_drift > 0]

        if volume_cvs:
            report.volume_cv_median = float(np.median(volume_cvs))
            report.volume_cv_mean = float(np.mean(volume_cvs))

        if r2_means:
            report.r2_mean = float(np.mean(r2_means))
            report.r2_median = float(np.median(r2_means))

        # LW stability aggregates
        if lw_f1_cvs:
            report.lw_f1_cv_median = float(np.median(lw_f1_cvs))
            report.n_lw_f1_stable = sum(1 for cv in lw_f1_cvs if cv <= self.LW_CV_WARN)
        if lw_f2_cvs:
            report.lw_f2_cv_median = float(np.median(lw_f2_cvs))
            report.n_lw_f2_stable = sum(1 for cv in lw_f2_cvs if cv <= self.LW_CV_WARN)

        # Count peaks with LW jumps
        report.n_lw_jumps = sum(
            1 for s in self.peak_stats.values()
            if s.lw_f1_jumps or s.lw_f2_jumps
        )

        # Position drift aggregates
        if pos_f1_drifts:
            report.pos_f1_drift_median = float(np.median(pos_f1_drifts))
            report.n_pos_f1_stable = sum(1 for d in pos_f1_drifts if d <= self.POS_F1_DRIFT_WARN)
        if pos_f2_drifts:
            report.pos_f2_drift_median = float(np.median(pos_f2_drifts))
            report.n_pos_f2_stable = sum(1 for d in pos_f2_drifts if d <= self.POS_F2_DRIFT_WARN)

        # Count cluster stable peaks
        report.n_cluster_stable = sum(
            1 for s in self.peak_stats.values() if s.is_cluster_stable
        )

        # Copy special categories
        report.dummy_peak_ids = self.dummy_peak_ids.copy()
        report.failed_peak_ids = self.failed_peak_ids.copy()

        # Copy peak stats
        report.peak_stats = self.peak_stats.copy()

        # Flag peaks based on thresholds and series type context
        # Volume/height CV only flagged for REPLICATE series
        flag_intensity_cv = (self.series_type == SeriesType.REPLICATE)

        flagged = set()
        for pid, stats in self.peak_stats.items():
            flags = []

            # Check volume CV (only for replicate series)
            if flag_intensity_cv and stats.volume_cv > self.VOLUME_CV_WARN:
                flags.append("high_vol_cv")

            # Check height CV (only for replicate series)
            if flag_intensity_cv and stats.height_cv > self.HEIGHT_CV_WARN:
                flags.append("high_height_cv")

            # Check R² quality (always)
            if stats.r2_mean < 0.5:
                flags.append("very_poor_r2")
            elif stats.r2_mean < self.R2_POOR:
                flags.append("poor_r2")

            # Check linewidth stability (always)
            if stats.lw_f1_cv > self.LW_CV_WARN or stats.lw_f2_cv > self.LW_CV_WARN:
                flags.append("lw_unstable")

            # Check LW jumps (always)
            if stats.lw_f1_jumps or stats.lw_f2_jumps:
                flags.append("lw_jump")

            # Check position drift (always)
            if stats.pos_f1_drift > self.POS_F1_DRIFT_WARN:
                flags.append("pos_f1_drift")
            if stats.pos_f2_drift > self.POS_F2_DRIFT_WARN:
                flags.append("pos_f2_drift")

            # Check cluster stability (always)
            if not stats.is_cluster_stable:
                flags.append("cluster_changed")

            stats.flags = flags
            if flags:
                flagged.add(pid)

        report.flagged_peaks = sorted(list(flagged))

        self.report = report
        return report

    def save_report(self, output_dir: str) -> str:
        """Save QC_report.json to output directory. Returns path."""
        if self.report is None:
            return ""

        output_path = Path(output_dir) / "QC_report.json"

        # Convert report to dict for JSON serialization
        report_dict = {
            "source_path": self.report.source_path,
            "n_peaks": self.report.n_peaks,
            "n_spectra": self.report.n_spectra,
            "spectrum_names": self.report.spectrum_names,
            "processing_mode": self.report.processing_mode,
            "nucleus_type": self.report.nucleus_type,
            "series_type": self.report.series_type.value,
            "volume_cv_median": self.report.volume_cv_median,
            "volume_cv_mean": self.report.volume_cv_mean,
            "r2_mean": self.report.r2_mean,
            "r2_median": self.report.r2_median,
            "n_cluster_stable": self.report.n_cluster_stable,
            "n_lw_stable": self.report.n_lw_stable,
            "fit_success_rate": self.report.fit_success_rate,
            "flagged_peaks": self.report.flagged_peaks,
            "dummy_peak_ids": self.report.dummy_peak_ids,
            "failed_peak_ids": self.report.failed_peak_ids,
            "peak_stats": {
                pid: {
                    "peak_id": s.peak_id,
                    "assignment": s.assignment,
                    "volume_cv": s.volume_cv,
                    "height_cv": s.height_cv,
                    "r2_mean": s.r2_mean,
                    "r2_min": s.r2_min,
                    "lw_f1_cv": s.lw_f1_cv,
                    "lw_f2_cv": s.lw_f2_cv,
                    "pos_f1_drift": s.pos_f1_drift,
                    "pos_f2_drift": s.pos_f2_drift,
                    "is_cluster_stable": s.is_cluster_stable,
                    "cluster_transitions": s.cluster_transitions,
                    "lw_f1_jumps": s.lw_f1_jumps,
                    "lw_f2_jumps": s.lw_f2_jumps,
                    "flags": s.flags,
                }
                for pid, s in self.report.peak_stats.items()
            },
        }

        with open(output_path, 'w') as f:
            json.dump(report_dict, f, indent=2)

        return str(output_path)


class SeriesComparator:
    """Compares multiple series integration runs."""

    def __init__(self, series_type: SeriesType = SeriesType.REPLICATE):
        self.series: Dict[str, SeriesQCReport] = {}
        self.series_type = series_type

    def add_series(self, name: str, training_data_path: str) -> bool:
        """Load and analyze a series.

        Args:
            name: Name to identify this series (e.g., "param_A")
            training_data_path: Path to training_data.json

        Returns:
            True if successfully loaded, False otherwise
        """
        analyzer = SeriesQCAnalyzer(series_type=self.series_type)
        if not analyzer.load_training_data(training_data_path):
            return False

        analyzer.compute_per_peak_stats()
        analyzer.analyze_cluster_stability()
        analyzer.detect_lw_jumps()
        report = analyzer.generate_report()

        if report is None:
            return False

        self.series[name] = report
        return True

    def get_common_peaks(self) -> List[str]:
        """Get peak IDs present in ALL series (intersection)."""
        if not self.series:
            return []

        # Start with peaks from first series
        series_list = list(self.series.values())
        common = set(series_list[0].peak_stats.keys())

        # Intersect with each subsequent series
        for report in series_list[1:]:
            common &= set(report.peak_stats.keys())

        return sorted(list(common))

    def get_series_only_peaks(self, series_name: str) -> List[str]:
        """Get peaks only in this series (not in others)."""
        if series_name not in self.series:
            return []

        this_series_peaks = set(self.series[series_name].peak_stats.keys())

        # Collect peaks from all other series
        other_peaks = set()
        for name, report in self.series.items():
            if name != series_name:
                other_peaks |= set(report.peak_stats.keys())

        # Return peaks in this series but not in others
        unique = this_series_peaks - other_peaks
        return sorted(list(unique))

    def compare(self) -> Dict:
        """Generate comparison report with winner per metric."""
        if not self.series:
            return {}

        result = {}

        # Compute metrics for each series
        metrics = {
            "volume_cv_median": ("lower_better", lambda r: r.volume_cv_median),
            "r2_mean": ("higher_better", lambda r: r.r2_mean),
            "n_flagged": ("lower_better", lambda r: len(r.flagged_peaks)),
            "n_cluster_stable": ("higher_better", lambda r: r.n_cluster_stable),
        }

        for metric_name, (direction, getter) in metrics.items():
            values = {}
            for name, report in self.series.items():
                values[name] = getter(report)

            # Determine winner
            if direction == "lower_better":
                winner = min(values.keys(), key=lambda k: values[k])
            else:
                winner = max(values.keys(), key=lambda k: values[k])

            result[metric_name] = {
                "values": values,
                "winner": winner,
            }

        # Add summary counts
        result["n_common_peaks"] = len(self.get_common_peaks())
        result["n_unique_peaks"] = {
            name: len(self.get_series_only_peaks(name))
            for name in self.series.keys()
        }

        return result

    def get_per_peak_winner(self) -> Dict[str, str]:
        """For each peak, which series gave best result (by R² mean)."""
        common_peaks = self.get_common_peaks()
        winners = {}

        for pid in common_peaks:
            best_series = None
            best_r2 = -1.0

            for name, report in self.series.items():
                if pid in report.peak_stats:
                    r2 = report.peak_stats[pid].r2_mean
                    if r2 > best_r2:
                        best_r2 = r2
                        best_series = name

            if best_series is not None:
                winners[pid] = best_series

        return winners

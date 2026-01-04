# ABOUTME: Extracts local spectral features at candidate peak locations for ML classification.
# ABOUTME: Used to distinguish real NMR peaks from noise based on shape, intensity, and symmetry.

"""
Peak Feature Extractor for ML-based peak detection.

Extracts 16 features from the local spectral region around a candidate location:
- Intensity features (peak value, local SNR)
- Width features (FWHM estimates in F1 and F2)
- Shape features (symmetry, prominence, sharpness)
- Context features (nearby maxima, curvature, neighborhood density)
- Peak-likeness features (template correlation, Gaussian fit R²)

These features are used to train a Random Forest classifier to distinguish
real peaks from noise.
"""

import numpy as np
from typing import Dict, Tuple, Optional, List
from dataclasses import dataclass
from scipy.ndimage import maximum_filter
from scipy.signal import find_peaks
from scipy.optimize import curve_fit


@dataclass
class PeakFeatures:
    """Container for extracted peak features."""

    # Intensity features
    intensity: float = 0.0
    local_snr: float = 0.0

    # Width features
    fwhm_f1: float = 0.0
    fwhm_f2: float = 0.0
    aspect_ratio: float = 0.0

    # Shape features
    symmetry_f1: float = 0.0
    symmetry_f2: float = 0.0
    prominence_f1: float = 0.0
    prominence_f2: float = 0.0

    # Context features
    local_max_count: int = 0
    gradient_magnitude: float = 0.0
    curvature: float = 0.0

    # Peak-likeness features
    template_correlation: float = 0.0
    gaussian_fit_r2: float = 0.0

    # Additional discriminative features
    peak_sharpness: float = 0.0      # 2nd derivative magnitude at apex
    neighborhood_peaks: float = 0.0  # Density of peaks in local neighborhood

    def to_array(self) -> np.ndarray:
        """Convert features to numpy array for ML model input."""
        return np.array([
            self.intensity,
            self.local_snr,
            self.fwhm_f1,
            self.fwhm_f2,
            self.aspect_ratio,
            self.symmetry_f1,
            self.symmetry_f2,
            self.prominence_f1,
            self.prominence_f2,
            self.local_max_count,
            self.gradient_magnitude,
            self.curvature,
            self.template_correlation,
            self.gaussian_fit_r2,
            self.peak_sharpness,
            self.neighborhood_peaks,
        ])

    @staticmethod
    def feature_names() -> List[str]:
        """Return list of feature names in order."""
        return [
            'intensity',
            'local_snr',
            'fwhm_f1',
            'fwhm_f2',
            'aspect_ratio',
            'symmetry_f1',
            'symmetry_f2',
            'prominence_f1',
            'prominence_f2',
            'local_max_count',
            'gradient_magnitude',
            'curvature',
            'template_correlation',
            'gaussian_fit_r2',
            'peak_sharpness',
            'neighborhood_peaks',
        ]


def gaussian_2d(coords, amplitude, x0, y0, sigma_x, sigma_y, offset):
    """2D Gaussian function for curve_fit."""
    x, y = coords
    return offset + amplitude * np.exp(
        -((x - x0)**2 / (2 * sigma_x**2) + (y - y0)**2 / (2 * sigma_y**2))
    )


class PeakTemplateGenerator:
    """
    Generate nucleus-specific 2D Gaussian templates for peak detection.

    Uses median linewidths from training data to create templates that
    match typical peak shapes for 15N-HSQC and 13C-HSQC spectra.
    """

    # Median linewidths from training data (in ppm)
    # These were computed from ml_training_data/adaptativev1/training_data.json
    LINEWIDTHS = {
        '15N': {
            'lw_f1': 0.224,  # 15N dimension median linewidth
            'lw_f2': 0.035,  # 1H dimension median linewidth
        },
        '13C': {
            'lw_f1': 0.137,  # 13C dimension median linewidth
            'lw_f2': 0.028,  # 1H dimension median linewidth
        },
    }

    def __init__(self, nucleus_type: str = '15N', template_size: int = 7):
        """
        Initialize template generator.

        Args:
            nucleus_type: '15N' or '13C' for appropriate linewidths
            template_size: Size of template in points (default 7x7)
        """
        self.nucleus_type = nucleus_type if nucleus_type in self.LINEWIDTHS else '15N'
        self.template_size = template_size
        self._template_cache = {}  # Cache templates by resolution

    def get_template(
        self,
        ppm_per_point_f1: float,
        ppm_per_point_f2: float,
    ) -> np.ndarray:
        """
        Get or generate 2D Gaussian template at given resolution.

        Args:
            ppm_per_point_f1: PPM per data point in F1 dimension
            ppm_per_point_f2: PPM per data point in F2 dimension

        Returns:
            2D numpy array with normalized template (max = 1)
        """
        # Cache key based on resolution (rounded to avoid floating point issues)
        cache_key = (round(ppm_per_point_f1, 6), round(ppm_per_point_f2, 6))

        if cache_key in self._template_cache:
            return self._template_cache[cache_key]

        # Generate template
        template = self._generate_template(ppm_per_point_f1, ppm_per_point_f2)
        self._template_cache[cache_key] = template
        return template

    def _generate_template(
        self,
        ppm_per_point_f1: float,
        ppm_per_point_f2: float,
    ) -> np.ndarray:
        """Generate 2D Gaussian template."""
        lw = self.LINEWIDTHS[self.nucleus_type]

        # Convert linewidths to sigma (FWHM = 2.355 * sigma)
        sigma_f1_ppm = lw['lw_f1'] / 2.355
        sigma_f2_ppm = lw['lw_f2'] / 2.355

        # Convert to points
        sigma_f1_pts = sigma_f1_ppm / max(abs(ppm_per_point_f1), 1e-10)
        sigma_f2_pts = sigma_f2_ppm / max(abs(ppm_per_point_f2), 1e-10)

        # Ensure minimum sigma
        sigma_f1_pts = max(sigma_f1_pts, 0.5)
        sigma_f2_pts = max(sigma_f2_pts, 0.5)

        # Create coordinate grid centered at (0, 0)
        half = self.template_size // 2
        y = np.arange(-half, half + 1)
        x = np.arange(-half, half + 1)
        xx, yy = np.meshgrid(x, y)

        # Generate 2D Gaussian
        template = np.exp(
            -(xx**2 / (2 * sigma_f2_pts**2) + yy**2 / (2 * sigma_f1_pts**2))
        )

        # Normalize to max = 1
        template = template / template.max()

        return template

    def correlate(self, local_region: np.ndarray, template: np.ndarray) -> float:
        """
        Compute normalized cross-correlation between local region and template.

        Args:
            local_region: Extracted spectrum region around candidate
            template: Pre-computed template from get_template()

        Returns:
            Correlation coefficient in range [-1, 1], typically [0, 1] for peaks
        """
        if local_region.shape != template.shape:
            return 0.0

        # Flatten for correlation
        region_flat = local_region.flatten()
        template_flat = template.flatten()

        # Center the data
        region_centered = region_flat - np.mean(region_flat)
        template_centered = template_flat - np.mean(template_flat)

        # Compute norms
        region_norm = np.linalg.norm(region_centered)
        template_norm = np.linalg.norm(template_centered)

        if region_norm < 1e-10 or template_norm < 1e-10:
            return 0.0

        # Normalized cross-correlation
        correlation = np.dot(region_centered, template_centered) / (region_norm * template_norm)

        return float(correlation)


class PeakFeatureExtractor:
    """
    Extracts features from spectrum at candidate peak locations.

    Features are designed to discriminate between real peaks and noise:
    - Real peaks have consistent FWHM, good symmetry, high prominence
    - Noise has irregular shape, poor symmetry, low prominence
    """

    def __init__(
        self,
        noise_level: float = 1.0,
        window_size_f1: int = 15,
        window_size_f2: int = 15,
        nucleus_type: str = '15N',
        ppm_per_point_f1: float = 0.1,
        ppm_per_point_f2: float = 0.01,
    ):
        """
        Initialize feature extractor.

        Args:
            noise_level: Estimated noise level for SNR calculation
            window_size_f1: Window size (points) in F1 dimension for local analysis
            window_size_f2: Window size (points) in F2 dimension for local analysis
            nucleus_type: '15N' or '13C' for nucleus-specific templates
            ppm_per_point_f1: PPM per data point in F1 dimension
            ppm_per_point_f2: PPM per data point in F2 dimension
        """
        self.noise_level = max(noise_level, 1e-10)
        self.window_size_f1 = window_size_f1
        self.window_size_f2 = window_size_f2
        self.nucleus_type = nucleus_type
        self.ppm_per_point_f1 = ppm_per_point_f1
        self.ppm_per_point_f2 = ppm_per_point_f2

        # Initialize template generator for peak-likeness features
        self.template_generator = PeakTemplateGenerator(
            nucleus_type=nucleus_type,
            template_size=7,
        )
        self._template = None  # Lazy-loaded template

    def extract_features(
        self,
        spectrum: np.ndarray,
        y_idx: int,
        x_idx: int,
    ) -> PeakFeatures:
        """
        Extract features at a candidate location.

        Args:
            spectrum: 2D spectrum data (y=F1, x=F2)
            y_idx: Row index (F1 dimension)
            x_idx: Column index (F2 dimension)

        Returns:
            PeakFeatures object with all extracted features
        """
        features = PeakFeatures()

        # Validate indices
        h, w = spectrum.shape
        if not (0 <= y_idx < h and 0 <= x_idx < w):
            return features

        # Extract 1D slices through candidate
        slice_f1 = spectrum[:, x_idx]  # Column (F1 dimension)
        slice_f2 = spectrum[y_idx, :]  # Row (F2 dimension)

        # Intensity features
        features.intensity = float(spectrum[y_idx, x_idx])
        features.local_snr = self._compute_local_snr(spectrum, y_idx, x_idx)

        # Width features (FWHM)
        features.fwhm_f1 = self._estimate_fwhm(slice_f1, y_idx)
        features.fwhm_f2 = self._estimate_fwhm(slice_f2, x_idx)
        if features.fwhm_f2 > 0:
            features.aspect_ratio = features.fwhm_f1 / features.fwhm_f2
        else:
            features.aspect_ratio = 1.0

        # Shape features
        features.symmetry_f1 = self._measure_symmetry(slice_f1, y_idx)
        features.symmetry_f2 = self._measure_symmetry(slice_f2, x_idx)
        features.prominence_f1 = self._compute_prominence(slice_f1, y_idx)
        features.prominence_f2 = self._compute_prominence(slice_f2, x_idx)

        # Context features
        features.local_max_count = self._count_local_maxima(spectrum, y_idx, x_idx)
        features.gradient_magnitude = self._compute_gradient(spectrum, y_idx, x_idx)
        features.curvature = self._compute_curvature(spectrum, y_idx, x_idx)

        # Peak-likeness features
        features.template_correlation = self._compute_template_correlation(spectrum, y_idx, x_idx)
        features.gaussian_fit_r2 = self._compute_gaussian_fit_r2(spectrum, y_idx, x_idx)

        # Additional discriminative features
        features.peak_sharpness = self._compute_peak_sharpness(spectrum, y_idx, x_idx)
        features.neighborhood_peaks = self._compute_neighborhood_peaks(spectrum, y_idx, x_idx)

        return features

    def extract_features_batch(
        self,
        spectrum: np.ndarray,
        candidates: List[Dict],
    ) -> np.ndarray:
        """
        Extract features for multiple candidates.

        Args:
            spectrum: 2D spectrum data
            candidates: List of dicts with 'y_idx' and 'x_idx' keys

        Returns:
            2D array of shape (n_candidates, n_features)
        """
        features_list = []
        for candidate in candidates:
            y_idx = candidate.get('y_idx', candidate.get('row', 0))
            x_idx = candidate.get('x_idx', candidate.get('col', 0))
            features = self.extract_features(spectrum, y_idx, x_idx)
            features_list.append(features.to_array())

        return np.array(features_list)

    def _compute_local_snr(
        self,
        spectrum: np.ndarray,
        y_idx: int,
        x_idx: int,
    ) -> float:
        """Compute local SNR at candidate location."""
        intensity = spectrum[y_idx, x_idx]

        # Also compute local noise from surrounding region
        h, w = spectrum.shape
        y_start = max(0, y_idx - self.window_size_f1)
        y_end = min(h, y_idx + self.window_size_f1 + 1)
        x_start = max(0, x_idx - self.window_size_f2)
        x_end = min(w, x_idx + self.window_size_f2 + 1)

        local_region = spectrum[y_start:y_end, x_start:x_end]
        local_noise = np.std(local_region)

        # Use the larger of local noise and global noise estimate
        effective_noise = max(local_noise, self.noise_level, 1e-10)

        return float(intensity / effective_noise)

    def _estimate_fwhm(self, profile: np.ndarray, center_idx: int) -> float:
        """
        Estimate FWHM from 1D profile.

        Uses half-maximum crossing points on both sides of peak.
        """
        if center_idx < 0 or center_idx >= len(profile):
            return 0.0

        peak_value = profile[center_idx]
        if peak_value <= 0:
            return 0.0

        half_max = peak_value / 2.0

        # Find left crossing
        left_idx = center_idx
        for i in range(center_idx - 1, -1, -1):
            if profile[i] < half_max:
                # Interpolate
                if profile[i + 1] != profile[i]:
                    frac = (half_max - profile[i]) / (profile[i + 1] - profile[i])
                    left_idx = i + frac
                else:
                    left_idx = i + 0.5
                break
        else:
            left_idx = 0

        # Find right crossing
        right_idx = center_idx
        for i in range(center_idx + 1, len(profile)):
            if profile[i] < half_max:
                # Interpolate
                if profile[i - 1] != profile[i]:
                    frac = (half_max - profile[i]) / (profile[i - 1] - profile[i])
                    right_idx = i - frac
                else:
                    right_idx = i - 0.5
                break
        else:
            right_idx = len(profile) - 1

        return float(right_idx - left_idx)

    def _measure_symmetry(self, profile: np.ndarray, center_idx: int) -> float:
        """
        Measure symmetry of profile around center.

        Returns value between 0 (asymmetric) and 1 (perfectly symmetric).
        Uses correlation between left and right sides.
        """
        n = len(profile)
        if center_idx < 1 or center_idx >= n - 1:
            return 0.0

        # Determine window size (limited by distance to edge)
        max_dist = min(center_idx, n - 1 - center_idx)
        window = min(max_dist, max(self.window_size_f1, self.window_size_f2) // 2)

        if window < 2:
            return 0.0

        # Extract left and right sides
        left = profile[center_idx - window:center_idx]
        right = profile[center_idx + 1:center_idx + 1 + window][::-1]  # Reverse right side

        # Ensure same length
        min_len = min(len(left), len(right))
        if min_len < 2:
            return 0.0
        left = left[-min_len:]
        right = right[:min_len]

        # Compute correlation
        left_centered = left - np.mean(left)
        right_centered = right - np.mean(right)

        left_norm = np.linalg.norm(left_centered)
        right_norm = np.linalg.norm(right_centered)

        if left_norm < 1e-10 or right_norm < 1e-10:
            return 0.0

        correlation = np.dot(left_centered, right_centered) / (left_norm * right_norm)

        # Convert to 0-1 range (correlation can be -1 to 1)
        return float((correlation + 1) / 2)

    def _compute_prominence(self, profile: np.ndarray, center_idx: int) -> float:
        """
        Compute peak prominence in 1D profile.

        Prominence is the height of peak above the highest connecting saddle point.
        """
        if center_idx < 0 or center_idx >= len(profile):
            return 0.0

        peak_value = profile[center_idx]

        # Find minimum on left side
        left_min = peak_value
        for i in range(center_idx - 1, -1, -1):
            if profile[i] < left_min:
                left_min = profile[i]
            if profile[i] > peak_value:
                break

        # Find minimum on right side
        right_min = peak_value
        for i in range(center_idx + 1, len(profile)):
            if profile[i] < right_min:
                right_min = profile[i]
            if profile[i] > peak_value:
                break

        # Prominence is height above the higher of the two saddle points
        saddle = max(left_min, right_min)
        prominence = peak_value - saddle

        # Normalize by peak value
        if peak_value > 0:
            return float(prominence / peak_value)
        return 0.0

    def _count_local_maxima(
        self,
        spectrum: np.ndarray,
        y_idx: int,
        x_idx: int,
        radius: int = 5,
    ) -> int:
        """Count local maxima within radius of candidate."""
        h, w = spectrum.shape

        y_start = max(0, y_idx - radius)
        y_end = min(h, y_idx + radius + 1)
        x_start = max(0, x_idx - radius)
        x_end = min(w, x_idx + radius + 1)

        local_region = spectrum[y_start:y_end, x_start:x_end]

        # Find local maxima using maximum filter
        max_filtered = maximum_filter(local_region, size=3)
        is_max = (local_region == max_filtered) & (local_region > 0)

        return int(np.sum(is_max))

    def _compute_gradient(
        self,
        spectrum: np.ndarray,
        y_idx: int,
        x_idx: int,
    ) -> float:
        """Compute gradient magnitude at candidate location."""
        h, w = spectrum.shape

        # Compute gradients using central differences
        if y_idx > 0 and y_idx < h - 1:
            grad_y = (spectrum[y_idx + 1, x_idx] - spectrum[y_idx - 1, x_idx]) / 2
        else:
            grad_y = 0.0

        if x_idx > 0 and x_idx < w - 1:
            grad_x = (spectrum[y_idx, x_idx + 1] - spectrum[y_idx, x_idx - 1]) / 2
        else:
            grad_x = 0.0

        return float(np.sqrt(grad_y**2 + grad_x**2))

    def _compute_curvature(
        self,
        spectrum: np.ndarray,
        y_idx: int,
        x_idx: int,
    ) -> float:
        """
        Compute curvature (Laplacian) at candidate location.

        Peaks have negative curvature (concave down), noise is variable.
        """
        h, w = spectrum.shape

        # Compute second derivatives
        if y_idx > 0 and y_idx < h - 1:
            d2y = (spectrum[y_idx + 1, x_idx] - 2 * spectrum[y_idx, x_idx] +
                   spectrum[y_idx - 1, x_idx])
        else:
            d2y = 0.0

        if x_idx > 0 and x_idx < w - 1:
            d2x = (spectrum[y_idx, x_idx + 1] - 2 * spectrum[y_idx, x_idx] +
                   spectrum[y_idx, x_idx - 1])
        else:
            d2x = 0.0

        # Laplacian (sum of second derivatives)
        laplacian = d2y + d2x

        # Normalize by intensity to make scale-invariant
        intensity = spectrum[y_idx, x_idx]
        if abs(intensity) > 1e-10:
            return float(laplacian / intensity)
        return 0.0

    def _compute_template_correlation(
        self,
        spectrum: np.ndarray,
        y_idx: int,
        x_idx: int,
    ) -> float:
        """
        Compute correlation between local region and nucleus-specific template.

        Uses a 7x7 region around the candidate and correlates with a pre-computed
        2D Gaussian template matching typical peak shapes for the nucleus type.

        Returns value in range [-1, 1], typically [0.5, 1.0] for real peaks.
        """
        h, w = spectrum.shape
        half = 3  # 7x7 template

        # Check bounds
        if (y_idx - half < 0 or y_idx + half >= h or
            x_idx - half < 0 or x_idx + half >= w):
            return 0.0

        # Extract local region
        local_region = spectrum[y_idx - half:y_idx + half + 1,
                                x_idx - half:x_idx + half + 1]

        # Get or generate template (lazy load)
        if self._template is None:
            self._template = self.template_generator.get_template(
                self.ppm_per_point_f1,
                self.ppm_per_point_f2,
            )

        # Compute correlation
        correlation = self.template_generator.correlate(local_region, self._template)

        return correlation

    def _compute_gaussian_fit_r2(
        self,
        spectrum: np.ndarray,
        y_idx: int,
        x_idx: int,
    ) -> float:
        """
        Fit 2D Gaussian to local region and return R² goodness-of-fit.

        Real peaks should fit well (R² > 0.8), while noise fits poorly.

        Returns value in range [0, 1], with higher values for better fits.
        """
        h, w = spectrum.shape
        half = 3  # 7x7 region

        # Check bounds
        if (y_idx - half < 0 or y_idx + half >= h or
            x_idx - half < 0 or x_idx + half >= w):
            return 0.0

        # Extract local region
        local_region = spectrum[y_idx - half:y_idx + half + 1,
                                x_idx - half:x_idx + half + 1]

        # Create coordinate grid
        size = 2 * half + 1
        y_coords = np.arange(size)
        x_coords = np.arange(size)
        xx, yy = np.meshgrid(x_coords, y_coords)
        coords = (xx.flatten(), yy.flatten())
        data = local_region.flatten()

        # Initial guesses
        amplitude = local_region[half, half] - local_region.min()
        x0, y0 = half, half
        sigma_x, sigma_y = 1.5, 1.5
        offset = local_region.min()

        if amplitude <= 0:
            return 0.0

        try:
            # Fit 2D Gaussian with bounds
            popt, _ = curve_fit(
                gaussian_2d,
                coords,
                data,
                p0=[amplitude, x0, y0, sigma_x, sigma_y, offset],
                bounds=(
                    [0, 0, 0, 0.3, 0.3, -np.inf],  # Lower bounds
                    [np.inf, size, size, size, size, np.inf],  # Upper bounds
                ),
                maxfev=200,
            )

            # Compute R²
            fitted = gaussian_2d(coords, *popt)
            ss_res = np.sum((data - fitted) ** 2)
            ss_tot = np.sum((data - np.mean(data)) ** 2)

            if ss_tot < 1e-10:
                return 0.0

            r_squared = 1 - (ss_res / ss_tot)
            return float(max(0.0, r_squared))

        except (RuntimeError, ValueError):
            # Fit failed
            return 0.0

    def _compute_peak_sharpness(
        self,
        spectrum: np.ndarray,
        y_idx: int,
        x_idx: int,
    ) -> float:
        """
        Compute peak sharpness using 2nd derivative magnitude at apex.

        Real peaks have characteristic sharpness (negative 2nd derivative),
        while noise tends to have irregular/inconsistent sharpness.

        Returns normalized sharpness value (higher = sharper peak).
        """
        h, w = spectrum.shape

        # Need at least 2 neighbors in each direction
        if y_idx < 2 or y_idx >= h - 2 or x_idx < 2 or x_idx >= w - 2:
            return 0.0

        intensity = spectrum[y_idx, x_idx]
        if intensity <= 0:
            return 0.0

        # Compute 2nd derivative in both dimensions using 5-point stencil
        # More accurate than 3-point: f''(x) ≈ (-f(x-2) + 16*f(x-1) - 30*f(x) + 16*f(x+1) - f(x+2)) / 12
        d2y = (-spectrum[y_idx - 2, x_idx] +
               16 * spectrum[y_idx - 1, x_idx] -
               30 * spectrum[y_idx, x_idx] +
               16 * spectrum[y_idx + 1, x_idx] -
               spectrum[y_idx + 2, x_idx]) / 12.0

        d2x = (-spectrum[y_idx, x_idx - 2] +
               16 * spectrum[y_idx, x_idx - 1] -
               30 * spectrum[y_idx, x_idx] +
               16 * spectrum[y_idx, x_idx + 1] -
               spectrum[y_idx, x_idx + 2]) / 12.0

        # Peak sharpness is the magnitude of the negative 2nd derivative
        # Real peaks have negative 2nd derivative (concave down at apex)
        sharpness = -1.0 * (d2y + d2x)

        # Normalize by intensity to be scale-invariant
        normalized_sharpness = sharpness / intensity

        # Clip to reasonable range and scale to [0, 1]
        # Typical values for real peaks are 0.1-2.0 after normalization
        return float(max(0.0, min(1.0, normalized_sharpness / 2.0)))

    def _compute_neighborhood_peaks(
        self,
        spectrum: np.ndarray,
        y_idx: int,
        x_idx: int,
        inner_radius: int = 3,
        outer_radius: int = 15,
    ) -> float:
        """
        Compute density of peaks in the local neighborhood.

        Real peaks in NMR spectra tend to appear in regions with other peaks
        (protein spectra have clustered peaks), while isolated noise is random.

        Args:
            inner_radius: Exclude this region around candidate (don't count self)
            outer_radius: Search radius for neighborhood peaks

        Returns:
            Normalized peak density in range [0, 1]
        """
        h, w = spectrum.shape

        # Define outer search region
        y_start = max(0, y_idx - outer_radius)
        y_end = min(h, y_idx + outer_radius + 1)
        x_start = max(0, x_idx - outer_radius)
        x_end = min(w, x_idx + outer_radius + 1)

        # Extract region
        region = spectrum[y_start:y_end, x_start:x_end].copy()

        # Create mask to exclude inner region (the candidate itself)
        local_y = y_idx - y_start
        local_x = x_idx - x_start

        # Zero out inner region so we don't count the candidate
        inner_y_start = max(0, local_y - inner_radius)
        inner_y_end = min(region.shape[0], local_y + inner_radius + 1)
        inner_x_start = max(0, local_x - inner_radius)
        inner_x_end = min(region.shape[1], local_x + inner_radius + 1)
        region[inner_y_start:inner_y_end, inner_x_start:inner_x_end] = 0

        # Find local maxima in the neighborhood
        max_filtered = maximum_filter(region, size=3)
        is_max = (region == max_filtered) & (region > 0)

        # Count peaks above noise threshold
        threshold = self.noise_level * 5  # Only count significant peaks
        is_significant = region > threshold
        neighbor_peaks = np.sum(is_max & is_significant)

        # Normalize by area (accounting for excluded region)
        total_area = region.size
        inner_area = (inner_y_end - inner_y_start) * (inner_x_end - inner_x_start)
        effective_area = total_area - inner_area

        if effective_area <= 0:
            return 0.0

        # Density: peaks per 100 points, capped at 1
        density = (neighbor_peaks / effective_area) * 100
        return float(min(1.0, density))


def ppm_to_index(ppm_axis: np.ndarray, ppm_value: float) -> int:
    """Convert ppm value to nearest index."""
    return int(np.argmin(np.abs(ppm_axis - ppm_value)))


def extract_features_at_ppm(
    spectrum: np.ndarray,
    ppm_x_axis: np.ndarray,
    ppm_y_axis: np.ndarray,
    pos_f1: float,
    pos_f2: float,
    noise_level: float = 1.0,
    nucleus_type: str = '15N',
) -> PeakFeatures:
    """
    Extract features at a position specified in ppm.

    Args:
        spectrum: 2D spectrum data
        ppm_x_axis: PPM values for x-axis (F2)
        ppm_y_axis: PPM values for y-axis (F1)
        pos_f1: Position in F1 dimension (ppm)
        pos_f2: Position in F2 dimension (ppm)
        noise_level: Estimated noise level
        nucleus_type: '15N' or '13C' for nucleus-specific templates

    Returns:
        PeakFeatures object
    """
    y_idx = ppm_to_index(ppm_y_axis, pos_f1)
    x_idx = ppm_to_index(ppm_x_axis, pos_f2)

    # Compute ppm per point from axes
    ppm_per_point_f1 = abs(ppm_y_axis[1] - ppm_y_axis[0]) if len(ppm_y_axis) > 1 else 0.1
    ppm_per_point_f2 = abs(ppm_x_axis[1] - ppm_x_axis[0]) if len(ppm_x_axis) > 1 else 0.01

    extractor = PeakFeatureExtractor(
        noise_level=noise_level,
        nucleus_type=nucleus_type,
        ppm_per_point_f1=ppm_per_point_f1,
        ppm_per_point_f2=ppm_per_point_f2,
    )
    return extractor.extract_features(spectrum, y_idx, x_idx)

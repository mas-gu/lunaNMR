#!/usr/bin/env python3
# ABOUTME: Voigt fitter held by the integrator: 1D Voigt profile, R², nucleus detection, local regions,
# ABOUTME: and enhanced_peak_fitting_parallel, the entry point that fits a peak list cluster by cluster in 2D.

import numpy as np
from scipy.special import wofz

from lunaNMR.utils.output_manager import log_progress, log_warning


class EnhancedVoigtFitter:
    """Voigt profile helpers and the cluster-by-cluster 2D fit of a peak list"""

    def __init__(self):
        self.fitting_parameters = {
            'max_iterations': 1000,
            'tolerance': 1e-8,
            'min_r_squared': 0.3,  # Lowered for realistic data
            'baseline_polynomial_degree': 1,
            'smoothing_window': 5,
            'peak_detection_prominence': 0.1,
            ##
            'multi_peak_r2_threshold': 0.7,          # R² threshold to trigger multi-peak detection
            'multi_peak_improvement_threshold': 0.1,  # Minimum improvement required from multi-peak fit
            'peak_detection_sensitivity': 1.5,        # Height threshold multiplier (lower = more sensitive)
            'overlap_detection_factor': 0.8,          # Separation factor for overlap detection
            'residual_analysis_threshold': 1.5      # Threshold for residual peak detection

        }

        # NMR-specific parameters (corrected ranges)
        self.nmr_ranges = {
            '1H': {'min': 5.5, 'max': 12.0, 'typical_width': 0.01},  #0.02    # ¹H range
            '15N': {'min': 100.0, 'max': 140.0, 'typical_width': 0.5}, #1.5  # ¹⁵N range
            '13C': {'min': 5, 'max': 50.0, 'typical_width': 0.6}  #1.0    # ¹³C range for completeness
        }

    def _calculate_gui_based_multiplier(self, nucleus_type, ppm_range, _data_length, fitted_width=None):
        """
        DEPRECATED: This function is NO LONGER USED.


        FIXED windows (NO multipliers, NO adaptation):
        - radF1 = 0.6 ppm (indirect dimension, 15N/13C)
        - radF2 = 0.06 ppm (direct dimension, 1H)

        Data selection uses elliptical boundary: radius² = (F1-pos)²/radF1² + (F2-pos)²/radF2² <= 1.0

        This function returns a dummy value for backward compatibility only.
        The actual window selection is handled by Ps2dDataSelector in core_integrator.py.

        Parameters:
        -----------
        nucleus_type : str
            '1H', '15N', or '13C' - determines which GUI parameter to use
        ppm_range : float
            Total ppm range of the data
        data_length : int
            Number of data points in the spectrum dimension
        fitted_width : float, optional
            Fitted peak linewidth for reference

        Returns:
        --------
        float : Returns 1.0 (dummy value - NOT USED)
        """
        # RETURN DUMMY VALUE - NOT USED ANYMORE
        return 1.0  # No multiplier needed - PS2D uses fixed radF1/radF2

    @staticmethod
    def voigt_profile_1d(x, amplitude, center, sigma, gamma, baseline=0):
        """
        1D Voigt profile using Faddeeva function

        Parameters:
        - x: frequency/chemical shift array
        - amplitude: peak amplitude
        - center: peak center position
        - sigma: Gaussian width (instrumental broadening)
        - gamma: Lorentzian width (natural line width)
        - baseline: baseline offset
        """
        try:
            # Avoid division by zero
            sigma = max(sigma, 1e-6)

            # Compute complex argument for Faddeeva function
            z = ((x - center) + 1j*gamma) / (sigma * np.sqrt(2))

            # Voigt profile using Faddeeva function
            voigt = amplitude * np.real(wofz(z)) / (sigma * np.sqrt(2*np.pi))

            return voigt + baseline
        except Exception:
            # Fallback to Gaussian if Voigt fails
            return amplitude * np.exp(-0.5 * ((x - center) / max(sigma, 1e-6))**2) + baseline

    @staticmethod
    def calculate_r_squared(y_true, y_pred):
        """Calculate R-squared value (moved here for better organization)"""
        ss_res = np.sum((y_true - y_pred) ** 2)
        ss_tot = np.sum((y_true - np.mean(y_true)) ** 2)
        return 1 - (ss_res / ss_tot) if ss_tot != 0 else 0

    def detect_nucleus_type(self, ppm_range):
        """
        Detect nucleus type based on ppm range
        """
        ppm_span = abs(ppm_range[1] - ppm_range[0])
        center_ppm = (ppm_range[0] + ppm_range[1]) / 2

        # Check ¹H range (5.5-12 ppm)
        if 5.0 <= center_ppm <= 13.0 and ppm_span < 20:
            return '1H'

        # Check ¹⁵N range (100-140 ppm)
        elif 90 <= center_ppm <= 150 and ppm_span < 100:
            return '15N'

        # Check ¹³C range (0-220 ppm)
        elif 5 <= center_ppm <= 50 and ppm_span < 50: #was 230
            return '13C'

        # Default to ¹H for small ranges
        return '1H'

    def adaptive_width_estimation(self, x_data, y_data, peak_center, nucleus_type=None):
        """
        Adaptive width estimation using multiple approaches
        """
        if nucleus_type is None:
            nucleus_type = self.detect_nucleus_type([x_data[0], x_data[-1]])

        typical_width = self.nmr_ranges.get(nucleus_type, self.nmr_ranges['1H'])['typical_width']

        try:
            # Method 1: FWHM-based estimation
            baseline = self.robust_baseline_estimation(x_data, y_data)
            peak_height = np.max(y_data) - baseline
            half_max = baseline + peak_height / 2

            above_half = y_data > half_max
            if np.sum(above_half) > 3:  # Need at least 3 points
                indices = np.where(above_half)[0]
                fwhm = x_data[indices[-1]] - x_data[indices[0]]
                width_fwhm = fwhm / 2.355  # Convert FWHM to sigma equivalent
            else:
                width_fwhm = typical_width

            # Method 2: Second moment estimation
            try:
                y_norm = y_data - baseline
                y_norm = np.maximum(y_norm, 0)  # Only positive values
                total_intensity = np.sum(y_norm)

                if total_intensity > 0:
                    # Calculate second moment
                    mean_x = np.sum(x_data * y_norm) / total_intensity
                    second_moment = np.sum(y_norm * (x_data - mean_x)**2) / total_intensity
                    width_moment = np.sqrt(second_moment)
                else:
                    width_moment = typical_width
            except Exception:
                width_moment = typical_width

            # Method 3: Gradient-based estimation
            try:
                # Find steepest gradients (peak edges)
                gradient = np.gradient(y_data, x_data)
                max_grad_left = np.argmax(gradient[:len(gradient)//2])
                max_grad_right = np.argmin(gradient[len(gradient)//2:]) + len(gradient)//2

                if max_grad_right > max_grad_left:
                    width_gradient = (x_data[max_grad_right] - x_data[max_grad_left]) / 4
                else:
                    width_gradient = typical_width
            except Exception:
                width_gradient = typical_width

            weights = [0.5, 0.3, 0.2]  # Prefer FWHM, then moment, then gradient
            combined_width = (weights[0] * width_fwhm +
                            weights[1] * width_moment +
                            weights[2] * width_gradient)

            min_width = typical_width * 0.1
            max_width = typical_width * 10

            final_width = np.clip(combined_width, min_width, max_width)

            return final_width

        except Exception:
            return typical_width

    def extract_local_peak_region(self, x_data, y_data, peak_center, fitted_width=None,
                                  nucleus_type=None, window_multiplier=None):
        """
        Extract local region around peak for quality assessment

        ENHANCED: Now supports GUI-based window sizing for consistent display

        This limits quality evaluation to the peak itself within reasonable range
        of its linewidth, excluding distant peaks that shouldn't affect local fit quality.

        Parameters:
        - x_data: full x-axis data
        - y_data: full intensity data
        - peak_center: fitted peak center
        - fitted_width: fitted peak width (sigma + gamma), or None for estimation
        - nucleus_type: nucleus type for typical width fallback
        - window_multiplier: how many widths to include on each side (None = use GUI params)

        Returns:
        - Dictionary with local region data and indices

        BACKWARD COMPATIBILITY:
        - If window_multiplier is provided: uses that value (original behavior)
        - If window_multiplier is None: calculates from GUI parameters (new behavior)
        """
        try:
            if fitted_width is None:
                fitted_width = self.adaptive_width_estimation(x_data, y_data, peak_center, nucleus_type)

            # ENHANCED: Calculate window multiplier from GUI parameters if not provided
            if window_multiplier is None:
                ppm_range = abs(x_data[-1] - x_data[0]) if len(x_data) > 1 else 1.0
                window_multiplier = self._calculate_gui_based_multiplier(
                    nucleus_type, ppm_range, len(x_data), fitted_width
                )

            # Define local window around peak (now using GUI-based multiplier)
            half_window = fitted_width * window_multiplier

            local_mask = (x_data >= peak_center - half_window) & (x_data <= peak_center + half_window)
            local_indices = np.where(local_mask)[0]

            if len(local_indices) < 5:  # Need minimum points for meaningful assessment
                # Expand to at least 5 points around peak center
                center_idx = np.argmin(np.abs(x_data - peak_center))
                min_points = 5
                start_idx = max(0, center_idx - min_points//2)
                end_idx = min(len(x_data), center_idx + min_points//2 + 1)
                local_indices = np.arange(start_idx, end_idx)
                local_mask = np.zeros(len(x_data), dtype=bool)
                local_mask[local_indices] = True

            x_local = x_data[local_indices]
            y_local = y_data[local_indices]

            # Validation warning
            if len(local_indices) < 3:
                log_warning(f"Very small local region ({len(local_indices)} points) for peak at {peak_center:.3f}")

            return {
                #'x_data': x_local,
                #'y_data': y_local,
                #'indices': local_indices,
                'x_local': x_local,
                'y_local': y_local,
                'local_indices': local_indices,
                'peak_center': peak_center,
                'window_size': half_window * 2,  # total window size
                'window_width': half_window * 2,
                'n_points': len(local_indices),
                'ppm_range': x_local[-1] - x_local[0] if len(x_local) > 1 else 0,
                'multiplier_used': window_multiplier,
                'gui_based': window_multiplier != 6.0  # flag for logging/debugging
            }

        except Exception as e:
            log_warning(f"Local region extraction failed: {e}")
            # Fallback to full data
            return {
                'x_local': x_data,
                'y_local': y_data,
                'local_indices': np.arange(len(x_data)),
                'local_mask': np.ones(len(x_data), dtype=bool),
                'window_width': abs(x_data[-1] - x_data[0]),
                'n_points': len(x_data)
            }

    def calculate_local_baseline(self, x_local, y_local):
        """
        Calculate local baseline within the peak region for better quality assessment
        """
        try:
            # For local regions, use edge points for baseline
            n_edge = max(1, len(y_local) // 6)  # Use outer 1/6 on each side

            if len(y_local) < 6:
                # For very small regions, just use minimum
                return np.min(y_local)

            left_edge = y_local[:n_edge]
            right_edge = y_local[-n_edge:]

            # Combine edges and use median for robustness
            edge_values = np.concatenate([left_edge, right_edge])
            local_baseline = np.median(edge_values)

            return local_baseline

        except Exception:
            return np.percentile(y_local, 20)

    def enhanced_peak_fitting_parallel(self, peak_list, progress_callback=None,
                                        parent_integrator=None, locked_clusters_by_assignment=None,
                                        pre_learned_statistics=None, reference_linewidths=None,
                                        skip_series_params=False, reference_positions=None,
                                        in_process=False, stop=None):
        """
        Fit a peak list with the cluster-based parallel workflow.

        Whether a list is worth parallelising is the caller's decision, so a list of any
        size arriving here is fitted in parallel.

        Args:
            peak_list: DataFrame with peak information or single peak coordinates
            progress_callback: Progress update callback function
            parent_integrator: Parent integrator with nmr_data (auto-detected if None)
            locked_clusters_by_assignment: Optional pre-computed clusters from reference spectrum
            pre_learned_statistics: Optional linewidth statistics from reference spectrum.
                If provided, PASS 1 learning is skipped and these values are used as initial guesses.
            reference_linewidths: Optional dict mapping assignment -> {lw_lor_f1, lw_gau_f1, etc.}
                for per-peak linewidth reuse. If provided, linewidths are fixed to these values.
            skip_series_params: If True, don't apply series_params from previous spectra.
                This forces fresh adaptive optimization (used for Independent mode).
            reference_positions: Optional dict mapping assignment -> (x_ppm, y_ppm), the
                original reference positions a cascade's absolute drift limit holds to.
            in_process: Run the same pipeline in this process rather than over a pool -
                sequential mode.
            stop: Optional threading.Event that ends the fit between clusters; a stopped
                fit raises lunaNMR.utils.cancellation.Cancelled.

        Returns:
            Tuple of (fitted_results, learned_statistics) where learned_statistics is the
            spectrum statistics learned during PASS 1 (None if skipped or sequential processing)
        """
        import pandas as pd

        if not isinstance(peak_list, pd.DataFrame):
            # Assume single peak with (peak_x, peak_y, assignment) format
            if isinstance(peak_list, (list, tuple)) and len(peak_list) >= 2:
                peak_x, peak_y = peak_list[0], peak_list[1]
                assignment = peak_list[2] if len(peak_list) > 2 else 'Single_Peak'

                peak_list = pd.DataFrame({
                    'Position_X': [peak_x],
                    'Position_Y': [peak_y],
                    'Assignment': [assignment]
                })
            else:
                raise ValueError("Invalid peak_list format. Expected DataFrame or (peak_x, peak_y, assignment) tuple")

        # Auto-detect parent integrator if not provided
        if parent_integrator is None:
            parent_integrator = getattr(self, 'parent', None)

        if parent_integrator is None:
            # Look for integrator in common places
            import inspect
            frame = inspect.currentframe()
            try:
                # Check calling context for integrator
                while frame:
                    frame_locals = frame.f_locals
                    if 'self' in frame_locals:
                        candidate = frame_locals['self']
                        if hasattr(candidate, 'nmr_data') and hasattr(candidate, 'enhanced_fitter'):
                            if candidate.enhanced_fitter is self:
                                parent_integrator = candidate
                                break
                    frame = frame.f_back
            finally:
                del frame

        if parent_integrator is None or not hasattr(parent_integrator, 'nmr_data'):
            raise ValueError("No parent integrator with nmr_data found - cannot run parallel processing")

        # Thread-safe setup of data context
        import threading
        with threading.Lock():
            self.nmr_data = parent_integrator.nmr_data
            self.ppm_x_axis = parent_integrator.ppm_x_axis
            self.ppm_y_axis = parent_integrator.ppm_y_axis
            # The workers read the noise from here. Without it every parallel fit judged
            # each peak's sign, and its detection sigma, against no noise at all.
            self.noise_level = getattr(parent_integrator, 'noise_level', None)

        from lunaNMR.core.parallel_voigt_processor import ParallelVoigtProcessor

        # Off unless asked for: on real titrations the search always chose the smallest
        # window (1.0 x linewidth), which cuts the peak at half height and broke CSP
        enable_adaptive = False
        if parent_integrator and hasattr(parent_integrator, 'gui_params'):
            enable_adaptive = parent_integrator.gui_params.get('use_adaptive_optimization', False)
        log_progress(f"Using parallel Voigt fitting for {len(peak_list)} peaks")
        parallel_processor = ParallelVoigtProcessor(self, enable_adaptive=enable_adaptive,
                                                    in_process=in_process, stop=stop)
        self.parallel_processor = parallel_processor

        # Pass series_params if available (for subsequent spectra in series)
        # Skip if skip_series_params=True (Independent mode wants fresh adaptive each spectrum)
        if not skip_series_params and hasattr(self, 'series_params') and self.series_params is not None:
            parallel_processor.set_series_params(self.series_params)

        # Pass locked clusters, pre-learned statistics, and reference linewidths
        results, learned_statistics = parallel_processor.fit_all_peaks_parallel(
            peak_list,
            progress_callback,
            locked_clusters_by_assignment=locked_clusters_by_assignment,
            pre_learned_statistics=pre_learned_statistics,
            reference_linewidths=reference_linewidths,
            reference_positions=reference_positions
        )

        if parallel_processor.series_params and not hasattr(self, 'series_params'):
            self.series_params = parallel_processor.series_params
        elif parallel_processor.series_params:
            self.series_params = parallel_processor.series_params

        if len(results) == 1 and len(peak_list) == 1:
            return (results[0], learned_statistics)
        else:
            return (results, learned_statistics)


#!/usr/bin/env python3
# ABOUTME: Voigt fitter held by the integrator: its fitting parameters and enhanced_peak_fitting_parallel,
# ABOUTME: the entry point that fits a peak list cluster by cluster in 2D through the parallel processor.

from lunaNMR.utils.output_manager import log_progress


class EnhancedVoigtFitter:
    """Fitting parameters and the cluster-by-cluster 2D fit of a peak list"""

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


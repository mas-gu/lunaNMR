"""
Single Spectrum Peak Processor
Handles batch fitting of peaks in a single NMR spectrum independently

Author: Guillaume Mas
Date: 2025
"""

from typing import List, Dict, Any, Optional, Tuple
import pandas as pd

from lunaNMR.core.ps2d_config import restore_base_ps2d_config
from lunaNMR.utils.output_manager import log_progress, log_info, log_warning, log_error
from lunaNMR.utils.cancellation import Cancelled

# Below this many peaks, distributing the work costs more than it saves.
MIN_PEAKS_FOR_PARALLEL = 3

# ML module - imported lazily to avoid circular imports and handle missing sklearn
_ml_manager = None


def _get_ml_manager():
    """Lazy initialization of ML manager singleton."""
    global _ml_manager
    if _ml_manager is None:
        try:
            from lunaNMR.ml import ModelManager
            _ml_manager = ModelManager()
        except ImportError:
            _ml_manager = False  # Marker that import failed
    return _ml_manager if _ml_manager is not False else None

class SingleSpectrumProcessor:
    """
    Processes peaks in a single NMR spectrum

    This class encapsulates all logic for fitting multiple peaks in a single
    spectrum, including sequential, parallel, and global optimization modes.
    It operates independently of the GUI and can be easily tested.
    """

    def __init__(self, integrator, parameter_manager):
        """
        Initialize single spectrum processor

        Args:
            integrator: Core integrator instance for this spectrum
            parameter_manager: Parameter manager for configuration
        """
        self.integrator = integrator
        self.parameter_manager = parameter_manager
        self.processing_active = False
        self.progress_callback = None

        # Processing statistics
        self.stats = {
            'total_processed': 0,
            'successful_fits': 0,
            'failed_fits': 0,
            'average_quality': 0.0
        }

    def process_peak_list(self, peak_list: pd.DataFrame,
                         processing_options: Optional[Dict] = None,
                         progress_callback: Optional[callable] = None,
                         locked_clusters_by_assignment: Optional[List[List[str]]] = None,
                         pre_learned_statistics: Optional[Dict] = None,
                         stop=None) -> Tuple[List[Dict], Optional[Dict]]:
        """
        Process all peaks in a peak list using specified options

        Args:
            peak_list: DataFrame with columns Position_X, Position_Y, Assignment
            processing_options: Dict with processing preferences
            progress_callback: Function to call for progress updates (progress, task, log_msg, failed)
            locked_clusters_by_assignment: Optional pre-computed clusters from reference spectrum.
                If provided, these clusters are used instead of recalculating, ensuring
                consistent peak grouping across all spectra in a series.
                Format: List of lists of assignment strings, e.g. [['Peak1', 'Peak2'], ['Peak3']]
            pre_learned_statistics: Optional linewidth statistics learned from reference spectrum.
                If provided, PASS 1 learning is skipped and these values are used as initial
                guesses (not fixed constraints). Format:
                {'lw_f1_median': float, 'lw_f2_median': float, 'lw_f1_mad': float,
                 'lw_f2_mad': float, 'n_samples': int, 'alpha': float}
            stop: Optional threading.Event; once set, the fit ends before its next cluster
                and raises Cancelled.

        Returns:
            Tuple of (fitted_results: List[Dict], learned_statistics: Optional[Dict])
            learned_statistics is the spectrum statistics learned during PASS 1 (None if skipped)
        """
        if processing_options is None:
            processing_options = {
                'mode': 'sequential',  # 'sequential', 'parallel', 'global_optimization'
                'use_parallel': self.parameter_manager.current_params.get('use_parallel_processing', False),
                'use_global_optimization': self.parameter_manager.current_params.get('use_global_optimization', False)
            }

        self.progress_callback = progress_callback
        self.processing_active = True
        self._locked_clusters_by_assignment = locked_clusters_by_assignment  # Store for internal methods
        self._pre_learned_statistics = pre_learned_statistics  # Store for parallel fitting
        # Store original reference positions for cascade mode absolute drift limiting
        self._original_reference_positions = processing_options.get('original_reference_positions') if processing_options else None
        # Store reference linewidths for per-peak linewidth reuse
        self._reference_linewidths = processing_options.get('reference_linewidths', {}) if processing_options else {}
        self._force_reference_linewidths = processing_options.get('force_reference_linewidths', False) if processing_options else False
        # Skip series_params to force fresh adaptive optimization (Independent mode)
        self._skip_series_params = processing_options.get('skip_series_params', False) if processing_options else False
        self._stop = stop

        try:
            # A fresh fit - nothing inherited from a reference spectrum - starts from the
            # base PS2D config, not from the radii a previous run's adaptive search left
            if locked_clusters_by_assignment is None and pre_learned_statistics is None:
                restore_base_ps2d_config()

            # Update integrator parameters from parameter manager
            self._sync_parameters_to_integrator()

            # Override adaptive optimization if skip_adaptive is requested
            if processing_options and processing_options.get('skip_adaptive', False):
                if hasattr(self.integrator, 'gui_params') and self.integrator.gui_params:
                    self.integrator.gui_params['use_adaptive_optimization'] = False

            # Apply ML predictions if available (provides initial parameter suggestions)
            # Skip during data collection to avoid confusing warnings
            if not (processing_options and processing_options.get('skip_ml_prediction', False)):
                ml_predictions = self._apply_ml_predictions(peak_list, processing_options)
            else:
                ml_predictions = None
            if ml_predictions is not None:
                # Store ML predictions in processing options for use by fitting methods
                if processing_options is None:
                    processing_options = {}
                processing_options['ml_predictions'] = ml_predictions

            # Choose processing strategy
            # Returns tuple: (fitted_results, learned_statistics)
            if processing_options.get('use_global_optimization', False):
                results = self._process_with_global_optimization(peak_list)
                return (results, None)  # Global optimization doesn't learn statistics
            elif processing_options.get('use_parallel', False):
                return self._process_with_parallel_fitting(peak_list)
            else:
                # Sequential is the same pipeline, run in this process
                return self._process_with_parallel_fitting(peak_list, in_process=True)

        except Cancelled:
            raise
        except Exception as e:
            log_error(f"Error in process_peak_list: {e}")
            if self.progress_callback:
                self.progress_callback(0, f"Processing failed: {e}", None, True)
            return ([], None)
        finally:
            self.processing_active = False

    def _sync_parameters_to_integrator(self):
        """Sync parameters from parameter manager to integrator"""

        params = self.parameter_manager.get_integrator_parameters()

        # Apply detection parameters
        detection_params = params['detection_params']
        self.integrator.set_search_window(
            detection_params['search_window_x'],
            detection_params['search_window_y']
        )
        self.integrator.set_threshold_multiplier(detection_params['noise_threshold'])

        # Apply fitting parameters
        self.integrator.fitting_parameters.update(params['fitting_params'])

        # Apply GUI parameters
        self.integrator.gui_params = params['gui_params']

        log_info("Single spectrum processor parameters synchronized")

    def _apply_ml_predictions(self, peak_list: pd.DataFrame,
                              processing_options: Optional[Dict] = None) -> Optional[Dict]:
        """
        Apply ML predictions to suggest optimal parameters.

        This method extracts spectrum features and uses the ML model (if available)
        to suggest optimal linewidths and integration windows. The predictions are
        used as initial guesses - PASS1 learning can still refine them.

        Parameters
        ----------
        peak_list : pd.DataFrame
            Peak list with positions
        processing_options : dict, optional
            Processing options that may include ML config

        Returns
        -------
        dict or None
            ML predictions if available and confident, None otherwise
        """
        ml_manager = _get_ml_manager()
        if ml_manager is None:
            return None

        # ML parameter prediction was never wired up (the model has no trained
        # weights and the feature extractor reads a spectrum attribute that does not
        # exist), so it is opt-in: skipped unless a caller explicitly enables it.
        ml_config = processing_options.get('ml_config', {}) if processing_options else {}
        if not ml_config.get('enabled', False):
            return None
        if not ml_config.get('use_predictions', False):
            return None

        try:
            from lunaNMR.ml import FeatureExtractor, SpectrumFeatures

            # Extract features from spectrum
            extractor = FeatureExtractor()

            # Get spectrum data from integrator
            spectrum_data = self.integrator.spectrum_data
            ppm_f1 = self.integrator.ppm_f1
            ppm_f2 = self.integrator.ppm_f2

            if spectrum_data is None or ppm_f1 is None or ppm_f2 is None:
                return None

            features = extractor.extract_from_spectrum(
                spectrum_data=spectrum_data,
                ppm_f1=ppm_f1,
                ppm_f2=ppm_f2,
                peak_count=len(peak_list),
            )

            # Get prediction
            prediction = ml_manager.predict(features)

            if prediction is None:
                return None

            # Check confidence threshold
            min_confidence = ml_config.get('ml_confidence_threshold', 0.5)
            if prediction.confidence < min_confidence:
                log_info(f"ML prediction confidence {prediction.confidence:.0%} below threshold")
                return None

            # Log prediction if configured
            if ml_config.get('show_predictions_in_log', True):
                log_info(
                    f"ML suggests: LW F1={prediction.lw_f1_median:.3f} ppm, "
                    f"LW F2={prediction.lw_f2_median:.3f} ppm "
                    f"(confidence: {prediction.confidence:.0%}, source: {prediction.source})"
                )

            return {
                'lw_f1_median': prediction.lw_f1_median,
                'lw_f2_median': prediction.lw_f2_median,
                'rad_f1': prediction.rad_f1,
                'rad_f2': prediction.rad_f2,
                'achievable_r2': prediction.achievable_r2,
                'confidence': prediction.confidence,
                'source': prediction.source,
            }

        except Exception as e:
            log_warning(f"ML prediction failed: {e}")
            return None

    def _process_with_parallel_fitting(self, peak_list: pd.DataFrame,
                                       in_process: bool = False) -> Tuple[List[Dict], Optional[Dict]]:
        """Fit a peak list with the cluster-based Voigt pipeline: PASS 1 on the isolated
        peaks, the adaptive radius search, PASS 2 on the clusters.

        in_process runs that pipeline in this process instead of over a pool - this is
        sequential mode, and the two modes therefore cannot disagree. A list too short to
        repay the cost of distributing it is fitted in process.

        Returns:
            Tuple of (fitted_results, learned_statistics)
        """

        if not in_process and len(peak_list) < MIN_PEAKS_FOR_PARALLEL:
            log_progress(f"Fitting {len(peak_list)} peak(s) sequentially: too few to distribute")
            in_process = True

        log_progress(f"Starting {'sequential' if in_process else 'parallel'} processing of "
                     f"{len(peak_list)} peaks")

        # Get locked clusters if available (for series integration)
        locked_clusters = getattr(self, '_locked_clusters_by_assignment', None)
        if locked_clusters is not None:
            log_info(f"Parallel mode: Using {len(locked_clusters)} locked clusters from reference")

        # Get pre-learned statistics if available (from reference spectrum)
        pre_learned_statistics = getattr(self, '_pre_learned_statistics', None)
        if pre_learned_statistics is not None:
            log_info("Parallel mode: Using pre-learned statistics (skipping PASS 1)")

        # Get reference linewidths if available (for per-peak linewidth reuse)
        reference_linewidths = getattr(self, '_reference_linewidths', {})
        force_reference_linewidths = getattr(self, '_force_reference_linewidths', False)
        if force_reference_linewidths and reference_linewidths:
            log_info(f"Parallel mode: Using {len(reference_linewidths)} reference linewidths")

        fallback_reason = "the cluster-based parallel fitter is unavailable"
        # Check if enhanced parallel fitting is available
        if (hasattr(self.integrator, 'enhanced_fitter') and
            hasattr(self.integrator.enhanced_fitter, 'enhanced_peak_fitting_parallel')):

            #print("✨ Using enhanced parallel Voigt fitting")

            # Define progress callback for parallel processing
            def parallel_progress_callback(progress, status, current_item):
                if self.progress_callback and self.processing_active:
                    self.progress_callback(progress, status, current_item)

            try:
                # Use new complete parallel implementation
                # Returns tuple: (fitted_results, learned_statistics)
                # skip_series_params forces fresh adaptive optimization (Independent mode)
                skip_series_params = getattr(self, '_skip_series_params', False)
                result = self.integrator.enhanced_fitter.enhanced_peak_fitting_parallel(
                    peak_list,
                    progress_callback=parallel_progress_callback,
                    locked_clusters_by_assignment=locked_clusters,
                    pre_learned_statistics=pre_learned_statistics,
                    reference_linewidths=reference_linewidths if force_reference_linewidths else None,
                    skip_series_params=skip_series_params,
                    reference_positions=getattr(self, '_original_reference_positions', None),
                    in_process=in_process,
                    stop=getattr(self, '_stop', None)
                )

                # Handle both old (list) and new (tuple) return formats
                if isinstance(result, tuple):
                    fitted_results, learned_statistics = result
                else:
                    fitted_results = result
                    learned_statistics = None
                
                # Ensure results is a list
                if not isinstance(fitted_results, list):
                    fitted_results = [fitted_results] if fitted_results else []
                
                # Add metadata to results
                for i, result in enumerate(fitted_results):
                    if result:
                        # A fallback's own mode is the one thing that says it happened
                        if 'fitting_fallback' not in result:
                            result['processing_mode'] = 'sequential' if in_process else 'enhanced_parallel'
                        if 'peak_number' not in result:
                            result['peak_number'] = i + 1

                log_info(f"Enhanced parallel processing completed: {len(fitted_results)} results")

                # CRITICAL: Copy clusters from enhanced_fitter.series_params for series propagation
                # This ensures the same cluster storage works for both parallel and sequential modes
                if (hasattr(self.integrator, 'enhanced_fitter') and
                    hasattr(self.integrator.enhanced_fitter, 'series_params') and
                    self.integrator.enhanced_fitter.series_params):
                    clusters_from_parallel = self.integrator.enhanced_fitter.series_params.get(
                        'locked_clusters_by_assignment'
                    )
                    if clusters_from_parallel:
                        self._computed_clusters_by_assignment = clusters_from_parallel
                        log_info(f"Stored {len(clusters_from_parallel)} clusters from parallel fitting for series propagation")

                return (fitted_results, learned_statistics)

            except Cancelled:
                raise
            except Exception as e:
                if in_process:
                    raise
                log_warning(f"Enhanced parallel processing failed: {e}")
                fallback_reason = f"{type(e).__name__}: {e}"
                # Fall through to existing parallel implementation

        log_warning(f"Falling back to the original parallel fitter, with no locked clusters "
                    f"or learned statistics, because {fallback_reason}")

        # Fallback to existing ParallelPeakFitter (unchanged)
        try:
            from lunaNMR.processors.parallel_fitting import ParallelPeakFitter

            log_info(f"Using original parallel fitting for {len(peak_list)} peaks")

            # Create parallel fitter
            parallel_fitter = ParallelPeakFitter(self.integrator)

            # Define progress callback for original parallel processing
            def original_parallel_progress_callback(completed, total, _current_assignment):
                if self.progress_callback and self.processing_active:
                    progress = (completed / total) * 100
                    self.progress_callback(
                        progress,
                        f"Original parallel fitting: {completed}/{total} completed",
                        "Processing peaks in parallel"
                    )

            # Run original parallel fitting
            fitted_results = parallel_fitter.fit_peaks_parallel(peak_list, original_parallel_progress_callback)

            # Add metadata to results
            for i, result in enumerate(fitted_results):
                if result:
                    result['processing_mode'] = 'original_parallel'
                    result['fitting_fallback'] = f"original parallel fitting ({fallback_reason})"
                    if 'peak_number' not in result:
                        result['peak_number'] = i + 1

            return (fitted_results, None)  # Original parallel doesn't learn statistics

        except Exception as e:
            log_error(f"Original parallel fitting also failed: {e}")
            log_info("Falling back to sequential processing")
            return self._process_with_parallel_fitting(peak_list, in_process=True)

    def _process_with_global_optimization(self, peak_list: pd.DataFrame) -> List[Dict]:
        """Process peaks using global optimization"""

        try:
            log_progress(f"Starting global optimization of {len(peak_list)} peaks")

            # Convert peak list to format expected by global optimization
            peak_tuples = []
            for i, (_, row) in enumerate(peak_list.iterrows()):
                peak_x = float(row['Position_X'])
                peak_y = float(row['Position_Y'])
                assignment = row.get('Assignment', f'Peak_{i+1}')
                peak_tuples.append((peak_x, peak_y, assignment))

            # Update progress
            if self.progress_callback:
                self.progress_callback(10, "Global optimization", "Analyzing peak patterns...")

            # Use global optimization from integrator
            optimization_report = self.integrator.optimize_peak_list_globally(
                peak_tuples,
                convergence_threshold=0.05,
                max_rounds=5
            )

            # Update progress
            if self.progress_callback:
                self.progress_callback(90, "Global optimization", "Finalizing results...")

            # Extract fitted results from optimization report
            fitted_results = []
            if 'fitted_results' in optimization_report:
                fitted_results = optimization_report['fitted_results']

                # Add metadata
                for i, result in enumerate(fitted_results):
                    if result:
                        result['processing_mode'] = 'global_optimization'
                        if 'peak_number' not in result:
                            result['peak_number'] = i + 1

            # Get summary statistics
            if 'optimization_summary' in optimization_report:
                summary = optimization_report['optimization_summary']
                success_rate = summary.get('final_success_rate', 0)
                optimization_rounds = summary.get('total_rounds', 0)

                log_info(f"Global optimization complete: {success_rate:.1f}% success rate after {optimization_rounds} rounds")

            return fitted_results

        except Exception as e:
            log_error(f"Global optimization failed: {e}")
            log_info("Falling back to sequential processing")
            return self._process_with_parallel_fitting(peak_list, in_process=True)[0]

    def get_processing_summary(self, fitted_results: List[Dict], total_peaks: int) -> Dict[str, Any]:
        """Generate comprehensive processing summary"""

        successful_count = len(fitted_results)
        success_rate = (successful_count / total_peaks * 100) if total_peaks > 0 else 0

        # Calculate quality statistics
        quality_scores = []
        for result in fitted_results:
            r_squared = result.get('avg_r_squared', 0)
            if r_squared > 0:
                quality_scores.append(r_squared)

        avg_quality = sum(quality_scores) / len(quality_scores) if quality_scores else 0

        # Categorize quality
        excellent_count = sum(1 for score in quality_scores if score >= 0.95)
        good_count = sum(1 for score in quality_scores if 0.85 <= score < 0.95)
        poor_count = sum(1 for score in quality_scores if score < 0.85)

        return {
            'total_peaks': total_peaks,
            'successful_peaks': successful_count,
            'failed_peaks': total_peaks - successful_count,
            'success_rate': success_rate,
            'average_r_squared': avg_quality,
            'quality_distribution': {
                'excellent': excellent_count,  # R² >= 0.95
                'good': good_count,           # 0.85 <= R² < 0.95
                'poor': poor_count            # R² < 0.85
            },
            'results': fitted_results,
            'processing_stats': self.stats.copy()
        }

    def cancel_processing(self):
        """Cancel current processing operation"""
        self.processing_active = False
        log_warning("Single spectrum processing cancelled")

    def extract_clusters_by_assignment(
        self,
        peak_list: pd.DataFrame
    ) -> List[List[str]]:
        """
        Extract clusters as assignment lists from a peak list.

        Call this after processing the reference spectrum to get clusters
        that can be locked for subsequent spectra.

        Args:
            peak_list: DataFrame with Position_X, Position_Y, Assignment columns

        Returns:
            List of clusters by assignment: [['Peak1', 'Peak2'], ['Peak3'], ...]
        """
        # Build all_peaks_context
        all_peaks_context = []
        for _, row in peak_list.iterrows():
            all_peaks_context.append({
                'assignment': str(row.get('Assignment', 'Unknown')),
                'x_ppm': float(row['Position_X']),
                'y_ppm': float(row['Position_Y']),
                'pos_x': float(row['Position_X']),
                'pos_y': float(row['Position_Y'])
            })

        # Get position-based clusters
        position_clusters = self.integrator.identify_overlap_clusters(all_peaks_context)

        # Build position → assignment mapping
        position_to_assignment = {}
        for peak in all_peaks_context:
            x = peak.get('x_ppm') or peak.get('pos_x')
            y = peak.get('y_ppm') or peak.get('pos_y')
            assignment = peak.get('assignment', 'Unknown')
            if x is not None and y is not None:
                position_to_assignment[(x, y)] = assignment

        # Convert position clusters to assignment clusters
        assignment_clusters = []
        for cluster_positions in position_clusters:
            cluster_assignments = []
            for pos in cluster_positions:
                assignment = position_to_assignment.get(pos)
                if assignment:
                    cluster_assignments.append(assignment)
            if cluster_assignments:
                assignment_clusters.append(cluster_assignments)

        log_info(f"Extracted {len(assignment_clusters)} clusters by assignment")
        for i, cluster in enumerate(assignment_clusters[:5]):  # Show first 5
            if len(cluster) > 1:
                log_info(f"  Cluster {i+1}: {cluster}")

        return assignment_clusters

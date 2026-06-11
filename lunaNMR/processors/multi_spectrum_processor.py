"""
Independent Multi-Spectrum Processor for NMR Series Analysis
Combines complete GUI decoupling with retrocompatibility features

Author: Guillaume Mas
Date: 2025
"""

import os
import json
import pandas as pd
import numpy as np
from datetime import datetime
from typing import List, Dict, Any, Optional

# Independent imports - no GUI dependencies
from lunaNMR.core.core_integrator import EnhancedVoigtIntegrator
from lunaNMR.utils.file_manager import NMRFileManager
from lunaNMR.processors.single_spectrum_processor import SingleSpectrumProcessor
from lunaNMR.utils.parameter_manager import NMRParameterManager
from lunaNMR.utils.output_manager import log_progress, log_info, log_warning, log_error

class MultiSpectrumProcessor:
    """
    Completely independent multi-spectrum processor.
    Has its own integrator instance and creates all necessary output files.
    """

    def __init__(self, voigt_params):
        """Initialize with independent integrator instance"""
        # Create OWN integrator instance - completely independent from GUI
        self.integrator = EnhancedVoigtIntegrator()
        self.file_manager = NMRFileManager()
        self.voigt_params = voigt_params
        self.reference_peaks = None
        self.processing_active = False
        self.progress_callback = None

        # Configure the internal integrator with provided parameters
        # OPTION C FIX: Extract gui_params correctly from nested structure
        if isinstance(self.voigt_params, dict) and 'gui_params' in self.voigt_params:
            # Nested structure from get_effective_parameters()
            self.integrator.gui_params = self.voigt_params['gui_params'].copy()
            self.integrator.fitting_parameters.update(
                self.voigt_params.get('fitting_params', {})
            )
        else:
            # Flat structure (backward compatibility)
            self.integrator.gui_params = self.voigt_params.copy()
            self.integrator.fitting_parameters.update(self.voigt_params)

        # Initialize output folder for retrocompatibility
        self.output_folder = None

        # PS2D linewidth reuse configuration
        self.use_ps2d_linewidth_reuse = self.voigt_params.get('use_ps2d_linewidth_reuse', False)
        if not self.use_ps2d_linewidth_reuse and 'gui_params' in self.voigt_params:
            # Try nested structure
            self.use_ps2d_linewidth_reuse = self.voigt_params['gui_params'].get('use_ps2d_linewidth_reuse', False)

        # Rerun adaptive optimization per spectrum (instead of reusing spectrum 1 stats)
        self.rerun_adaptive_per_spectrum = self.voigt_params.get('rerun_adaptive_per_spectrum', False)
        if not self.rerun_adaptive_per_spectrum and 'processing_options' in self.voigt_params:
            self.rerun_adaptive_per_spectrum = self.voigt_params['processing_options'].get('rerun_adaptive_per_spectrum', False)

        # Use original reference for detection (independent mode - no N-1 position propagation)
        self.use_original_reference_for_detection = self.voigt_params.get('use_original_reference_for_detection', False)
        if not self.use_original_reference_for_detection and 'processing_options' in self.voigt_params:
            self.use_original_reference_for_detection = self.voigt_params['processing_options'].get('use_original_reference_for_detection', False)

        # Lock cluster assignments from reference spectrum (for Independent mode)
        self.lock_cluster_assignments = self.voigt_params.get('lock_cluster_assignments', False)
        if not self.lock_cluster_assignments and 'processing_options' in self.voigt_params:
            self.lock_cluster_assignments = self.voigt_params['processing_options'].get('lock_cluster_assignments', False)

        self.reference_linewidths = {}  # Stores {assignment: {x_sigma, x_gamma, y_sigma, y_gamma}}

        # Cascade mode state (for sequential peak position propagation)
        self.peak_source_mode = None
        self.previous_fitted_results = None  # Stores fitted results from spectrum N for cascade to N+1

        # Cluster locking: Store reference spectrum clusters by assignment
        # This ensures consistent peak grouping across all spectra in series
        # CRITICAL: Initialize from voigt_params if clusters were passed from reference fit
        self.reference_clusters_by_assignment = self.voigt_params.get('reference_clusters_by_assignment', None)
        if self.reference_clusters_by_assignment:
            log_info(f"Using {len(self.reference_clusters_by_assignment)} clusters from reference fit")

        # Learned statistics from reference spectrum (PASS 1 learning)
        # Used as initial guesses for spectrum 2+ (not fixed, just better starting point)
        # Format: {'lw_f1_median': float, 'lw_f2_median': float, 'lw_f1_mad': float,
        #          'lw_f2_mad': float, 'n_samples': int, 'alpha': float}
        self.reference_spectrum_statistics = None

        # ML Training Data Collection (initialized in process_nmr_series with output_folder)
        self.training_collector = None
        self.collect_training_data = self.voigt_params.get('collect_training_data', False)
        if 'gui_params' in self.voigt_params:
            self.collect_training_data = self.voigt_params['gui_params'].get('collect_training_data', False)

    def process_nmr_series(self, nmr_files: List[str], reference_peaks: pd.DataFrame,
                          output_folder: str, peak_source_mode: str = 'reference',
                          progress_callback: Optional[callable] = None,
                          extract_delays: bool = False,
                          series_mode: str = "time",
                          pre_detected_peaks: Optional[List[Dict]] = None,
                          sn_from_gui_locked: Optional[Dict] = None,
                          peak_list_contour_min: Optional[float] = None,
                          lock_detection_threshold: bool = True) -> Dict[str, Any]:
        """
        Main entry point for processing NMR series.
        Returns comprehensive results in both new and legacy-compatible formats.

        Args:
            nmr_files: List of NMR spectrum file paths
            reference_peaks: DataFrame with peak positions
            output_folder: Folder for output files
            peak_source_mode: 'reference', 'detected', or 'cascade'
            progress_callback: Optional callback for progress updates
            extract_delays: If True, extract per-spectrum values from filenames and use as headers
            series_mode: 'time' (read _50ms/_1s delays) or 'titration' (read _1o0/_0o5 points)
            pre_detected_peaks: Optional list of already-detected peaks from GUI (skips detection for spectrum 1)
            sn_from_gui_locked: Optional dict with {'absolute_threshold': float, 'contour_min': float}
                               If provided, all spectra use same absolute threshold (no rescaling)
            peak_list_contour_min: Optional float for Peak List mode detection threshold
                               Each spectrum computes: threshold = max_intensity * contour_min * 1.05
            lock_detection_threshold: If True, use spectrum 1 threshold for all spectra (default True)
        """
        self.reference_peaks = reference_peaks.copy()
        self.pre_detected_peaks = pre_detected_peaks  # Store for cascade mode spectrum 1
        self.output_folder = output_folder
        self.progress_callback = progress_callback
        self.processing_active = True
        self.peak_source_mode = peak_source_mode  # Store for cascade mode
        self.previous_fitted_results = None  # Initialize cascade tracking
        self.extract_delays = extract_delays  # Store for output file generation
        self.series_mode = series_mode  # 'time' or 'titration' for value extraction
        self.nmr_files = nmr_files  # Store for ML training data collection
        self.sn_from_gui_locked = sn_from_gui_locked  # Store locked threshold for series
        self.peak_list_contour_min = peak_list_contour_min  # Store for Peak List mode detection
        self.lock_detection_threshold = lock_detection_threshold  # Lock threshold from spectrum 1
        self.locked_contour_threshold = None  # Will store spectrum 1 threshold if locking enabled

        # Initialize ML training collector with series output folder
        if self.collect_training_data and output_folder:
            self._init_training_collector(output_folder)

        # Build original reference position dictionary for cascade mode drift limiting
        # Maps assignment -> (x_ppm, y_ppm) for absolute bounds enforcement
        self.original_reference_positions = {}
        for _, row in reference_peaks.iterrows():
            assignment = row.get('Assignment', '')
            if assignment:
                self.original_reference_positions[str(assignment)] = (
                    row['Position_X'],
                    row['Position_Y']
                )

        # Build delay/column mapping if extract_delays is enabled
        # Maps spectrum_name -> unique column name (e.g., "50", "50_2" for duplicates)
        self.delay_mapping = {}
        if extract_delays:
            from lunaNMR.utils.delay_extractor import DelayExtractor
            extractor = DelayExtractor(mode=series_mode)

            # Get filenames
            filenames = [os.path.basename(f) for f in nmr_files]

            # Build column mapping (handles duplicates automatically)
            self.delay_mapping = extractor.build_column_mapping(filenames)



        # Track processing duration
        start_time = datetime.now()

        try:
            # Initialize comprehensive batch results structure
            batch_results = self._initialize_comprehensive_batch_results(
                nmr_files, output_folder, peak_source_mode
            )

            total_spectra = len(nmr_files)
            successful_spectra = 0
            total_peaks_processed = 0
            total_successful_fits = 0

            # Process each spectrum independently
            for i, nmr_file in enumerate(nmr_files):
                if not self.processing_active:
                    break

                spectrum_name = os.path.basename(nmr_file)
                spectrum_key = os.path.splitext(spectrum_name)[0]

                # Update progress
                progress = (i / total_spectra) * 100
                if self.progress_callback:
                    self.progress_callback(
                        progress,
                        f"Processing spectrum {i+1}/{total_spectra}",
                        f"Loading {spectrum_name}"
                    )

                try:
                    # Process single spectrum using independent integrator
                    spectrum_result = self._process_single_spectrum_independent(
                        nmr_file, spectrum_name, i+1, total_spectra
                    )

                    if spectrum_result['success']:
                        successful_spectra += 1
                        total_successful_fits += spectrum_result['successful_fits']

                    total_peaks_processed += spectrum_result['total_peaks']
                    batch_results['results'][spectrum_name] = spectrum_result

                except Exception as e:
                    log_error(f"Failed to process {spectrum_name}: {e}")
                    batch_results['results'][spectrum_name] = {
                        'success': False,
                        'error': str(e),
                        'fitted_results': [],
                        'total_peaks': len(self.reference_peaks),
                        'successful_fits': 0,
                        'spectrum_file': spectrum_name
                    }

            # Calculate processing duration
            end_time = datetime.now()
            duration = end_time - start_time
            duration_str = str(duration).split('.')[0]  # Remove microseconds

            # Finalize results with comprehensive statistics
            batch_results['summary'] = {
                'total_spectra': total_spectra,
                'successful': successful_spectra,
                'failed': total_spectra - successful_spectra,
                'success_rate': (successful_spectra / total_spectra * 100) if total_spectra > 0 else 0,
                'total_peaks_processed': total_peaks_processed,
                'total_successful_fits': total_successful_fits,
                'overall_detection_rate': (total_successful_fits / total_peaks_processed * 100) if total_peaks_processed > 0 else 0,
                'duration': duration_str,
                'processing_mode': 'series_integration'
            }

            # Create all output files for downstream compatibility
            self._create_comprehensive_output_files(batch_results)

            # Convert dict to BatchResults object for GUI compatibility
            batch_results_obj = self._convert_to_batch_results(batch_results, start_time, end_time)

            return batch_results_obj

        except Exception as e:
            log_error(f"Multi-spectrum processing failed: {e}")
            # Return BatchResults object even on error for consistent interface
            error_result = BatchResults()
            error_result.errors.append({'error': str(e)})
            return error_result
        finally:
            self.processing_active = False

    def _process_single_spectrum_independent(self, nmr_file: str, spectrum_name: str,
                                           spectrum_number: int, total_spectra: int) -> Dict[str, Any]:
        """Process a single spectrum using the independent integrator"""

        # REFINED STATE MANAGEMENT: Preserve essential integrator functionality
        # Save critical components before clearing
        essential_components = {}
        preserve_attrs = ['enhanced_fitter', 'fitting_parameters', 'gui_params', 'threshold_multiplier', 'processing_mode']

        for attr in preserve_attrs:
            if hasattr(self.integrator, attr):
                essential_components[attr] = getattr(self.integrator, attr)

        # Clear only problematic state attributes
        problematic_attrs = ['peak_list_path', 'nmr_file_path']
        for attr in problematic_attrs:
            if hasattr(self.integrator, attr):
                delattr(self.integrator, attr)

        # Clear NMR data only (let load_nmr_file handle proper loading)
        self.integrator.nmr_data = None
        self.integrator.nmr_dict = None
        self.integrator.ppm_x_axis = None
        self.integrator.ppm_y_axis = None

        # Always reset peak_list to None to ensure clean loading
        self.integrator.peak_list = None

        # Restore essential components
        for attr, value in essential_components.items():
            setattr(self.integrator, attr, value)

        # First attempt: Use integrator's load method
        # OPTIMIZATION: Skip nucleus auto-detection for spectra 2+ (use reference config)
        skip_detection = (spectrum_number > 1)
        load_success = False
        try:
            load_success = self.integrator.load_nmr_file(nmr_file, skip_nucleus_detection=skip_detection)
        except Exception:
            pass  # Will try fallback

        # Fallback: Direct loading only if regular loading fails
        if not load_success or self.integrator.nmr_data is None:
            try:
                import nmrglue as ng
                import os

                # Detect format and load accordingly
                basename = os.path.basename(nmr_file)
                bruker_pdata_names = ('2rr', '2ri', '2ir', '2ii', '1r', '1i',
                                      '3rrr', '3rri', '3rir', '3rii', '3irr', '3iri', '3iir', '3iii')

                if basename in bruker_pdata_names:
                    # Bruker processed data - use parent directory
                    pdata_dir = os.path.dirname(nmr_file)
                    self.integrator.nmr_dict, self.integrator.nmr_data = ng.bruker.read_pdata(dir=pdata_dir)
                    self.integrator._nmr_format = 'bruker'
                else:
                    # NMRPipe format (default)
                    self.integrator.nmr_dict, self.integrator.nmr_data = ng.pipe.read(nmr_file)
                    self.integrator._nmr_format = 'pipe'

                # Only calculate axes if not already done
                if not hasattr(self.integrator, 'ppm_x_axis') or self.integrator.ppm_x_axis is None:
                    self.integrator._calculate_ppm_axes()

                # Auto-detect nucleus type (fallback path) - ONLY for spectrum 1
                if not skip_detection and hasattr(self.integrator, '_detect_nucleus_type'):
                    detected_nucleus = self.integrator._detect_nucleus_type()
                    if detected_nucleus:
                        self.integrator.auto_detected_nucleus = detected_nucleus

                # Only estimate noise if not already done
                if not hasattr(self.integrator, 'noise_level'):
                    self.integrator._estimate_noise_level()

                self.integrator.nmr_file_path = nmr_file

            except Exception as e:
                log_error(f"Direct loading failed: {e}")
                raise IOError(f"All loading methods failed for: {nmr_file}")
##


        # VERIFICATION: Ensure NMR data is properly loaded
        if (self.integrator.nmr_data is None or
            self.integrator.ppm_x_axis is None or
            self.integrator.ppm_y_axis is None):
            log_error(f"NMR data verification failed for {spectrum_name}")
            raise IOError(f"NMR data not properly loaded after direct loading: {nmr_file}")

        # APPLY LOCKED THRESHOLD FROM GUI (if in from-GUI mode)
        if self.sn_from_gui_locked:
            self.integrator.sn_from_gui_mode = True
            self.integrator.sn_absolute_threshold = self.sn_from_gui_locked['absolute_threshold']
            self.integrator.sn_contour_min_value = self.sn_from_gui_locked['contour_min']

            if spectrum_number == 1:
                log_info(f"Spectrum {spectrum_number} (reference): Using from-GUI threshold "
                         f"{self.integrator.sn_absolute_threshold:.2e}")
            else:
                log_info(f"Spectrum {spectrum_number}: Using locked threshold from reference "
                         f"(absolute={self.integrator.sn_absolute_threshold:.2e})")

        # APPLY CONTOUR THRESHOLD FOR PEAK LIST MODE DETECTION
        if self.peak_list_contour_min is not None:
            if spectrum_number == 1:
                # Spectrum 1: Compute and optionally store threshold
                max_intensity = np.max(self.integrator.nmr_data)
                safety_factor = 1.05  # 5% margin to ensure visible peaks are detected
                self.integrator.peak_list_contour_threshold = max_intensity * self.peak_list_contour_min * safety_factor
                log_info(f"Spectrum {spectrum_number}: Peak List contour threshold={self.integrator.peak_list_contour_threshold:.2e} "
                         f"(contour_min={self.peak_list_contour_min:.5f}, max={max_intensity:.2e})")
                # Store for locking if enabled
                if self.lock_detection_threshold:
                    self.locked_contour_threshold = self.integrator.peak_list_contour_threshold
                    log_info(f"Detection threshold locked at {self.locked_contour_threshold:.2e}")
            elif self.lock_detection_threshold and self.locked_contour_threshold is not None:
                # Spectrum 2+: Use locked threshold from spectrum 1
                self.integrator.peak_list_contour_threshold = self.locked_contour_threshold
                log_info(f"Spectrum {spectrum_number}: Using locked threshold={self.locked_contour_threshold:.2e}")
            else:
                # Spectrum 2+ without locking: Compute per-spectrum threshold
                max_intensity = np.max(self.integrator.nmr_data)
                safety_factor = 1.05
                self.integrator.peak_list_contour_threshold = max_intensity * self.peak_list_contour_min * safety_factor
                log_info(f"Spectrum {spectrum_number}: Peak List contour threshold={self.integrator.peak_list_contour_threshold:.2e} "
                         f"(contour_min={self.peak_list_contour_min:.5f}, max={max_intensity:.2e})")

        # DETECTION-BASED MATCHING: Choose reference source based on mode
        if self.peak_source_mode == 'reference':
            # REFERENCE MODE: Always use original reference peaks
            self.integrator.peak_list = self.reference_peaks.copy()

        elif self.peak_source_mode == 'independent':
            # INDEPENDENT MODE: Always use original reference peaks for detection
            # Full detection + fitting pipeline runs for each spectrum
            log_info(f"Spectrum {spectrum_number}: Independent mode - using original reference for detection")
            self.integrator.peak_list = self.reference_peaks.copy()

        elif self.peak_source_mode in ['cascade', 'detected']:
            # CASCADE/DETECTED MODE: Use n-1 fitted results for spectrum 2+
            # UNLESS use_original_reference_for_detection is enabled
            if spectrum_number > 1 and self.previous_fitted_results is not None:
                if self.use_original_reference_for_detection:
                    log_info(f"Spectrum {spectrum_number}: Using original reference (override enabled)")
                    self.integrator.peak_list = self.reference_peaks.copy()
                else:
                    # Convert previous fitted results to reference DataFrame format
                    detection_reference = self._convert_fitted_results_to_reference_df(self.previous_fitted_results)
                    self.integrator.peak_list = detection_reference.copy()
            else:
                # Spectrum 1: Use provided reference peaks
                self.integrator.peak_list = self.reference_peaks.copy()
        else:
            # Fallback for unknown mode
            log_warning(f"Unknown peak_source_mode '{self.peak_source_mode}', using reference peaks")
            self.integrator.peak_list = self.reference_peaks.copy()

        # OPTION A: SKIP DETECTION FOR SPECTRUM 2+ IN CASCADE MODE
        # In cascade mode, spectrum 1 runs full detection to find peak positions.
        # Spectrum 2+ skips detection entirely and uses spectrum N-1's fitted positions directly.
        # This prevents detection from finding different positions in different spectra,
        # which was causing peaks to jump dramatically (e.g., peak 168.0 jumping 0.48 ppm in 15N).
        # EXCEPTION: If use_original_reference_for_detection is enabled, run detection anyway.

        if self.peak_source_mode == 'cascade' and spectrum_number > 1 and not self.use_original_reference_for_detection:
            # CASCADE MODE SPECTRUM 2+: Skip detection, use reference positions directly
            # peak_list already contains positions from previous spectrum's fitted results

            # INTENSITY SAMPLING: Read actual intensities from current spectrum
            self._sample_intensities_for_peak_list()

        elif self.peak_source_mode == 'cascade' and spectrum_number == 1 and self.pre_detected_peaks:
            # CASCADE MODE SPECTRUM 1 WITH PRE-DETECTED PEAKS:
            # User already ran detection in GUI - use those peaks directly, skip detection
            log_info(f"Using {len(self.pre_detected_peaks)} pre-detected peaks from GUI (skipping detection)")

            # Convert pre-detected peaks to DataFrame format
            detected_peak_data = []
            for peak in self.pre_detected_peaks:
                detected_peak_data.append({
                    'Assignment': peak.get('assignment', peak.get('Assignment', 'Unknown')),
                    'Position_X': peak.get('ppm_x', peak.get('peak_x', peak.get('Position_X', 0))),
                    'Position_Y': peak.get('ppm_y', peak.get('peak_y', peak.get('Position_Y', 0))),
                    'Height': peak.get('intensity', peak.get('Height', 0)),
                    'Intensity': peak.get('intensity', peak.get('Intensity', 0))
                })

            detected_df = pd.DataFrame(detected_peak_data)
            self.integrator.peak_list = detected_df

            # Sample intensities from current spectrum at these positions
            self._sample_intensities_for_peak_list()

        else:
            # SPECTRUM 1 WITHOUT PRE-DETECTED PEAKS OR NON-CASCADE MODE: Run full detection

            # Ensure noise level is estimated
            if not hasattr(self.integrator, 'noise_level') or self.integrator.noise_level is None:
                self.integrator._estimate_noise_level()

            # Only add dummy peaks for FIRST spectrum in series
            # Subsequent spectra should use the same peak list (including dummies from spectrum 1)
            original_auto_add_dummy = getattr(self.integrator, 'auto_add_dummy_peaks', False)
            if spectrum_number > 1:
                self.integrator.auto_add_dummy_peaks = False

            # Call the same detection method used by GUI - this finds actual peak positions
            detected_peaks = self.integrator._detect_peaks_reference_based()

            # Restore original setting
            self.integrator.auto_add_dummy_peaks = original_auto_add_dummy

            if detected_peaks and len(detected_peaks) > 0:
                log_info(f"Detection complete: {len(detected_peaks)} peaks matched and recentered")

                # Convert detected peaks (fitted_peaks) to DataFrame format for fitting
                # This is the same conversion that GUI does (line 3283-3286 in main_gui.py)
                detected_peak_data = []
                for peak in detected_peaks:
                    detected_peak_data.append({
                        'Assignment': peak.get('assignment', 'Unknown'),
                        'Position_X': peak.get('ppm_x', 0),
                        'Position_Y': peak.get('ppm_y', 0),
                        'Height': peak.get('intensity', 0),
                        'Intensity': peak.get('intensity', 0)
                    })

                detected_df = pd.DataFrame(detected_peak_data)

                # CRITICAL: Replace peak_list with detected positions (not reference positions)
                # This ensures fitting uses DETECTED/RECENTERED positions, just like GUI
                self.integrator.peak_list = detected_df

            else:
                log_warning("No peaks detected - using reference positions as fallback")
                # Keep the reference peak_list (already set above)

        # USE SHARED SingleSpectrumProcessor LOGIC (NO DUPLICATION)
        # Create a temporary parameter manager with current voigt_params
        param_manager = NMRParameterManager()

        # OPTION C FIX: Set current_params correctly based on structure
        if 'gui_params' in self.voigt_params:
            # Nested structure - flatten for parameter manager
            param_manager.current_params = {}
            param_manager.current_params.update(self.voigt_params.get('detection_params', {}))
            param_manager.current_params.update(self.voigt_params.get('fitting_params', {}))
            param_manager.current_params.update(self.voigt_params.get('gui_params', {}))
            param_manager.current_params.update(self.voigt_params.get('processing_options', {}))
        else:
            # Flat structure
            param_manager.current_params = self.voigt_params.copy()

        # Create SingleSpectrumProcessor instance
        single_processor = SingleSpectrumProcessor(self.integrator, param_manager)

        # Set processing active flag
        single_processor.processing_active = self.processing_active

        # Define series-specific progress callback wrapper
        def series_progress_callback(progress, task, log_msg=None, failed=False):
            if self.progress_callback:
                # Adjust progress to account for multi-spectrum context
                adjusted_progress = ((spectrum_number - 1) / total_spectra + progress / 100 / total_spectra) * 100
                self.progress_callback(
                    adjusted_progress,
                    f"Spectrum {spectrum_number}/{total_spectra}",
                    log_msg or task
                )

        # OPTION C FIX: Use PUBLIC API - same path as "Fit All Peaks"
        log_progress(f"Processing spectrum {spectrum_number}/{total_spectra}")

        # Extract parameters from appropriate level
        if 'gui_params' in self.voigt_params:
            use_parallel = self.voigt_params['gui_params'].get('use_parallel_processing', False)
            use_voigt = self.voigt_params.get('use_voigt_fitting', True)
        else:
            use_parallel = self.voigt_params.get('use_parallel_processing', False)
            use_voigt = self.voigt_params.get('use_voigt_fitting', True)

        # Get cascade drift limit setting (default ON)
        if 'processing_options' in self.voigt_params:
            enable_cascade_drift_limit = self.voigt_params['processing_options'].get('enable_cascade_drift_limit', True)
        else:
            enable_cascade_drift_limit = self.voigt_params.get('enable_cascade_drift_limit', True)

        # Only pass reference positions if cascade/independent mode AND drift limit enabled
        use_drift_limit = (self.peak_source_mode in ['cascade', 'independent'] and enable_cascade_drift_limit)

        processing_options = {
            'use_parallel': use_parallel,
            'use_global_optimization': False,
            'use_voigt_fitting': use_voigt,
            # Pass original reference positions for drift limiting
            # Only if cascade/independent mode AND drift limit enabled
            'original_reference_positions': self.original_reference_positions if use_drift_limit else None,
            'enable_cascade_drift_limit': enable_cascade_drift_limit,
            # Skip series_params to allow fresh adaptive optimization per spectrum
            # (used for Independent mode where each spectrum runs full pipeline)
            'skip_series_params': self.rerun_adaptive_per_spectrum
        }

        # PER-PEAK LINEWIDTH REUSE: Pass reference linewidths for spectrum 2+
        # Spectrum 1: linewidths are extracted and stored in self.reference_linewidths
        # Spectrum 2+: use stored linewidths as initial guesses (optimizer can adjust)
        if self.use_ps2d_linewidth_reuse and spectrum_number > 1 and self.reference_linewidths:
            processing_options['reference_linewidths'] = self.reference_linewidths
            processing_options['force_reference_linewidths'] = True  # Use as initial guess
            log_info(f"Spectrum {spectrum_number}: Using {len(self.reference_linewidths)} reference linewidths as initial guess")

        if self.peak_source_mode == 'cascade':
            pass  # Cascade mode uses drift limiting

        # CLUSTER LOCKING: Use reference spectrum clusters for consistent grouping
        # Spectrum 1: compute clusters fresh, then store them
        # Spectrum 2+: use locked clusters from spectrum 1 (except in independent mode without lock)
        locked_clusters = None
        if spectrum_number > 1 and self.reference_clusters_by_assignment is not None:
            if self.peak_source_mode == 'independent':
                if self.lock_cluster_assignments:
                    # INDEPENDENT MODE with lock: Use reference clusters for consistent grouping
                    locked_clusters = self.reference_clusters_by_assignment
                    log_info(f"Spectrum {spectrum_number}: Using locked clusters from reference (Independent mode)")
                else:
                    # INDEPENDENT MODE: Each spectrum computes its own clusters (true GUI behavior)
                    log_info(f"Spectrum {spectrum_number}: Independent mode - computing fresh clusters")
            else:
                # CASCADE/REFERENCE MODE: Use reference clusters for consistent grouping
                locked_clusters = self.reference_clusters_by_assignment
                log_info(f"Spectrum {spectrum_number}: Using {len(locked_clusters)} locked clusters from spectrum 1")


        # LEARNED STATISTICS: Use reference spectrum statistics as initial guesses
        # Spectrum 1: run full learning cycle (PASS 1 → PASS 1-bis → PASS 2)
        # Spectrum 2+: skip PASS 1 learning, use reference stats as initial guesses
        #              UNLESS rerun_adaptive_per_spectrum is enabled
        pre_learned_statistics = None
        if spectrum_number > 1 and self.reference_spectrum_statistics is not None:
            if self.rerun_adaptive_per_spectrum:
                log_info(f"Spectrum {spectrum_number}: Rerunning PASS 1 + adaptive (per-spectrum mode)")
            else:
                pre_learned_statistics = self.reference_spectrum_statistics

        # This calls _sync_parameters_to_integrator() internally, ensuring proper parameter flow
        # CRITICAL: Use integrator.peak_list which now contains DETECTED positions (not reference_peaks)
        fitted_results, learned_statistics = single_processor.process_peak_list(
            self.integrator.peak_list,
            processing_options,
            series_progress_callback,
            locked_clusters_by_assignment=locked_clusters,
            pre_learned_statistics=pre_learned_statistics
        )

        # Extract and store clusters from reference spectrum (spectrum 1)
        # CRITICAL: Use the ACTUAL clusters from fitting, not recomputed clusters
        if spectrum_number == 1 and self.reference_clusters_by_assignment is not None:
            # Clusters already set from reference fit passed via voigt_params
            log_info(f"Spectrum 1: Using {len(self.reference_clusters_by_assignment)} pre-set clusters from reference fit")
        elif spectrum_number == 1 and self.reference_clusters_by_assignment is None:
            # Try multiple sources for clusters, in order of preference:
            clusters_found = None

            # Source 1: SingleSpectrumProcessor's computed clusters (works for sequential AND parallel)
            if hasattr(single_processor, '_computed_clusters_by_assignment') and single_processor._computed_clusters_by_assignment is not None and len(single_processor._computed_clusters_by_assignment) > 0:
                clusters_found = single_processor._computed_clusters_by_assignment
                log_info(f"Using ACTUAL clusters from single_processor: {len(clusters_found)} clusters")

            # Source 2: ParallelVoigtProcessor's series_params (parallel mode only)
            if clusters_found is None:
                if (hasattr(self.integrator, 'enhanced_fitter') and
                    hasattr(self.integrator.enhanced_fitter, 'series_params') and
                    self.integrator.enhanced_fitter.series_params):
                    clusters_found = self.integrator.enhanced_fitter.series_params.get(
                        'locked_clusters_by_assignment'
                    )
                    if clusters_found:
                        log_info(f"Using clusters from enhanced_fitter.series_params: {len(clusters_found)} clusters")

            # Source 3: Fallback - recompute (should not happen anymore)
            if clusters_found is None:
                log_warning("Fallback: No clusters found from fitting, recomputing (may differ!)")
                from lunaNMR.core.ps2d_config import get_ps2d_config
                config = get_ps2d_config()
                if learned_statistics is not None:
                    overlap_x = learned_statistics.get('overlap_threshold_x')
                    overlap_y = learned_statistics.get('overlap_threshold_y')
                    if overlap_x is not None and overlap_y is not None:
                        config.overlap_threshold_x = overlap_x
                        config.overlap_threshold_y = overlap_y

                clusters_found = single_processor.extract_clusters_by_assignment(
                    self.integrator.peak_list
                )

            self.reference_clusters_by_assignment = clusters_found
            log_info(f"Spectrum 1: Stored {len(clusters_found) if clusters_found else 0} clusters for series propagation")

        # Store learned statistics from reference spectrum (spectrum 1)
        if spectrum_number == 1 and learned_statistics is not None:
            self.reference_spectrum_statistics = learned_statistics

        # Add dummy peaks to original_reference_positions (spectrum 1 only)
        # This enables drift limiting for auto-detected dummy peaks in subsequent spectra
        if spectrum_number == 1:
            dummy_count = 0
            for result in fitted_results:
                if result and result.get('success', False):
                    assignment = result.get('assignment', '')
                    if isinstance(assignment, str) and assignment.startswith('dummy_'):
                        # Use detected position as reference (not fitted position)
                        pos_x = result.get('peak_x', result.get('center_x', 0))
                        pos_y = result.get('peak_y', result.get('center_y', 0))
                        if pos_x != 0 or pos_y != 0:
                            self.original_reference_positions[assignment] = (pos_x, pos_y)
                            dummy_count += 1
            if dummy_count > 0:
                log_info(f"Added {dummy_count} dummy peaks to drift limit reference positions")

        # Post-process results: Add spectrum metadata and handle linewidth reuse
        successful_fits = 0
        for result in fitted_results:
            if result and result.get('success', False):
                # Add spectrum-specific metadata
                result['spectrum_file'] = spectrum_name
                result['processing_mode'] = 'series_integration'
                successful_fits += 1

                # Extract linewidths from reference spectrum (first spectrum only)
                if spectrum_number == 1 and self.use_ps2d_linewidth_reuse:
                    assignment = result.get('assignment', 'Unknown')
                    self._extract_and_store_linewidth(result, assignment)

        total_peaks = len(self.integrator.peak_list)
        success_rate = (successful_fits / total_peaks * 100) if total_peaks > 0 else 0
        log_progress(f"Spectrum {spectrum_number} complete: {successful_fits}/{total_peaks} successful ({success_rate:.1f}%)")

        # Store fitted results for next spectrum (used by all modes now)
        self.previous_fitted_results = fitted_results

        # Collect ML training data if enabled
        if self.collect_training_data and self.training_collector is not None:
            self._collect_spectrum_training_data(
                spectrum_name=spectrum_name,
                spectrum_number=spectrum_number,
                fitted_results=fitted_results,
                learned_statistics=learned_statistics,
                total_spectra=len(self.nmr_files),
            )

        return {
            'success': successful_fits > 0,
            'fitted_results': fitted_results,
            'total_peaks': total_peaks,
            'successful_fits': successful_fits,
            'success_rate': success_rate,
            'spectrum_file': spectrum_name,
            'integration_results': self._convert_to_integration_format(fitted_results)
        }

    def _convert_fitted_results_to_reference_df(self, fitted_results):
        """
        Convert fitted results to reference peak DataFrame format.

        This creates a reference peak list from fitted results that can be used
        as input for peak detection in the next spectrum.

        IMPORTANT: Uses DETECTED positions (peak_x/peak_y) instead of FITTED positions
        (center_x/center_y) to prevent position drift accumulation in cascade mode.
        Fitted positions can drift during 2D optimization, especially for overlapping
        clusters, causing errors to compound across spectra.

        Args:
            fitted_results: List of fitting result dictionaries

        Returns:
            DataFrame with columns: Assignment, Position_X, Position_Y, Height, Intensity
        """
        if fitted_results is None or len(fitted_results) == 0:
            return self.reference_peaks.copy()

        peak_data = []
        for result in fitted_results:
            if result and result.get('success', False):
                # CRITICAL FIX: Use DETECTED position (peak_x, peak_y) NOT fitted position
                # (center_x, center_y). Fitted positions can drift during 2D multi-peak
                # optimization, causing cascade errors to compound.
                # peak_x/peak_y = where detection found the peak in the current spectrum
                # center_x/center_y = where 2D fitter placed the peak (can be wrong)
                pos_x = result.get('peak_x', result.get('center_x', 0))
                pos_y = result.get('peak_y', result.get('center_y', 0))

                peak_data.append({
                    'Assignment': result.get('assignment', 'Unknown'),
                    'Position_X': pos_x,
                    'Position_Y': pos_y,
                    'Height': result.get('height', result.get('amplitude', 0)),
                    'Intensity': result.get('volume', result.get('intensity', 0))
                })
            else:
                # Failed fit - use original reference position or previous position for dummies
                assignment = result.get('assignment', 'Unknown') if result else 'Unknown'

                # For dummy peaks, use position from result (previous fitted position)
                # since they don't exist in original reference_peaks
                if isinstance(assignment, str) and assignment.startswith('dummy_'):
                    pos_x = result.get('peak_x', result.get('center_x', 0)) if result else 0
                    pos_y = result.get('peak_y', result.get('center_y', 0)) if result else 0
                    if pos_x != 0 or pos_y != 0:
                        peak_data.append({
                            'Assignment': assignment,
                            'Position_X': pos_x,
                            'Position_Y': pos_y,
                            'Height': result.get('height', 0) if result else 0,
                            'Intensity': result.get('intensity', 0) if result else 0
                        })
                    continue

                # Try to find in original reference peaks
                matching_ref = self.reference_peaks[
                    self.reference_peaks['Assignment'] == assignment
                ]

                if not matching_ref.empty:
                    peak_data.append({
                        'Assignment': assignment,
                        'Position_X': matching_ref['Position_X'].iloc[0],
                        'Position_Y': matching_ref['Position_Y'].iloc[0],
                        'Height': matching_ref.get('Height', pd.Series([0])).iloc[0] if 'Height' in matching_ref.columns else 0,
                        'Intensity': matching_ref.get('Intensity', pd.Series([0])).iloc[0] if 'Intensity' in matching_ref.columns else 0
                    })

        if not peak_data:
            return self.reference_peaks.copy()

        return pd.DataFrame(peak_data)

    def _create_tidy_results_file(self, batch_results: Dict[str, Any]):
        """
        Creates a single, tidy CSV file from all series results.
        This format is ideal for external analysis and plotting.

        When extract_delays is enabled, the spectrum_name column uses delay values (e.g., "50")
        instead of filenames. Duplicate delays get unique suffixes (e.g., "50_2" for the second
        spectrum at 50ms). This improves compatibility with DynamiXs.
        """
        if not batch_results.get('results'):
            return

        tidy_data = []
        for spectrum_name, result_data in batch_results['results'].items():
            if not result_data.get('success', False):
                continue

            # Determine column header: use delay column name if extract_delays enabled
            if getattr(self, 'extract_delays', False) and hasattr(self, 'delay_mapping'):
                column_name = self.delay_mapping.get(spectrum_name)
                if column_name is not None:
                    # Use unique column name (e.g., "50" or "50_2" for duplicates)
                    spectrum_header = column_name
                else:
                    # Fallback to filename without extension
                    spectrum_header = os.path.splitext(spectrum_name)[0]
            else:
                # Use filename without extension
                spectrum_header = os.path.splitext(spectrum_name)[0]

            # Use the standardized integration results
            integration_results = result_data.get('integration_results', [])
            for peak in integration_results:
                row = {
                    'spectrum_name': spectrum_header,
                    'assignment': peak.get('assignment', 'Unknown'),
                    'peak_number': peak.get('peak_number', 0),
                    'ppm_x': peak.get('ppm_x', 0.0),
                    'ppm_y': peak.get('ppm_y', 0.0),
                    'height': peak.get('height', 0.0),
                    'volume': peak.get('volume', 0.0),
                    'snr': peak.get('snr', 0.0),
                    'quality': peak.get('quality', 'Unknown'),
                    'r_squared': peak.get('r_squared', 0.0)
                }
                tidy_data.append(row)

        if not tidy_data:
            return

        tidy_df = pd.DataFrame(tidy_data)
        tidy_file = os.path.join(self.output_folder, "series_analysis_tidy.csv")
        # Format ppm_x and ppm_y with 3 decimals (position data), other numeric columns with 1 decimal
        tidy_df_formatted = tidy_df.copy()
        for col in tidy_df_formatted.columns:
            if not pd.api.types.is_numeric_dtype(tidy_df_formatted[col]):
                continue
            # Position columns get 3 decimal places
            if col in ['ppm_x', 'ppm_y']:
                tidy_df_formatted[col] = tidy_df_formatted[col].apply(lambda x: f'{x:.3f}' if pd.notna(x) else '')
            # All other numeric columns get 1 decimal place
            else:
                tidy_df_formatted[col] = tidy_df_formatted[col].apply(lambda x: f'{x:.1f}' if pd.notna(x) else '')
        tidy_df_formatted.to_csv(tidy_file, index=False)

        log_info(f"Created tidy results file: {tidy_file}")

    def _convert_to_batch_results(self, batch_results_dict: Dict[str, Any],
                                   start_time, end_time) -> 'BatchResults':
        """Convert dict batch_results to BatchResults object for GUI compatibility.

        The dialogs expect a BatchResults object with .results attribute and get_summary() method,
        but internal processing uses dicts for flexibility. This converts at the end.

        Args:
            batch_results_dict: Dictionary with 'results' and 'summary' keys
            start_time: Processing start time
            end_time: Processing end time

        Returns:
            BatchResults object with same data
        """
        batch_obj = BatchResults()

        # Set metadata
        batch_obj.metadata['start_time'] = start_time
        batch_obj.metadata['end_time'] = end_time
        batch_obj.metadata['processing_mode'] = 'series_integration'
        batch_obj.metadata['total_spectra'] = batch_results_dict.get('summary', {}).get('total_spectra', 0)

        # Copy data_folder, output_folder, csv_path from batch_results_dict metadata
        source_metadata = batch_results_dict.get('metadata', {})
        if 'data_folder' in source_metadata:
            batch_obj.metadata['data_folder'] = source_metadata['data_folder']
        if 'output_folder' in source_metadata:
            batch_obj.metadata['output_folder'] = source_metadata['output_folder']
        if 'csv_path' in source_metadata:
            batch_obj.metadata['csv_path'] = source_metadata['csv_path']

        # Copy results - mark status as 'success' or 'failed' for each spectrum
        for spectrum_name, result in batch_results_dict.get('results', {}).items():
            # Ensure result has a 'status' field for the dialogs
            if 'status' not in result:
                result['status'] = 'success' if result.get('success', False) else 'failed'
            batch_obj.add_result(spectrum_name, result)

        # Copy statistics if available
        if 'statistics' in batch_results_dict:
            batch_obj.statistics = batch_results_dict['statistics']

        return batch_obj

    def _create_comprehensive_output_files(self, batch_results: Dict[str, Any]):
        """
        Creates all output files, including the new tidy format.
        (This is an update to the existing method).
        """
        if not os.path.exists(self.output_folder):
            os.makedirs(self.output_folder)

        # --- Existing code ---
        tracking_df = self._create_peak_tracking_table(batch_results)
        if not tracking_df.empty:
            tracking_file = os.path.join(self.output_folder, "comprehensive_peak_tracking.csv")
            # Format with 3 decimals for Reference_X/Y, 1 decimal for other columns
            tracking_df_formatted = self._format_dataframe_for_csv(tracking_df)
            tracking_df_formatted.to_csv(tracking_file, index=False)
            self._create_intensity_matrix(tracking_df)
            self._create_detection_matrix(tracking_df)
            self._create_summary_statistics_file(batch_results)

        # --- ADD THIS CALL ---
        self._create_tidy_results_file(batch_results)

## gm added

    def _convert_to_integration_format(self, fitted_results: List[Dict]) -> List[Dict]:
        """
        Convert detailed Voigt fitting results to the new standardized peak data format.
        This format is used by all downstream components, including the spectrum viewer.
        FIXED: Uses reference coordinates when Voigt coordinates are missing/zero.
        """
        standardized_results = []
        for i, fit_result in enumerate(fitted_results):
            if not isinstance(fit_result, dict):
                continue

            # Get assignment for reference lookup
            assignment = fit_result.get('assignment', f'Peak_{i+1}')

            # Find reference coordinates from the original peak list
            ref_x, ref_y = 0.0, 0.0
            if self.reference_peaks is not None:
                try:
                    # CRITICAL: Convert assignment to string for consistent comparison
                    # (same fix as Peak Navigator - assignments can be float or string)
                    assignment_str = str(assignment)

                    # Also convert reference peaks assignments to strings for comparison
                    matching_ref = self.reference_peaks[
                        self.reference_peaks['Assignment'].astype(str) == assignment_str
                    ]
                    if not matching_ref.empty:
                        ref_x = float(matching_ref['Position_X'].iloc[0])
                        ref_y = float(matching_ref['Position_Y'].iloc[0])
                except Exception as e:
                    log_warning(f"Error getting reference coordinates for {assignment}: {e}")

            # Initialize with reference coordinates as defaults
            peak_data = {
                'assignment': assignment,
                'peak_number': fit_result.get('peak_number', i + 1),
                'ppm_x': ref_x,  # Use reference coordinates as default
                'ppm_y': ref_y,  # Use reference coordinates as default
                'height': 0.0,
                'volume': 0.0,
                'snr': 0.0,
                'quality': 'Failed',
                'r_squared': 0.0,
                'detected': True,   # Always mark as detected since it was processed
                'fitted': False,
                'voigt_fit_data': None,

                # Add legacy format fields for backward compatibility
                'Assignment': assignment,
                'Position_X': ref_x,
                'Position_Y': ref_y,
                'Height': 0.0,
                'Volume': 0.0,
                'SNR': 0.0,
                'Quality': 'Failed',
                'R_Squared': 0.0,
                'Peak_Number': fit_result.get('peak_number', i + 1)
            }

            #success = fit_result.get('success', False)
            # Check multiple success indicators - be more permissive
            success = (fit_result.get('success', False) or
                      fit_result.get('fitted', False) or
                      bool(fit_result.get('x_fit')) or
                      bool(fit_result.get('y_fit')) or
                      fit_result.get('snr', 0) > 0 or
                      fit_result.get('avg_r_squared', 0) > 0)
            if success:
                # Extract data from the complex fit_result dictionary
                x_fit = fit_result.get('x_fit', {})
                y_fit = fit_result.get('y_fit', {})

                # Get fitted coordinates with robust fallback
                fitted_x = None
                fitted_y = None

                # Try multiple ways to get fitted coordinates
                if x_fit and 'center' in x_fit and x_fit['center'] is not None:
                    fitted_x = float(x_fit['center'])
                elif 'ppm_x' in fit_result and fit_result['ppm_x'] is not None:
                    fitted_x = float(fit_result['ppm_x'])
                elif 'peak_position' in fit_result and fit_result['peak_position']:
                    try:
                        fitted_x = float(fit_result['peak_position'][0])
                    except Exception:
                        pass

                if y_fit and 'center' in y_fit and y_fit['center'] is not None:
                    fitted_y = float(y_fit['center'])
                elif 'ppm_y' in fit_result and fit_result['ppm_y'] is not None:
                    fitted_y = float(fit_result['ppm_y'])
                elif 'peak_position' in fit_result and fit_result['peak_position']:
                    try:
                        fitted_y = float(fit_result['peak_position'][1])
                    except Exception:
                        pass

                # Use fitted coordinates if available and valid, otherwise keep reference
                final_x = fitted_x if fitted_x is not None and fitted_x != 0.0 else ref_x
                final_y = fitted_y if fitted_y is not None and fitted_y != 0.0 else ref_y

                # Height extraction - try multiple sources
                # IMPORTANT: For PS2D fits, 'height' is the calculated peak height at center
                # For 1D fits, 'amplitude' is the peak height
                height = 0.0
                if 'height' in fit_result:
                    # PS2D fitter returns 'height' (calculated from intensity and linewidths)
                    height = float(fit_result['height'])
                elif 'amplitude' in fit_result:
                    # Direct amplitude field
                    height = float(fit_result['amplitude'])
                elif x_fit and 'amplitude' in x_fit:
                    # 1D cross-section fit
                    height = float(x_fit['amplitude'])
                elif y_fit and 'amplitude' in y_fit:
                    # 1D cross-section fit
                    height = float(y_fit['amplitude'])
                elif 'peak_intensity' in fit_result:
                    height = float(fit_result['peak_intensity'])

                # Extract detected intensity (from peak detection, not fitting)
                detected_intensity = 0.0
                if 'detected_intensity' in fit_result and fit_result['detected_intensity'] is not None:
                    detected_intensity = float(fit_result['detected_intensity'])
                # Fallback: use intensity from all_peaks_context if available
                elif self.reference_peaks is not None:
                    try:
                        matching_peak = self.reference_peaks[
                            self.reference_peaks['Assignment'] == assignment
                        ]
                        if not matching_peak.empty:
                            if 'Height' in matching_peak.columns:
                                detected_intensity = float(matching_peak['Height'].iloc[0])
                            elif 'Intensity' in matching_peak.columns:
                                detected_intensity = float(matching_peak['Intensity'].iloc[0])
                    except Exception:
                        pass

                # Extract linewidths (FWHM in ppm)
                # NOTE: sigma/gamma in x_fit/y_fit are ALREADY FWHM values (not true sigma/gamma)
                # So we just add them directly (matching Peak Parameters display)
                lw_x, lw_y = 0.0, 0.0
                if x_fit:
                    x_sigma = x_fit.get('sigma', 0)  # Actually Gaussian FWHM
                    x_gamma = x_fit.get('gamma', 0)  # Actually Lorentzian FWHM
                    lw_x = x_sigma + x_gamma  # Simple sum of FWHM components

                if y_fit:
                    y_sigma = y_fit.get('sigma', 0)  # Actually Gaussian FWHM
                    y_gamma = y_fit.get('gamma', 0)  # Actually Lorentzian FWHM
                    lw_y = y_sigma + y_gamma  # Simple sum of FWHM components

                # SNR extraction - try multiple sources
                snr = 0.0
                if 'snr' in fit_result:
                    snr = float(fit_result['snr'])
                elif 'peak_snr' in fit_result:
                    snr = float(fit_result['peak_snr'])
                elif height > 0 and 'noise_level' in fit_result:
                    snr = height / fit_result['noise_level']

                # R-squared extraction - try multiple sources
                r_squared = 0.0
                if 'avg_r_squared' in fit_result:
                    r_squared = float(fit_result['avg_r_squared'])
                elif 'r_squared' in fit_result:
                    r_squared = float(fit_result['r_squared'])
                elif x_fit and 'r_squared' in x_fit:
                    r_squared = float(x_fit['r_squared'])
                elif y_fit and 'r_squared' in y_fit:
                    r_squared = float(y_fit['r_squared'])

                # Quality determination with fallbacks
                quality = 'Failed'
                if r_squared > 0.90:
                    quality = 'Excellent'
                elif r_squared > 0.80:
                    quality = 'Good'
                elif r_squared > 0.60:
                    quality = 'Fair'
                elif snr > 5.0:
                    quality = 'Poor'
                elif 'fitting_quality' in fit_result:
                    quality = str(fit_result['fitting_quality'])

                # Volume extraction - CRITICAL: Use volume from fit_result if available
                # For PS2D fits, 'volume' is the fitted intensity parameter
                # For 1D fits, calculate from height and linewidths
                volume = 0.0
                if 'volume' in fit_result:
                    # PS2D fitter returns 'volume' (fitted intensity parameter)
                    volume = float(fit_result['volume'])
                elif 'intensity' in fit_result:
                    # Alternative name for volume
                    volume = float(fit_result['intensity'])
                elif x_fit and y_fit and height > 0:
                    # Calculate volume from Voigt parameters for 1D fits
                    x_sigma = x_fit.get('sigma', 0)
                    x_gamma = x_fit.get('gamma', 0)
                    y_sigma = y_fit.get('sigma', 0)
                    y_gamma = y_fit.get('gamma', 0)

                    if (x_sigma > 0 or x_gamma > 0) and (y_sigma > 0 or y_gamma > 0):
                        x_fwhm = 2 * np.sqrt(2 * np.log(2)) * x_sigma + 2 * x_gamma
                        y_fwhm = 2 * np.sqrt(2 * np.log(2)) * y_sigma + 2 * y_gamma
                        volume = abs(height) * x_fwhm * y_fwhm

                # Update peak data with fitted values
                peak_data.update({
                    'ppm_x': final_x,
                    'ppm_y': final_y,
                    'height': height,
                    'volume': volume,
                    'snr': snr,
                    'quality': quality,
                    'r_squared': r_squared,
                    'fitted': True,
                    'detected_intensity': detected_intensity,
                    'lw_x': lw_x,
                    'lw_y': lw_y,
                    'voigt_fit_data': {
                        'x_fit': x_fit,
                        'y_fit': y_fit
                    },

                    # Update legacy fields as well
                    'Position_X': final_x,
                    'Position_Y': final_y,
                    'Height': height,
                    'Volume': volume,  # FIXED: Now correctly uses PS2D volume
                    'SNR': snr,
                    'Quality': quality,
                    'R_Squared': r_squared
                })

                # CRITICAL: Preserve 2D visualization data for Voigt analysis tabs
                # These fields are needed by VoigtAnalysisPlotter for 2D/3D visualization
                if 'method' in fit_result:
                    peak_data['method'] = fit_result['method']
                if 'region_2d' in fit_result:
                    peak_data['region_2d'] = fit_result['region_2d']
                if 'fitted_2d_surface' in fit_result:
                    peak_data['fitted_2d_surface'] = fit_result['fitted_2d_surface']
                if 'individual_surfaces' in fit_result:
                    peak_data['individual_surfaces'] = fit_result['individual_surfaces']
                if 'baseline' in fit_result:
                    peak_data['baseline'] = fit_result['baseline']
                if 'all_peaks' in fit_result:
                    peak_data['all_peaks'] = fit_result['all_peaks']

            standardized_results.append(peak_data)

        return standardized_results

##

##

    def _initialize_comprehensive_batch_results(self, nmr_files: List[str],
                                              output_folder: str, peak_source_mode: str) -> Dict[str, Any]:
        """Initialize comprehensive batch results structure"""
        # Derive data folder from first NMR file path
        data_folder = os.path.dirname(nmr_files[0]) if nmr_files else ''

        return {
            'metadata': {
                'processing_started': datetime.now().isoformat(),
                'peak_source_mode': peak_source_mode,
                'total_spectra': len(nmr_files),
                'output_folder': output_folder,
                'csv_path': os.path.join(output_folder, "series_analysis_tidy.csv"),
                'data_folder': data_folder,  # Store folder containing NMR files
                'reference_peaks_count': len(self.reference_peaks)
            },
            'results': {},
            'summary': {}
        }

    def _create_comprehensive_output_files(self, batch_results: Dict[str, Any]):
        """Create all output files needed for downstream GUI components"""
        if not os.path.exists(self.output_folder):
            os.makedirs(self.output_folder)

        # Create tracking DataFrame
        tracking_df = self._create_peak_tracking_table(batch_results)

        if not tracking_df.empty:
            # Main comprehensive tracking file
            tracking_file = os.path.join(self.output_folder, "comprehensive_peak_tracking.csv")
            # Format with 3 decimals for Reference_X/Y, 1 decimal for other columns
            tracking_df_formatted = self._format_dataframe_for_csv(tracking_df)
            tracking_df_formatted.to_csv(tracking_file, index=False)
            log_info(f"Created comprehensive peak tracking: {tracking_file}")

            # Create matrix files
            self._create_intensity_matrix(tracking_df)
            self._create_detected_intensity_matrix(tracking_df)
            self._create_volume_matrix(tracking_df)
            self._create_detection_matrix(tracking_df)

            # Create per-spectrum detailed files
            self._create_per_spectrum_files(batch_results)

            # Create summary statistics file
            self._create_summary_statistics_file(batch_results)

    def _create_peak_tracking_table(self, batch_results: Dict[str, Any]) -> pd.DataFrame:
        """Create comprehensive peak tracking table.

        When extract_delays is enabled, column names use delay values (e.g., "50_Height")
        instead of filenames (e.g., "T1_50ms_Height"). Duplicate delays get unique
        suffixes (e.g., "50_2_Height" for the second spectrum at 50ms).
        """
        all_peaks_data = []

        for peak_idx, ref_peak in self.reference_peaks.iterrows():
            peak_row = {
                'Peak_Number': peak_idx + 1,
                'Assignment': ref_peak.get('Assignment', f'Peak_{peak_idx + 1}'),
                'Reference_X': ref_peak['Position_X'],
                'Reference_Y': ref_peak['Position_Y']
            }

            # Add data for each spectrum
            for spectrum_name, result_data in batch_results['results'].items():
                # Determine column key: use delay column name if extract_delays enabled
                if getattr(self, 'extract_delays', False) and hasattr(self, 'delay_mapping'):
                    column_name = self.delay_mapping.get(spectrum_name)
                    if column_name is not None:
                        spectrum_key = column_name
                    else:
                        spectrum_key = os.path.splitext(spectrum_name)[0]
                else:
                    spectrum_key = os.path.splitext(spectrum_name)[0]
                found_peak = False

                if result_data.get('success', False):
                    integration_results = result_data.get('integration_results', [])
                    for peak in integration_results:
                        # CRITICAL FIX: Convert both sides to string to avoid type mismatch (int vs str)
                        # Bug: peak.get('Assignment') might be "142" (str) while peak_row['Assignment'] is 142 (int)
                        # Additional fix: int64 converts to "142.0" with str(), need to handle properly
                        peak_assignment_raw = peak.get('Assignment', '')
                        row_assignment_raw = peak_row['Assignment']

                        # Normalize both to strings, handling int/float conversions
                        try:
                            # Try to convert to int first to strip any .0 suffix
                            peak_assignment = str(int(float(peak_assignment_raw)))
                        except (ValueError, TypeError):
                            # If conversion fails, use string directly
                            peak_assignment = str(peak_assignment_raw)

                        try:
                            # Same for row assignment
                            row_assignment = str(int(float(row_assignment_raw)))
                        except (ValueError, TypeError):
                            row_assignment = str(row_assignment_raw)

                        # DEBUG: Print matching attempt
                        # print(f"      Matching: peak '{peak_assignment}' vs row '{row_assignment}' -> {peak_assignment == row_assignment}")

                        if peak_assignment == row_assignment:
                            peak_row[f'{spectrum_key}_Detected'] = True
                            peak_row[f'{spectrum_key}_Height'] = peak.get('Height', 0.0)
                            peak_row[f'{spectrum_key}_Volume'] = peak.get('Volume', 0.0)
                            peak_row[f'{spectrum_key}_SNR'] = peak.get('SNR', 0.0)
                            peak_row[f'{spectrum_key}_Quality'] = peak.get('Quality', 'Poor')
                            peak_row[f'{spectrum_key}_Position_X'] = peak.get('Position_X', 0.0)
                            peak_row[f'{spectrum_key}_Position_Y'] = peak.get('Position_Y', 0.0)
                            peak_row[f'{spectrum_key}_R_Squared'] = peak.get('R_Squared', 0.0)

                            # Add detected intensity and linewidth information
                            peak_row[f'{spectrum_key}_Detected_Intensity'] = peak.get('detected_intensity', 0.0)
                            peak_row[f'{spectrum_key}_LW_X'] = peak.get('lw_x', 0.0)
                            peak_row[f'{spectrum_key}_LW_Y'] = peak.get('lw_y', 0.0)

                            found_peak = True
                            break

                if not found_peak:
                    peak_row[f'{spectrum_key}_Detected'] = False
                    peak_row[f'{spectrum_key}_Height'] = 0.0
                    peak_row[f'{spectrum_key}_Volume'] = 0.0
                    peak_row[f'{spectrum_key}_SNR'] = 0.0
                    peak_row[f'{spectrum_key}_Quality'] = 'Not Detected'
                    peak_row[f'{spectrum_key}_Position_X'] = 0.0
                    peak_row[f'{spectrum_key}_Position_Y'] = 0.0
                    peak_row[f'{spectrum_key}_R_Squared'] = 0.0
                    peak_row[f'{spectrum_key}_Detected_Intensity'] = 0.0
                    peak_row[f'{spectrum_key}_LW_X'] = 0.0
                    peak_row[f'{spectrum_key}_LW_Y'] = 0.0

            all_peaks_data.append(peak_row)

        return pd.DataFrame(all_peaks_data)

    def _format_dataframe_for_csv(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Format dataframe columns with specific decimal places:
        - Reference_X, LW_X: 4 decimal places (1H dimension)
        - Reference_Y, LW_Y: 3 decimal places (15N dimension)
        - R_Squared: 3 decimal places
        - Integer columns: preserved as integers (no decimal point)
        - All other numeric columns: 1 decimal place

        Returns a copy with formatted values.
        """
        df_formatted = df.copy()

        for col in df_formatted.columns:
            # Skip non-numeric columns
            if not pd.api.types.is_numeric_dtype(df_formatted[col]):
                continue

            # Check if column is integer type (preserve as integer without decimals)
            if pd.api.types.is_integer_dtype(df_formatted[col]):
                df_formatted[col] = df_formatted[col].apply(lambda x: f'{int(x)}' if pd.notna(x) else '')
            # Format LW_X (1H linewidth) with 4 decimal places
            elif col == 'LW_X' or col.endswith('_LW_X'):
                df_formatted[col] = df_formatted[col].apply(lambda x: f'{x:.4f}' if pd.notna(x) else '')
            # Format positions and LW_Y (15N linewidth) with 3 decimal places
            elif col in ['Reference_X', 'Reference_Y', 'LW_Y'] or col.endswith('_LW_Y'):
                df_formatted[col] = df_formatted[col].apply(lambda x: f'{x:.3f}' if pd.notna(x) else '')
            # Format R_Squared with 3 decimal places
            elif col == 'R_Squared' or col.endswith('_R_Squared'):
                df_formatted[col] = df_formatted[col].apply(lambda x: f'{x:.3f}' if pd.notna(x) else '')
            # Format all other numeric columns with 1 decimal place
            else:
                df_formatted[col] = df_formatted[col].apply(lambda x: f'{x:.1f}' if pd.notna(x) else '')

        return df_formatted

    def _create_intensity_matrix(self, tracking_df: pd.DataFrame):
        """Create peak intensity matrix file"""
        intensity_data = tracking_df[['Peak_Number', 'Assignment', 'Reference_X', 'Reference_Y']].copy()

        for col in tracking_df.columns:
            if col.endswith('_Height'):
                spectrum_name = col.replace('_Height', '')
                intensity_data[spectrum_name] = tracking_df[col]

        intensity_file = os.path.join(self.output_folder, "peak_intensity_matrix.csv")
        # Format with 3 decimals for Reference_X/Y, 1 decimal for other columns
        intensity_data_formatted = self._format_dataframe_for_csv(intensity_data)
        intensity_data_formatted.to_csv(intensity_file, index=False)

    def _create_detection_matrix(self, tracking_df: pd.DataFrame):
        """Create peak detection matrix file"""
        detection_data = tracking_df[['Peak_Number', 'Assignment', 'Reference_X', 'Reference_Y']].copy()

        for col in tracking_df.columns:
            if col.endswith('_Detected'):
                spectrum_name = col.replace('_Detected', '')
                detection_data[spectrum_name] = tracking_df[col].astype(int)

        detection_file = os.path.join(self.output_folder, "peak_detection_matrix.csv")
        # Format Reference_X/Y with 3 decimals (detection columns are integers, keep as is)
        detection_data['Reference_X'] = detection_data['Reference_X'].apply(lambda x: f'{x:.3f}' if pd.notna(x) else '')
        detection_data['Reference_Y'] = detection_data['Reference_Y'].apply(lambda x: f'{x:.3f}' if pd.notna(x) else '')
        detection_data.to_csv(detection_file, index=False)

    def _create_detected_intensity_matrix(self, tracking_df: pd.DataFrame):
        """Create peak detected intensity matrix file (from peak detection, not fitting)"""
        intensity_detected_data = tracking_df[['Peak_Number', 'Assignment', 'Reference_X', 'Reference_Y']].copy()

        for col in tracking_df.columns:
            if col.endswith('_Detected_Intensity'):
                spectrum_name = col.replace('_Detected_Intensity', '')
                intensity_detected_data[spectrum_name] = tracking_df[col]

        detected_intensity_file = os.path.join(self.output_folder, "peak_intensity_detected_matrix.csv")
        # Format with 3 decimals for Reference_X/Y, 1 decimal for other columns
        intensity_detected_data_formatted = self._format_dataframe_for_csv(intensity_detected_data)
        intensity_detected_data_formatted.to_csv(detected_intensity_file, index=False)
        log_info(f"Created detected intensity matrix: {detected_intensity_file}")

    def _create_volume_matrix(self, tracking_df: pd.DataFrame):
        """Create peak volume matrix file"""
        volume_data = tracking_df[['Peak_Number', 'Assignment', 'Reference_X', 'Reference_Y']].copy()

        for col in tracking_df.columns:
            if col.endswith('_Volume'):
                spectrum_name = col.replace('_Volume', '')
                volume_data[spectrum_name] = tracking_df[col]

        volume_file = os.path.join(self.output_folder, "peak_volume_matrix.csv")
        # Format with 3 decimals for Reference_X/Y, 1 decimal for other columns
        volume_data_formatted = self._format_dataframe_for_csv(volume_data)
        volume_data_formatted.to_csv(volume_file, index=False)
        log_info(f"Created volume matrix: {volume_file}")

    def _create_per_spectrum_files(self, batch_results: Dict[str, Any]):
        """
        Create individual CSV file for each spectrum with comprehensive peak information

        Columns: peak_number, assignment, Reference_X, Reference_Y,
                 Detected_Intensity, Height, Volume, LW_X, LW_Y, R_Squared, Quality
        """
        per_spectrum_folder = os.path.join(self.output_folder, "per_spectrum_results")
        if not os.path.exists(per_spectrum_folder):
            os.makedirs(per_spectrum_folder)

        for spectrum_name, result_data in batch_results['results'].items():
            if not result_data.get('success', False):
                continue

            spectrum_key = os.path.splitext(spectrum_name)[0]
            integration_results = result_data.get('integration_results', [])

            if not integration_results:
                continue

            # Build per-spectrum DataFrame
            spectrum_peaks = []
            for peak in integration_results:
                # Find reference coordinates
                ref_x, ref_y = 0.0, 0.0
                assignment = peak.get('assignment', peak.get('Assignment', 'Unknown'))

                if self.reference_peaks is not None:
                    try:
                        matching_ref = self.reference_peaks[
                            self.reference_peaks['Assignment'] == assignment
                        ]
                        if not matching_ref.empty:
                            ref_x = float(matching_ref['Position_X'].iloc[0])
                            ref_y = float(matching_ref['Position_Y'].iloc[0])
                    except Exception:
                        pass

                peak_row = {
                    'Peak_Number': peak.get('peak_number', peak.get('Peak_Number', 0)),
                    'Assignment': assignment,
                    'Reference_X': ref_x,
                    'Reference_Y': ref_y,
                    'Detected_Intensity': peak.get('detected_intensity', 0.0),
                    'Height': peak.get('Height', peak.get('height', 0.0)),
                    'Volume': peak.get('Volume', peak.get('volume', 0.0)),
                    'LW_X': peak.get('lw_x', 0.0),
                    'LW_Y': peak.get('lw_y', 0.0),
                    'R_Squared': peak.get('R_Squared', peak.get('r_squared', 0.0)),
                    'Quality': peak.get('Quality', peak.get('quality', 'Unknown'))
                }
                spectrum_peaks.append(peak_row)

            if spectrum_peaks:
                spectrum_df = pd.DataFrame(spectrum_peaks)
                spectrum_file = os.path.join(per_spectrum_folder, f"{spectrum_key}.csv")
                # Format with 3 decimals for Reference_X/Y, 1 decimal for other columns
                spectrum_df_formatted = self._format_dataframe_for_csv(spectrum_df)
                spectrum_df_formatted.to_csv(spectrum_file, index=False)

        log_info(f"Created per-spectrum CSV files in: {per_spectrum_folder}")

    def _create_summary_statistics_file(self, batch_results: Dict[str, Any]):
        """Create summary statistics file"""
        summary_data = []

        for spectrum_name, result_data in batch_results['results'].items():
            summary_row = {
                'Spectrum': spectrum_name,
                'Success': result_data.get('success', False),
                'Total_Peaks': result_data.get('total_peaks', 0),
                'Successful_Fits': result_data.get('successful_fits', 0),
                'Success_Rate': result_data.get('success_rate', 0.0)
            }
            summary_data.append(summary_row)

        summary_df = pd.DataFrame(summary_data)
        summary_file = os.path.join(self.output_folder, "processing_summary.csv")
        summary_df.to_csv(summary_file, index=False, float_format='%.2f')

    def cancel_processing(self):
        """Cancel current processing operation"""
        self.processing_active = False
        log_info("Multi-spectrum processing cancelled")

    def _extract_and_store_linewidth(self, result, assignment):
        """
        Extract linewidths from reference spectrum fit results.

        Handles both 1D (x_fit/y_fit subdicts) and 2D PS2D (top-level sigma/gamma) formats.
        Stores in PS2D-compatible format (lw_lor_f1, lw_gau_f1, etc.) for reuse.

        Parameters
        ----------
        result : dict
            Fitting result from reference spectrum
        assignment : str
            Peak assignment identifier
        """
        try:
            # Check for top-level sigma/gamma fields (from PS2D 2D fit results)
            if 'sigma_x' in result or 'sigma_y' in result:
                # PS2D top-level format: sigma_x = 1H Gaussian, gamma_x = 1H Lorentzian
                self.reference_linewidths[assignment] = {
                    'lw_gau_f2': result.get('sigma_x', 0.015),   # 1H Gaussian (F2)
                    'lw_lor_f2': result.get('gamma_x', 0.001),   # 1H Lorentzian (F2)
                    'lw_gau_f1': result.get('sigma_y', 0.3),     # 15N Gaussian (F1)
                    'lw_lor_f1': result.get('gamma_y', 0.0001),  # 15N Lorentzian (F1)
                    'r_squared': result.get('r_squared', result.get('avg_r_squared', 0))
                }
            else:
                # 1D format with x_fit/y_fit subdicts
                x_fit = result.get('x_fit', {})
                y_fit = result.get('y_fit', {})

                self.reference_linewidths[assignment] = {
                    'lw_gau_f2': x_fit.get('sigma', 0.015),   # 1H Gaussian (F2)
                    'lw_lor_f2': x_fit.get('gamma', 0.001),   # 1H Lorentzian (F2)
                    'lw_gau_f1': y_fit.get('sigma', 0.3),     # 15N Gaussian (F1)
                    'lw_lor_f1': y_fit.get('gamma', 0.0001),  # 15N Lorentzian (F1)
                    'r_squared': max(x_fit.get('r_squared', 0), y_fit.get('r_squared', 0))
                }

            ref_lw = self.reference_linewidths[assignment]
            log_info(f"Stored linewidths for {assignment}: "
                     f"F1(Gau={ref_lw['lw_gau_f1']:.4f}, Lor={ref_lw['lw_lor_f1']:.5f}), "
                     f"F2(Gau={ref_lw['lw_gau_f2']:.5f}, Lor={ref_lw['lw_lor_f2']:.6f})")

        except Exception as e:
            log_warning(f"Could not extract linewidths for {assignment}: {e}")

    def _init_training_collector(self, output_folder: str):
        """Initialize the ML training data collector with series output folder."""
        try:
            from lunaNMR.ml.comprehensive_collector import ComprehensiveTrainingCollector
            from pathlib import Path

            # Use series_results folder for data collection
            storage_dir = Path(output_folder) / "data_collection"

            self.training_collector = ComprehensiveTrainingCollector(
                storage_dir=storage_dir,
                min_r2=0.70,
                auto_save=True,
            )
            log_info(f"ML training data collector initialized: {storage_dir}")
        except Exception as e:
            log_warning(f"Could not initialize ML training collector: {e}")
            self.training_collector = None
            self.collect_training_data = False

    def _collect_spectrum_training_data(
        self,
        spectrum_name: str,
        spectrum_number: int,
        fitted_results: List[Dict],
        learned_statistics: Optional[Dict],
        total_spectra: int,
    ):
        """Collect training data for a single spectrum in series processing.

        Parameters
        ----------
        spectrum_name : str
            Name of the spectrum file
        spectrum_number : int
            Position in series (1-indexed)
        fitted_results : list
            Fit results for all peaks
        learned_statistics : dict, optional
            PASS1 statistics
        total_spectra : int
            Total number of spectra in series
        """
        if self.training_collector is None:
            return

        try:
            # Get spectrum data from integrator
            spectrum_data = self.integrator.nmr_data
            if spectrum_data is None:
                return

            # Get peak list
            peak_list = self.integrator.peak_list
            if peak_list is None or len(peak_list) == 0:
                return
            if hasattr(peak_list, 'to_dict'):
                peak_list_dicts = peak_list.to_dict('records')
            else:
                peak_list_dicts = peak_list

            # Enhance learned_statistics with additional stats from fit results
            if learned_statistics is None:
                learned_statistics = {}

            # Calculate additional statistics from fit results if not present
            if fitted_results and 'lw_f1_std' not in learned_statistics:
                lw_f1_values = []
                lw_f2_values = []
                r2_values = []
                for r in fitted_results:
                    if r.get('success', False):
                        # Get linewidths from fit results
                        # PS2D 2D fits: lw_lor_f1/lw_gau_f1 or gamma_y/sigma_y
                        # 1D fits: linewidths in nested y_fit/x_fit dicts
                        y_fit = r.get('y_fit', {})
                        x_fit = r.get('x_fit', {})

                        # F1 linewidth (Y dimension)
                        lw_lor_f1 = r.get('lw_lor_f1', r.get('gamma_y', y_fit.get('gamma', 0)))
                        lw_gau_f1 = r.get('lw_gau_f1', r.get('sigma_y', y_fit.get('sigma', 0)))
                        lw_f1 = lw_lor_f1 + lw_gau_f1

                        # F2 linewidth (X dimension)
                        lw_lor_f2 = r.get('lw_lor_f2', r.get('gamma_x', r.get('gamma', x_fit.get('gamma', 0))))
                        lw_gau_f2 = r.get('lw_gau_f2', r.get('sigma_x', r.get('sigma', x_fit.get('sigma', 0))))
                        lw_f2 = lw_lor_f2 + lw_gau_f2

                        if lw_f1 > 0:
                            lw_f1_values.append(lw_f1)
                        if lw_f2 > 0:
                            lw_f2_values.append(lw_f2)
                        r2 = r.get('r_squared', r.get('avg_r_squared', 0))
                        if r2 > 0:
                            r2_values.append(r2)

                if lw_f1_values:
                    learned_statistics['lw_f1_std'] = float(np.std(lw_f1_values))
                    learned_statistics['lw_f1_min'] = float(np.min(lw_f1_values))
                    learned_statistics['lw_f1_max'] = float(np.max(lw_f1_values))
                if lw_f2_values:
                    learned_statistics['lw_f2_std'] = float(np.std(lw_f2_values))
                    learned_statistics['lw_f2_min'] = float(np.min(lw_f2_values))
                    learned_statistics['lw_f2_max'] = float(np.max(lw_f2_values))
                if r2_values:
                    learned_statistics['n_good_peaks'] = len(r2_values)
                    learned_statistics['mean_r_squared'] = float(np.mean(r2_values))

            # Get PPM ranges from integrator's ppm axes
            ppm_ranges = {}
            if hasattr(self.integrator, 'ppm_y_axis') and self.integrator.ppm_y_axis is not None:
                ppm_y = self.integrator.ppm_y_axis
                ppm_ranges['f1'] = (float(min(ppm_y)), float(max(ppm_y)))
            if hasattr(self.integrator, 'ppm_x_axis') and self.integrator.ppm_x_axis is not None:
                ppm_x = self.integrator.ppm_x_axis
                ppm_ranges['f2'] = (float(min(ppm_x)), float(max(ppm_x)))

            # Detect nucleus type from F1 range
            nucleus_type = "15N"
            if ppm_ranges.get('f1'):
                f1_max = ppm_ranges['f1'][1]
                if f1_max > 150:
                    nucleus_type = "13C"

            # Get field strength from metadata
            field_strength = 600.0
            if hasattr(self.integrator, 'dic') and self.integrator.dic:
                field_strength = float(
                    self.integrator.dic.get('FDF2OBS',
                    self.integrator.dic.get('SFO1',
                    self.integrator.dic.get('bf1', 600.0)))
                )

            # Get intermediate results from parallel processor (if available)
            adaptive_results = None
            timing_info = None
            cluster_info = None
            pass1_results = None
            pass1bis_results = None
            pass2_results = None

            if hasattr(self.integrator, 'enhanced_fitter'):
                fitter = self.integrator.enhanced_fitter
                if hasattr(fitter, 'parallel_processor') and fitter.parallel_processor:
                    pp = fitter.parallel_processor
                    if hasattr(pp, 'optimal_params') and pp.optimal_params:
                        adaptive_results = pp.optimal_params
                    if hasattr(pp, 'timing_info') and pp.timing_info:
                        timing_info = pp.timing_info
                    if hasattr(pp, 'cluster_info') and pp.cluster_info:
                        cluster_info = pp.cluster_info
                    if hasattr(pp, 'pass1_results'):
                        pass1_results = pp.pass1_results
                    if hasattr(pp, 'pass1bis_results'):
                        pass1bis_results = pp.pass1bis_results
                    if hasattr(pp, 'pass2_results'):
                        pass2_results = pp.pass2_results

            # Get detection parameters and statistics from integrator
            detection_params = {
                'search_window_x': getattr(self.integrator, 'search_window_x', 0.01),
                'search_window_y': getattr(self.integrator, 'search_window_y', 0.05),
                'noise_threshold': getattr(self.integrator, 'threshold_multiplier', 3.0),
                'noise_level': getattr(self.integrator, 'noise_level', 0),
            }
            detection_statistics = None
            if hasattr(self.integrator, 'get_detection_statistics'):
                detection_statistics = self.integrator.get_detection_statistics()

            # Collect training data
            from pathlib import Path
            success = self.training_collector.collect_spectrum_processing(
                spectrum_data=spectrum_data,
                peak_list=peak_list_dicts,
                fit_results=fitted_results,
                pass1_statistics=learned_statistics,
                adaptive_results=adaptive_results,
                pass1_results=pass1_results,
                pass1bis_results=pass1bis_results,
                pass2_results=pass2_results,
                cluster_info=cluster_info,
                timing_info=timing_info,
                spectrum_name=Path(spectrum_name).name,
                nucleus_type=nucleus_type,
                field_strength_mhz=field_strength,
                ppm_ranges=ppm_ranges,
                processing_mode='parallel',
                is_reference=(spectrum_number == 1),
                series_index=spectrum_number - 1,
                detection_params=detection_params,
                detection_statistics=detection_statistics,
            )

            if success:
                log_info(f"Collected ML training data: {Path(spectrum_name).name} ({spectrum_number}/{total_spectra})")

        except Exception as e:
            log_warning(f"Failed to collect training data for {spectrum_name}: {e}")

    def _sample_intensities_for_peak_list(self):
        """
        Sample intensities from current spectrum at each peak position.

        This is critical for cascade mode spectrum N+1: positions come from
        spectrum N-1, but intensities must come from the CURRENT spectrum
        for proper fitter initialization.

        Without this, T1/T2 experiments fail because initial intensity guess
        (from spectrum 1) is very different from actual intensity (spectrum N).
        """
        if self.integrator.nmr_data is None:
            log_warning("Cannot sample intensities: nmr_data is None")
            return

        if self.integrator.ppm_x_axis is None or self.integrator.ppm_y_axis is None:
            log_warning("Cannot sample intensities: ppm axes not available")
            return

        ppm_x = self.integrator.ppm_x_axis
        ppm_y = self.integrator.ppm_y_axis
        nmr_data = self.integrator.nmr_data

        sampled_count = 0
        for idx, row in self.integrator.peak_list.iterrows():
            try:
                pos_x = row['Position_X']
                pos_y = row['Position_Y']

                # Convert ppm to array indices
                x_idx = np.argmin(np.abs(ppm_x - pos_x))
                y_idx = np.argmin(np.abs(ppm_y - pos_y))

                # Ensure indices are within bounds
                x_idx = np.clip(x_idx, 0, nmr_data.shape[1] - 1)
                y_idx = np.clip(y_idx, 0, nmr_data.shape[0] - 1)

                # Sample intensity at peak position
                intensity = float(nmr_data[y_idx, x_idx])

                # Update peak_list with sampled intensity
                self.integrator.peak_list.at[idx, 'Height'] = intensity
                self.integrator.peak_list.at[idx, 'Intensity'] = intensity

                sampled_count += 1

            except Exception as e:
                assignment = row.get('Assignment', 'Unknown')
                log_warning(f"Could not sample intensity for {assignment}: {e}")


## GM added

class BatchResults:
    """Container for batch processing results with analysis capabilities"""

    def __init__(self):
        self.results = {}
        self.metadata = {
            'start_time': None,
            'end_time': None,
            'total_spectra': 0,
            'successful_spectra': 0,
            'failed_spectra': 0,
            'processing_mode': 'unknown'
        }
        self.statistics = {}
        self.errors = []

    def add_result(self, spectrum_name, result_data):
        """Add a single spectrum result"""
        self.results[spectrum_name] = result_data

        # Update metadata
        if result_data.get('status') == 'success':
            self.metadata['successful_spectra'] += 1
        else:
            self.metadata['failed_spectra'] += 1
            if 'error' in result_data:
                self.errors.append({
                    'spectrum': spectrum_name,
                    'error': result_data['error']
                })

    def get_summary(self):
        """Get processing summary"""
        total = len(self.results)
        success = self.metadata['successful_spectra']
        failed = self.metadata['failed_spectra']

        if self.metadata['start_time'] and self.metadata['end_time']:
            duration = self.metadata['end_time'] - self.metadata['start_time']
            duration_str = str(duration).split('.')[0]  # Remove microseconds
        else:
            duration_str = "Unknown"

        summary = {
            'total_spectra': total,
            'successful': success,
            'failed': failed,
            'success_rate': (success / total * 100) if total > 0 else 0,
            'duration': duration_str,
            'processing_mode': self.metadata['processing_mode'],
            'error_count': len(self.errors)
        }

        return summary

    def export_results(self, output_folder):
        """Export all results to files"""
        if not os.path.exists(output_folder):
            os.makedirs(output_folder)

        # Export summary
        summary_file = os.path.join(output_folder, "batch_summary.json")
        with open(summary_file, 'w') as f:
            json.dump({
                'metadata': self.metadata,
                'summary': self.get_summary(),
                'statistics': self.statistics,
                'errors': self.errors
            }, f, indent=2, default=str)

        # Export individual results
        if self.results:
            # Create consolidated integration results
            all_integrations = []
            detection_stats = []

            for spectrum_name, result in self.results.items():
                # Detection statistics
                detection_stats.append({
                    'spectrum': spectrum_name,
                    'status': result.get('status', 'unknown'),
                    'detected_peaks': result.get('detected_peaks', 0),
                    'total_peaks': result.get('total_peaks', 0),
                    'detection_rate': result.get('detection_rate', 0.0),
                    'noise_level': result.get('noise_level', 0.0),
                    'processing_time': result.get('processing_time', 0.0)
                })

                # Individual integrations
                if result.get('integration_results'):
                    for integration in result['integration_results']:
                        integration_with_spectrum = integration.copy()
                        integration_with_spectrum['spectrum'] = spectrum_name
                        all_integrations.append(integration_with_spectrum)

            # Save detection statistics
            if detection_stats:
                stats_file = os.path.join(output_folder, "detection_statistics.csv")
                pd.DataFrame(detection_stats).to_csv(stats_file, index=False)

            # Save consolidated integrations
            if all_integrations:
                integration_file = os.path.join(output_folder, "all_integrations.csv")
                integration_df = pd.DataFrame(all_integrations)
                # Format ppm_x and ppm_y with 3 decimals (position data), other numeric columns with 1 decimal
                integration_df_formatted = integration_df.copy()
                for col in integration_df_formatted.columns:
                    if not pd.api.types.is_numeric_dtype(integration_df_formatted[col]):
                        continue
                    # Position and linewidth columns get 3 decimal places
                    if col in ['ppm_x', 'ppm_y', 'lw_x', 'lw_y']:
                        integration_df_formatted[col] = integration_df_formatted[col].apply(lambda x: f'{x:.3f}' if pd.notna(x) else '')
                    # All other numeric columns get 1 decimal place
                    else:
                        integration_df_formatted[col] = integration_df_formatted[col].apply(lambda x: f'{x:.1f}' if pd.notna(x) else '')
                integration_df_formatted.to_csv(integration_file, index=False)

        return output_folder

    def to_tidy_dataframe(self) -> pd.DataFrame:
        """Generate tidy DataFrame from BatchResults for T1/T2 fitting.

        This method extracts the same data structure that would be in
        series_analysis_tidy.csv, allowing DynamiXs to consume results
        directly from memory without requiring a CSV file.

        Returns:
            DataFrame with columns: spectrum_name, assignment, peak_number,
            ppm_x, ppm_y, height, volume, snr, quality, r_squared
        """
        tidy_data = []

        for spectrum_name, result_data in self.results.items():
            # Skip failed spectra
            if not result_data.get('success', False) and result_data.get('status') != 'success':
                continue

            # Use filename without extension as spectrum header
            spectrum_header = os.path.splitext(spectrum_name)[0]

            # Extract integration results
            integration_results = result_data.get('integration_results', [])
            for peak in integration_results:
                row = {
                    'spectrum_name': spectrum_header,
                    'assignment': peak.get('assignment', ''),
                    'peak_number': peak.get('peak_number', 0),
                    'ppm_x': peak.get('ppm_x', 0.0),
                    'ppm_y': peak.get('ppm_y', 0.0),
                    'height': peak.get('height', 0.0),
                    'volume': peak.get('volume', 0.0),
                    'snr': peak.get('snr', 0.0),
                    'quality': peak.get('quality', 'Unknown'),
                    'r_squared': peak.get('r_squared', 0.0)
                }
                tidy_data.append(row)

        if not tidy_data:
            return pd.DataFrame(columns=['spectrum_name', 'assignment', 'peak_number',
                                        'ppm_x', 'ppm_y', 'height', 'volume',
                                        'snr', 'quality', 'r_squared'])

        return pd.DataFrame(tidy_data)

    def to_fitting_dataframe(self, value_column: str = 'volume') -> pd.DataFrame:
        """Generate pivot DataFrame for T1/T2 fitting analysis.

        Creates the wide format expected by the fitting module:
        - Row per peak (identified by assignment)
        - Column per spectrum (numeric delay time)
        - Values are volumes (or heights)

        The fitting module expects numeric column names (e.g., "10", "50", "100")
        representing delay times in ms. This method extracts delays from spectrum
        names like "T1_10ms" -> "10".

        Args:
            value_column: Which value to use for fitting ('volume' or 'height')

        Returns:
            DataFrame with columns: Assignment, <delay1>, <delay2>, ...
            Suitable for direct use with dynamiXs_T1_T2 fitting module.
        """
        from lunaNMR.utils.delay_extractor import DelayExtractor

        # First get the tidy DataFrame
        tidy_df = self.to_tidy_dataframe()

        if tidy_df.empty:
            return pd.DataFrame(columns=['Assignment'])

        # Extract values from spectrum names and create mapping
        extractor = DelayExtractor(mode=getattr(self, 'series_mode', 'time'))
        spectrum_names = tidy_df['spectrum_name'].unique()

        # Build mapping: spectrum_name -> numeric delay column name
        delay_mapping = {}
        delay_counts = {}  # Track duplicate delays

        for spec_name in spectrum_names:
            delay_ms = extractor.extract_value(spec_name)
            if delay_ms is not None:
                # Format as integer if whole number, else with decimal
                delay_str = str(int(delay_ms)) if delay_ms == int(delay_ms) else str(delay_ms)

                # Handle duplicate delays (e.g., two spectra at 50ms)
                if delay_str in delay_counts:
                    delay_counts[delay_str] += 1
                    delay_mapping[spec_name] = f"{delay_str}_{delay_counts[delay_str]}"
                else:
                    delay_counts[delay_str] = 1
                    delay_mapping[spec_name] = delay_str
            else:
                # Fallback: use spectrum name as-is (may cause fitting issues)
                delay_mapping[spec_name] = spec_name

        # Apply mapping to spectrum_name column
        tidy_df = tidy_df.copy()
        tidy_df['delay_column'] = tidy_df['spectrum_name'].map(delay_mapping)

        # Pivot: rows=assignment, columns=delay_column, values=volume or height
        pivot_df = tidy_df.pivot_table(
            index='assignment',
            columns='delay_column',
            values=value_column,
            aggfunc='first'  # In case of duplicates, take first
        )

        # Reset index to make Assignment a column
        pivot_df = pivot_df.reset_index()
        pivot_df = pivot_df.rename(columns={'assignment': 'Assignment'})

        # Sort columns: Assignment first, then numeric delay columns in order
        cols = pivot_df.columns.tolist()
        delay_cols = [c for c in cols if c != 'Assignment']

        # Sort delay columns numerically
        def sort_key(col):
            try:
                # Handle "50" or "50_2" format
                base = col.split('_')[0]
                return float(base)
            except (ValueError, IndexError):
                return float('inf')

        delay_cols.sort(key=sort_key)
        pivot_df = pivot_df[['Assignment'] + delay_cols]

        return pivot_df
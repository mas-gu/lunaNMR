# ABOUTME: Cluster-based parallel processing distributor using multiprocessing.Pool for 2.7× speedup.
# ABOUTME: Distributes entire overlap clusters (not individual peaks) across workers with deterministic clustering.
#!/usr/bin/env python3
"""
Complete Parallel Voigt Fitting Implementation

This module provides comprehensive parallel processing for the entire
Voigt fitting pipeline while maintaining complete compatibility with
the existing EnhancedVoigtFitter workflow.

Key Features:
- Complete parallel processing of ALL Voigt fitting steps
- Baseline correction (ArPLS) parallelization
- Parameter estimation parallelization  
- Voigt profile optimization parallelization
- Quality assessment parallelization
- Shared memory for efficient data transfer
- Maintains 100% compatibility with existing interfaces

Author: Guillaume Mas
Date: 2025
"""

import multiprocessing as mp
from multiprocessing import Pool, cpu_count, shared_memory
import numpy as np
import time
import traceback
import sys
import os


from lunaNMR.utils.output_manager import log_progress, log_info, log_warning, log_error

# Adaptive optimization imports
from lunaNMR.core.adaptive_optimizer import (
    AdaptiveOptimizer,
    estimate_noise_level,
    create_series_params,
    check_lw_deviation
)

# CRITICAL FIX for Linux deadlock: Set multiprocessing start method to 'spawn'
# On Linux, default 'fork' causes pipe buffer deadlock when workers print verbose output
# On Mac, 'spawn' is already the default, so this has no effect there
try:
    mp.set_start_method('spawn', force=False)
except RuntimeError:
    # Already set (can only be called once)
    current_method = mp.get_start_method()

class ParallelVoigtProcessor:
    """
    Master coordinator for parallel Voigt fitting that maintains 
    complete compatibility with current EnhancedVoigtFitter logic.
    
    This class orchestrates the complete parallel execution of:
    - Spectral region extraction
    - Baseline correction (ArPLS)
    - Initial parameter estimation  
    - Multi-peak detection
    - Voigt profile optimization
    - Quality assessment
    - Result compilation
    """
    
    def __init__(self, enhanced_voigt_fitter, max_workers=None, enable_adaptive=True):
        """
        Initialize parallel Voigt processor

        Args:
            enhanced_voigt_fitter: Original EnhancedVoigtFitter instance
            max_workers: Maximum number of worker processes (default: 75% of cores)
            enable_adaptive: Enable adaptive parameter optimization (default: True)
        """
        self.original_fitter = enhanced_voigt_fitter
        self.shared_memory_blocks = []  # Track for cleanup
        self.enable_adaptive = enable_adaptive

        # Adaptive optimization state
        self.optimal_params = None  # Stores grid search optimization results
        self.series_params = None   # Locked params for series integration

        # Intermediate results for ML training data collection
        self.timing_info = {}       # Timing for each processing phase
        self.cluster_info = {}      # Clustering information
        self.pass1_results = None   # PASS1 (isolated) results
        self.pass1bis_results = None  # PASS1-bis (refit) results
        self.pass2_results = None   # PASS2 (multi-peak) results

        # Configure worker processes
        if max_workers is None:
            self.max_workers = max(1, int(cpu_count() * 0.75))
        else:
            self.max_workers = max_workers
            
    def set_series_params(self, series_params):
        """
        Set locked series parameters from reference spectrum.

        Used for series integration to ensure consistent parameters
        and enable linewidth deviation detection.

        Args:
            series_params: Dict with radF1, radF2, overlap thresholds,
                          cluster_assignments, reference_LW values
        """
        self.series_params = series_params
        if series_params:
            radF1 = series_params.get('radF1')
            radF2 = series_params.get('radF2')
            radF1_str = f"{radF1:.4f}" if radF1 is not None else "N/A"
            radF2_str = f"{radF2:.5f}" if radF2 is not None else "N/A"
            log_info(f"Series parameters loaded from reference: radF1={radF1_str}, radF2={radF2_str}")

    def fit_all_peaks_parallel(self, peak_list, progress_callback=None, locked_clusters_by_assignment=None,
                                pre_learned_statistics=None, reference_linewidths=None):
        """
        Complete parallel workflow using IDENTICAL clustering algorithm as sequential mode.

        CRITICAL: This now uses identify_overlap_clusters() ONCE before distributing tasks,
        ensuring identical results to sequential mode (just faster with multiple cores).

        Args:
            peak_list: DataFrame with peak information
            progress_callback: Optional progress callback function
            locked_clusters_by_assignment: Optional pre-computed clusters from reference spectrum.
                If provided, these clusters are used instead of recalculating.
                Format: List of lists of assignment strings, e.g. [['Peak1', 'Peak2'], ['Peak3']]
            pre_learned_statistics: Optional linewidth statistics from reference spectrum.
                If provided, PASS 1 learning is SKIPPED and these values are used as initial
                guesses for linewidth bounds. This speeds up fitting for spectrum 2+ in series.
                Format: {'lw_f1_median': float, 'lw_f2_median': float, 'lw_f1_mad': float,
                         'lw_f2_mad': float, 'n_samples': int, 'alpha': float}
            reference_linewidths: Optional dict mapping assignment -> {lw_lor_f1, lw_gau_f1, etc.}
                for per-peak linewidth reuse. If provided, linewidths are fixed to these values.

        Returns:
            Tuple of (fitted_results, learned_statistics) where:
                - fitted_results: List of fitting results (same format as sequential)
                - learned_statistics: Dict with linewidth statistics learned from PASS 1,
                  or None if pre_learned_statistics was provided (PASS 1 skipped)
        """
        if len(peak_list) == 0:
            return ([], None)

        # Store reference linewidths for use in shared context
        self._reference_linewidths = reference_linewidths

        log_progress(f"Starting cluster-based parallel fitting of {len(peak_list)} peaks")
        if reference_linewidths:
            log_info(f"Per-peak linewidth reuse enabled: {len(reference_linewidths)} reference linewidths")
        start_time = time.time()

        # Reset intermediate results for ML training data collection
        self.timing_info = {'total': 0, 'pass1': 0, 'adaptive': 0, 'pass1bis': 0, 'pass2': 0}
        self.cluster_info = {}
        self.pass1_results = None
        self.pass1bis_results = None
        self.pass2_results = None

        try:
            # STEP 1: Build all_peaks_context (identical to sequential mode)
            # Include detection fields from fitted_peaks if available
            fitted_peaks_by_pos = {}
            if hasattr(self.original_fitter, 'parent') and hasattr(self.original_fitter.parent, 'fitted_peaks'):
                for fp in self.original_fitter.parent.fitted_peaks:
                    key = (fp.get('ppm_x', 0), fp.get('ppm_y', 0))
                    fitted_peaks_by_pos[key] = fp

            all_peaks_context = []
            for _, row in peak_list.iterrows():
                intensity = row.get('Height', row.get('Intensity', None))
                x_ppm = float(row['Position_X'])
                y_ppm = float(row['Position_Y'])

                peak_ctx = {
                    'assignment': str(row.get('Assignment', 'Unknown')),
                    'x_ppm': x_ppm,
                    'y_ppm': y_ppm,
                    'pos_x': x_ppm,
                    'pos_y': y_ppm,
                    'intensity': intensity
                }

                # Merge detection fields from fitted_peaks if available
                # Match by position with tolerance (detected peaks can drift from reference)
                for (fp_x, fp_y), fp_data in fitted_peaks_by_pos.items():
                    if abs(fp_x - x_ppm) < 0.001 and abs(fp_y - y_ppm) < 0.01:
                        peak_ctx['detected'] = fp_data.get('detected', True)
                        peak_ctx['detection_quality'] = fp_data.get('detection_quality', 'Matched')
                        peak_ctx['distance_from_reference'] = fp_data.get('distance_from_reference', 0)
                        peak_ctx['distance_from_reference_x'] = fp_data.get('distance_from_reference_x', 0)
                        peak_ctx['distance_from_reference_y'] = fp_data.get('distance_from_reference_y', 0)
                        peak_ctx['distance_from_reference_elliptical'] = fp_data.get('distance_from_reference_elliptical', 0)
                        peak_ctx['reference_retained'] = fp_data.get('reference_retained', False)
                        # Copy detected intensity if not already in peak_list
                        if peak_ctx['intensity'] is None:
                            peak_ctx['intensity'] = fp_data.get('intensity')
                        break

                all_peaks_context.append(peak_ctx)
            log_info(f"2D overlap detection enabled with {len(all_peaks_context)} peaks context")

            # STEP 2: Get clusters - either locked or computed
            # Priority: explicit parameter > series_params > fresh calculation
            effective_locked_clusters = locked_clusters_by_assignment

            # Fallback to series_params if explicit parameter not provided
            if effective_locked_clusters is None and self.series_params:
                effective_locked_clusters = self.series_params.get('locked_clusters_by_assignment')
                if effective_locked_clusters:
                    log_info(f"Using locked clusters from series_params ({len(effective_locked_clusters)} clusters)")

            if effective_locked_clusters is not None:
                # Use locked clusters from reference spectrum
                log_info(f"Using LOCKED clusters from reference spectrum ({len(effective_locked_clusters)} clusters)")
                clusters = self._convert_locked_clusters_to_positions(
                    effective_locked_clusters, all_peaks_context
                )
                # Apply locked thresholds from series_params if available
                if self.series_params:
                    from lunaNMR.core.ps2d_config import get_ps2d_config
                    config = get_ps2d_config()
                    # Only apply if keys exist (minimal series_params may not have these)
                    if 'radF1' in self.series_params:
                        config.radF1 = self.series_params['radF1']
                    if 'radF2' in self.series_params:
                        config.radF2 = self.series_params['radF2']
                    if 'overlap_threshold_x' in self.series_params:
                        config.overlap_threshold_x = self.series_params['overlap_threshold_x']
                    if 'overlap_threshold_y' in self.series_params:
                        config.overlap_threshold_y = self.series_params['overlap_threshold_y']
            else:
                # First spectrum - calculate fresh clusters
                clusters = self.original_fitter.parent.identify_overlap_clusters(all_peaks_context)
                log_info(f"Computed {len(clusters)} fresh clusters")

            # Build peak metadata mapping (peak position → metadata)
            peak_metadata = {}
            for i, (peak_idx, peak_row) in enumerate(peak_list.iterrows()):
                peak_x = float(peak_row['Position_X'])
                peak_y = float(peak_row['Position_Y'])
                assignment = peak_row.get('Assignment', f'Peak_{i+1}')
                peak_metadata[(peak_x, peak_y)] = {
                    'peak_number': i + 1,
                    'assignment': assignment,
                    'list_index': i
                }

            # STEP 3: Separate clusters into isolated peaks vs multi-peak clusters
            isolated_clusters = [c for c in clusters if len(c) == 1]
            multi_peak_clusters = [c for c in clusters if len(c) > 1]

            log_info(f"Cluster breakdown: {len(isolated_clusters)} isolated, {len(multi_peak_clusters)} multi-peak")

            # ========== PASS 1: Fit isolated peaks to collect linewidth statistics ==========
            # Track whether we learned new statistics (for return value)
            learned_statistics = None  # Will be set if we learn new stats

            # Check if pre-learned statistics are provided (from reference spectrum)
            if pre_learned_statistics is not None:
                # USE PRE-LEARNED STATISTICS: Skip PASS 1 learning entirely
                log_info("USING PRE-LEARNED STATISTICS from reference spectrum (skipping PASS 1 learning)")
                log_info(f"LW F1: {pre_learned_statistics['lw_f1_median']:.4f} ± {pre_learned_statistics['lw_f1_mad']:.4f} ppm")
                log_info(f"LW F2: {pre_learned_statistics['lw_f2_median']:.5f} ± {pre_learned_statistics['lw_f2_mad']:.5f} ppm")
                log_info(f"(From {pre_learned_statistics['n_samples']} fits in reference)")

                spectrum_statistics = pre_learned_statistics
                isolated_results = []  # No PASS 1 results - will fit isolated peaks in PASS 2

                # For pre-learned mode, all peaks go to PASS 2 (no PASS 1 or PASS 1-bis)
                # Isolated peaks will be fitted with pre-learned statistics as priors
            else:
                # LEARN NEW STATISTICS: Run full PASS 1 cycle
                spectrum_statistics = None

                if len(isolated_clusters) > 0:
                    log_progress(f"PASS 1: Fitting {len(isolated_clusters)} isolated peaks to learn linewidth statistics...")
                    pass1_start = time.time()

                    # Prepare shared context WITHOUT spectrum_statistics (Pass 1)
                    shared_context_pass1 = self._prepare_shared_context(all_peaks_context, spectrum_statistics=None,
                                                                        reference_linewidths=self._reference_linewidths)

                    # Create tasks for isolated peaks only
                    isolated_tasks = self._create_cluster_tasks(isolated_clusters, peak_metadata, peak_list, shared_context_pass1)

                    # Execute Pass 1 fitting in parallel
                    isolated_results = self._execute_parallel_cluster_fitting(isolated_tasks, progress_callback, pass_name="PASS 1")

                    # Collect spectrum-wide linewidth statistics from Pass 1 results
                    spectrum_statistics = self._collect_linewidth_statistics(isolated_results)

                    pass1_time = time.time() - pass1_start
                    self.timing_info['pass1'] = pass1_time
                    self.pass1_results = isolated_results  # Store PASS1 results

                    if spectrum_statistics:
                        log_info(f"PASS 1 completed in {pass1_time:.1f}s")
                        log_info(f"Learned statistics from {spectrum_statistics['n_samples']} good fits:")
                        log_info(f"  F1: {spectrum_statistics['lw_f1_median']:.4f} ± {spectrum_statistics['lw_f1_mad']:.4f} ppm")
                        log_info(f"  F2: {spectrum_statistics['lw_f2_median']:.5f} ± {spectrum_statistics['lw_f2_mad']:.5f} ppm")
                        # Store for return (this is what we learned)
                        learned_statistics = spectrum_statistics.copy()
                    else:
                        log_warning("PASS 1 completed but insufficient data for statistics (fallback to config)")
                else:
                    log_warning("No isolated peaks found - skipping Pass 1, using config defaults")
                    isolated_results = []

            # ========== DEVIATION CHECK (for series integration) ==========
            if self.series_params and spectrum_statistics and spectrum_statistics.get('n_samples', 0) >= 5:
                has_deviation, warning_msg = check_lw_deviation(
                    current_lw_f1=spectrum_statistics['lw_f1_median'],
                    current_lw_f2=spectrum_statistics['lw_f2_median'],
                    series_params=self.series_params
                )
                if has_deviation:
                    log_warning(f"LINEWIDTH DEVIATION WARNING: {warning_msg}")

            # ========== ADAPTIVE OPTIMIZATION (if enabled and not already have series_params) ==========
            optimal_params = None
            original_multi_peak_clusters = multi_peak_clusters  # Save for re-clustering

            # Only run optimization on reference spectrum (when series_params is None)
            # For subsequent spectra in series, use locked parameters from reference
            if self.enable_adaptive and not self.series_params and spectrum_statistics and spectrum_statistics.get('n_samples', 0) >= 10:
                log_progress("ADAPTIVE OPTIMIZATION: Learning optimal parameters...")
                adaptive_start = time.time()

                # Estimate noise level from spectrum
                noise_level = estimate_noise_level(self.original_fitter.parent.nmr_data)

                # Get PS2D config for bounds
                from lunaNMR.core.ps2d_config import get_ps2d_config
                config = get_ps2d_config()

                # Create optimizer
                optimizer = AdaptiveOptimizer(noise_level=noise_level, config=config)

                # Extract good isolated peaks for optimization
                good_isolated_peaks = self._extract_good_peaks_for_optimization(isolated_results)

                if len(good_isolated_peaks) >= 10:
                    # Prepare shared context for parallel adaptive optimization
                    # This includes spectrum_statistics for proper worker initialization
                    shared_context_adaptive = self._prepare_shared_context(
                        all_peaks_context=all_peaks_context,
                        spectrum_statistics=spectrum_statistics,
                        reference_linewidths=self._reference_linewidths
                    )
                    # Define fit function for grid search evaluation
                    _debug_shown = [False]  # Use list to allow modification in nested function

                    def fit_single_peak(peak, radF1, radF2):
                        """Fit a single peak with given radF1/radF2 for optimization."""
                        # Save original config values BEFORE try block (guarantees they're defined)
                        old_radF1 = config.radF1
                        old_radF2 = config.radF2
                        result = None
                        try:
                            # Temporarily update config for this fit
                            config.radF1 = radF1
                            config.radF2 = radF2

                            # Fit the peak using full PS2D 2D fitting
                            result = self.original_fitter.parent.enhanced_peak_fitting(
                                peak.get('x_ppm'), peak.get('y_ppm'),
                                all_peaks_context=all_peaks_context
                            )

                        except Exception:
                            pass  # Silently handle exceptions during optimization
                        finally:
                            # Always restore config (guaranteed to have old values)
                            config.radF1 = old_radF1
                            config.radF2 = old_radF2

                        return result

                    # Run grid search optimization in PARALLEL mode
                    # Distributes peaks across workers for 5-7× speedup
                    optimal_params = optimizer.optimize(
                        isolated_peaks=good_isolated_peaks,
                        median_lw_f1=spectrum_statistics['lw_f1_median'],
                        median_lw_f2=spectrum_statistics['lw_f2_median'],
                        fit_function=fit_single_peak,  # Fallback for sequential mode
                        progress_callback=progress_callback,
                        use_parallel=True,  # ENABLED: Parallel adaptive optimization
                        shared_context=shared_context_adaptive,  # Spectrum in shared memory
                        all_peaks_context=all_peaks_context
                    )

                    if optimal_params and optimal_params.get('success'):
                        # Apply optimal parameters to config
                        config.radF1 = optimal_params['radF1']
                        config.radF2 = optimal_params['radF2']
                        config.overlap_threshold_x = optimal_params['overlap_threshold_x']
                        config.overlap_threshold_y = optimal_params['overlap_threshold_y']

                        self.optimal_params = optimal_params

                        adaptive_time = time.time() - adaptive_start
                        self.timing_info['adaptive'] = adaptive_time
                        log_info(f"Adaptive optimization completed in {adaptive_time:.1f}s")
                        log_info(f"  → Selected: radF1={optimal_params['radF1']:.4f}, radF2={optimal_params['radF2']:.5f} (mult: {optimal_params['multiplier_f1']:.1f}×, {optimal_params['multiplier_f2']:.1f}×)")

                        # Add optimal overlap thresholds to learned_statistics for main process
                        # (multiprocessing doesn't propagate config changes back to main)
                        if learned_statistics is None:
                            learned_statistics = {}
                        learned_statistics['overlap_threshold_x'] = optimal_params['overlap_threshold_x']
                        learned_statistics['overlap_threshold_y'] = optimal_params['overlap_threshold_y']
                        learned_statistics['radF1'] = optimal_params['radF1']
                        learned_statistics['radF2'] = optimal_params['radF2']

                        # ========== RE-CLUSTER originally clustered peaks with optimal thresholds ==========
                        # CRITICAL: Skip re-clustering when clusters are locked from reference spectrum.
                        # Re-clustering would change cluster membership, breaking volume consistency.
                        if effective_locked_clusters is None and len(original_multi_peak_clusters) > 0:
                            log_info(f"RE-CLUSTERING: Applying optimal overlap thresholds to {sum(len(c) for c in original_multi_peak_clusters)} clustered peaks...")

                            # Extract only originally clustered peak positions
                            originally_clustered_peaks = []
                            for cluster in original_multi_peak_clusters:
                                for pos in cluster:
                                    originally_clustered_peaks.append({
                                        'x_ppm': pos[0],
                                        'y_ppm': pos[1],
                                        'pos_x': pos[0],
                                        'pos_y': pos[1]
                                    })

                            # Re-cluster with optimal thresholds
                            new_clusters = self.original_fitter.parent.identify_overlap_clusters(
                                originally_clustered_peaks,
                                overlap_threshold_x=optimal_params['overlap_threshold_x'],
                                overlap_threshold_y=optimal_params['overlap_threshold_y']
                            )

                            # Separate new isolated and multi-peak
                            new_isolated_from_recluster = [c for c in new_clusters if len(c) == 1]
                            multi_peak_clusters = [c for c in new_clusters if len(c) > 1]

                            log_info(f"Re-clustering result: {len(new_isolated_from_recluster)} now isolated, {len(multi_peak_clusters)} multi-peak")

                            # Add newly isolated peaks to isolated_clusters for PASS 1-bis
                            isolated_clusters = isolated_clusters + new_isolated_from_recluster

                            # CRITICAL FIX: Update clusters to reflect final state after re-clustering
                            # This ensures stored locked_clusters_by_assignment matches what was actually fitted
                            clusters = isolated_clusters + multi_peak_clusters
                            log_info(f"After re-clustering: {len(clusters)} total clusters")
                        elif effective_locked_clusters is not None and len(original_multi_peak_clusters) > 0:
                            log_info(f"RE-CLUSTERING SKIPPED: Using {len(multi_peak_clusters)} locked clusters from reference")
                    else:
                        log_warning("Adaptive optimization failed - using default parameters")
                else:
                    log_warning(f"Insufficient good peaks for adaptive optimization ({len(good_isolated_peaks)} < 10)")
            elif self.enable_adaptive:
                log_warning("Adaptive optimization skipped (insufficient PASS 1 statistics)")

            # ========== PASS 1-bis: Re-fit isolated peaks with optimal parameters ==========
            if optimal_params and optimal_params.get('success') and len(isolated_clusters) > 0:
                log_progress(f"PASS 1-bis: Re-fitting {len(isolated_clusters)} isolated peaks with optimal parameters...")
                pass1bis_start = time.time()

                # Prepare shared context with optimal parameters
                shared_context_pass1bis = self._prepare_shared_context(all_peaks_context, spectrum_statistics=None,
                                                                      reference_linewidths=self._reference_linewidths)

                # Create tasks for all isolated peaks (including newly isolated from re-clustering)
                isolated_tasks_bis = self._create_cluster_tasks(isolated_clusters, peak_metadata, peak_list, shared_context_pass1bis)

                # Execute Pass 1-bis fitting in parallel
                isolated_results = self._execute_parallel_cluster_fitting(isolated_tasks_bis, progress_callback, pass_name="PASS 1-bis")

                # Update statistics with new fits
                spectrum_statistics = self._collect_linewidth_statistics(isolated_results)

                pass1bis_time = time.time() - pass1bis_start
                self.timing_info['pass1bis'] = pass1bis_time
                self.pass1bis_results = isolated_results  # Store PASS1-bis results
                log_info(f"PASS 1-bis completed in {pass1bis_time:.1f}s")

            # ========== PASS 2: Fit remaining clusters using learned/pre-learned statistics ==========
            multi_peak_results = []

            # When using pre-learned statistics, isolated peaks weren't fitted in PASS 1/1-bis
            # So we need to fit them here in PASS 2
            if pre_learned_statistics is not None and len(isolated_results) == 0:
                # Fit ALL clusters (isolated + multi-peak) with pre-learned stats
                all_pass2_clusters = isolated_clusters + multi_peak_clusters
                cluster_type_desc = "all peaks (isolated + multi-peak)"
            else:
                # Normal mode: only multi-peak clusters in PASS 2
                all_pass2_clusters = multi_peak_clusters
                cluster_type_desc = "multi-peak clusters"

            if len(all_pass2_clusters) > 0:
                log_progress(f"PASS 2: Fitting {len(all_pass2_clusters)} {cluster_type_desc} with {'pre-learned' if pre_learned_statistics else 'learned'} statistics...")
                pass2_start = time.time()

                # Prepare shared context WITH spectrum_statistics (Pass 2)
                shared_context_pass2 = self._prepare_shared_context(all_peaks_context, spectrum_statistics=spectrum_statistics,
                                                                    reference_linewidths=self._reference_linewidths)

                # Create tasks for clusters
                pass2_tasks = self._create_cluster_tasks(all_pass2_clusters, peak_metadata, peak_list, shared_context_pass2)

                # Execute Pass 2 fitting in parallel
                multi_peak_results = self._execute_parallel_cluster_fitting(pass2_tasks, progress_callback, pass_name="PASS 2")

                pass2_time = time.time() - pass2_start
                self.timing_info['pass2'] = pass2_time
                self.pass2_results = multi_peak_results  # Store PASS2 results
                log_info(f"PASS 2 completed in {pass2_time:.1f}s")
            else:
                log_warning("No clusters to fit in Pass 2")

            # ========== Consolidate results from both passes ==========
            all_cluster_results = isolated_results + multi_peak_results
            consolidated_results = self._consolidate_cluster_results(all_cluster_results, peak_list, peak_metadata)

            # Store cluster info for ML training
            self.cluster_info = {
                'n_clusters': len(clusters),
                'n_isolated_clusters': len(isolated_clusters),
                'n_multi_peak_clusters': len(multi_peak_clusters),
                'cluster_sizes': [len(c) for c in clusters],
            }

            # ========== Create series_params if adaptive optimization succeeded ==========
            if self.optimal_params and self.optimal_params.get('success'):
                # Build cluster assignments map (position-based, for logging)
                cluster_assignments = {}
                for i, cluster in enumerate(isolated_clusters):
                    for pos in cluster:
                        peak_id = f"peak_{pos[0]:.4f}_{pos[1]:.4f}"
                        cluster_assignments[peak_id] = 'isolated'
                for i, cluster in enumerate(multi_peak_clusters):
                    for pos in cluster:
                        peak_id = f"peak_{pos[0]:.4f}_{pos[1]:.4f}"
                        cluster_assignments[peak_id] = f'cluster_{i}'

                # Build locked_clusters_by_assignment (assignment-based, for cluster locking)
                # Convert position clusters to assignment clusters using peak_metadata
                locked_clusters_by_assignment = []
                for cluster in clusters:
                    cluster_assignments_list = []
                    for pos in cluster:
                        meta = peak_metadata.get(pos)
                        if meta:
                            cluster_assignments_list.append(meta['assignment'])
                    if cluster_assignments_list:
                        locked_clusters_by_assignment.append(cluster_assignments_list)

                self.series_params = create_series_params(
                    optimization_result=self.optimal_params,
                    cluster_assignments=cluster_assignments,
                    reference_spectrum='current'  # Will be updated by GUI with actual filename
                )
                # Store locked clusters for series integration
                self.series_params['locked_clusters_by_assignment'] = locked_clusters_by_assignment
                log_info(f"Series parameters saved with {len(locked_clusters_by_assignment)} locked clusters")

            # Store locked clusters even if adaptive optimization didn't run
            # This ensures series integration always has clusters available
            if not self.series_params:
                # Build locked_clusters_by_assignment from current clusters
                locked_clusters_by_assignment = []
                for cluster in clusters:
                    cluster_assignments_list = []
                    for pos in cluster:
                        meta = peak_metadata.get(pos)
                        if meta:
                            cluster_assignments_list.append(meta['assignment'])
                    if cluster_assignments_list:
                        locked_clusters_by_assignment.append(cluster_assignments_list)

                # Create minimal series_params with just cluster info
                self.series_params = {
                    'locked_clusters_by_assignment': locked_clusters_by_assignment,
                    'success': True
                }
                log_info(f"Stored {len(locked_clusters_by_assignment)} clusters for series integration (no adaptive optimization)")

            elapsed_time = time.time() - start_time
            self.timing_info['total'] = elapsed_time
            log_info(f"Two-pass parallel fitting completed in {elapsed_time:.1f}s")
            log_info(f"Results: {len(consolidated_results)} successful fits from {len(clusters)} clusters")
            if spectrum_statistics:
                log_info("Improvement: Multi-peak clusters fitted with data-driven linewidth priors")

            # Return tuple: (results, learned_statistics)
            # learned_statistics is None if we used pre-learned stats (nothing new was learned)
            return (consolidated_results, learned_statistics)

        except Exception as e:
            log_error(f"Parallel Voigt fitting failed: {e}")
            log_info("Falling back to sequential processing")
            import traceback
            traceback.print_exc()
            results = self._sequential_fallback(peak_list, progress_callback)
            return (results, None)  # Sequential fallback doesn't learn statistics

        finally:
            # Always cleanup shared memory and force resource cleanup
            self._cleanup_shared_memory()
            
    def _cleanup_shared_memory(self):
        """Enhanced cleanup to prevent resource leaks"""
        for memory_block in self.shared_memory_blocks:
            try:
                memory_block.close()
                memory_block.unlink()
            except Exception:
                pass  # Already cleaned or invalid
        self.shared_memory_blocks.clear()

        # Force garbage collection to clear references
        import gc
        gc.collect()

    @staticmethod
    def _normalize_assignment(assignment):
        """Canonicalize a peak assignment so numeric values compare equal regardless
        of how they were stored: 3, 3.0, numpy.int64(3), '3', '3.0' all -> '3'.

        Locked clusters keep the raw peak-list value (int/float/numpy) while the
        current-spectrum context stringifies it, so the two sides must be normalized
        before matching. Non-numeric residue labels (e.g. 'A12N-H') pass through stripped.
        """
        text = str(assignment).strip()
        try:
            value = float(text)
        except (ValueError, TypeError):
            return text
        return str(int(value)) if value.is_integer() else str(value)

    def _convert_locked_clusters_to_positions(self, locked_clusters_by_assignment, all_peaks_context):
        """
        Convert assignment-based locked clusters to position-based clusters.

        Args:
            locked_clusters_by_assignment: List of clusters, each cluster is a list
                of assignment strings, e.g. [['Peak1', 'Peak2'], ['Peak3']]
            all_peaks_context: Current spectrum's peak context with positions

        Returns:
            List of clusters in position format: [[(x1, y1), (x2, y2)], [(x3, y3)], ...]
        """
        # Build assignment → position mapping from current spectrum
        assignment_to_position = {}
        for peak in all_peaks_context:
            assignment = self._normalize_assignment(peak.get('assignment', ''))
            x = peak.get('x_ppm') or peak.get('pos_x')
            y = peak.get('y_ppm') or peak.get('pos_y')
            if assignment and x is not None and y is not None:
                assignment_to_position[assignment] = (x, y)

        # Convert each cluster from assignments to positions
        position_clusters = []
        for cluster_assignments in locked_clusters_by_assignment:
            cluster_positions = []
            for assignment in cluster_assignments:
                key = self._normalize_assignment(assignment)
                if key in assignment_to_position:
                    cluster_positions.append(assignment_to_position[key])
                else:
                    log_warning(f"Locked cluster assignment '{assignment}' not found in current spectrum")

            if cluster_positions:
                position_clusters.append(cluster_positions)

        return position_clusters

    def _collect_linewidth_statistics(self, fitting_results):
        """
        Collect spectrum-wide linewidth statistics from isolated peak fits.

        Used in Pass 1 to learn typical linewidths for Pass 2 cluster fitting.

        Args:
            fitting_results: List of fitting results from isolated peaks

        Returns:
            dict with keys:
                'lw_f1_median': Median total linewidth F1 (ppm)
                'lw_f2_median': Median total linewidth F2 (ppm)
                'lg_f1_median': Median L/G ratio F1
                'lg_f2_median': Median L/G ratio F2
                'n_samples': Number of good fits used for statistics
                'lw_f1_mad': Median absolute deviation F1 (uncertainty)
                'lw_f2_mad': Median absolute deviation F2 (uncertainty)
            or None if insufficient data
        """
        import numpy as np

        linewidth_stats = {
            'lw_f1': [],  # Total linewidth F1 (Lorentz + Gauss)
            'lw_f2': [],  # Total linewidth F2
            'lg_f1': [],  # L/G ratio F1
            'lg_f2': []   # L/G ratio F2
        }

        successful_fits = 0

        # Extract linewidths from successful fits
        # fitting_results is a list of cluster results, each containing peak_results array
        for cluster_result in fitting_results:
            if not cluster_result or not cluster_result.get('success'):
                continue

            # Extract individual peak results from cluster
            for peak_result in cluster_result.get('peak_results', []):
                if not peak_result or not peak_result.get('success'):
                    continue

                # Quality threshold: only use fits with good R²
                r_squared = peak_result.get('r_squared', 0)
                if r_squared < 0.85:
                    continue

                successful_fits += 1

                # Extract linewidth components from x_fit and y_fit sub-dicts
                # (parallel worker stores linewidths in nested structure)
                x_fit = peak_result.get('x_fit', {})
                y_fit = peak_result.get('y_fit', {})

                lw_gau_f2 = x_fit.get('sigma', 0)  # F2 Gaussian (1H)
                lw_lor_f2 = x_fit.get('gamma', 0)  # F2 Lorentzian (1H)
                lw_gau_f1 = y_fit.get('sigma', 0)  # F1 Gaussian (15N)
                lw_lor_f1 = y_fit.get('gamma', 0)  # F1 Lorentzian (15N)

                # Total linewidth = Lorentz + Gauss (both contribute to FWHM)
                lw_f1_total = lw_lor_f1 + lw_gau_f1
                lw_f2_total = lw_lor_f2 + lw_gau_f2

                # ================================================================
                # NUCLEUS-ADAPTIVE SANITY CHECKS
                # ================================================================
                # Use config thresholds instead of hardcoded values
                # This ensures 13C peaks (min 0.025 ppm) are not rejected
                # ================================================================
                from lunaNMR.core.ps2d_config import get_ps2d_config
                config = get_ps2d_config()

                if lw_f1_total > config.min_linewidth_f1:  # Nucleus-adaptive (15N: 0.05, 13C: 0.025)
                    linewidth_stats['lw_f1'].append(lw_f1_total)
                    min_gau_f1 = config.min_linewidth_f1 * 0.2  # 20% of min to avoid division by zero
                    if lw_gau_f1 > min_gau_f1:
                        linewidth_stats['lg_f1'].append(lw_lor_f1 / lw_gau_f1)

                if lw_f2_total > config.min_linewidth_f2:  # Nucleus-adaptive (15N: 0.001, 13C: 0.005)
                    linewidth_stats['lw_f2'].append(lw_f2_total)
                    min_gau_f2 = config.min_linewidth_f2 * 0.2  # 20% of min to avoid division by zero
                    if lw_gau_f2 > min_gau_f2:
                        linewidth_stats['lg_f2'].append(lw_lor_f2 / lw_gau_f2)

        # ================================================================
        # SAMPLE SIZE CHECK AND ADAPTIVE ALPHA
        # ================================================================
        # MAD reliability depends on sample size:
        #   n < 3:  Unreliable, return None (fall back to config)
        #   n < 10: Use conservative alpha=7 (wider bounds)
        #   n >= 10: Use standard alpha=5
        # ================================================================
        MIN_SAMPLES = 3
        RELIABLE_SAMPLES = 10

        n_samples = len(linewidth_stats['lw_f1'])
        if n_samples < MIN_SAMPLES:
            log_warning(f"Insufficient isolated peaks for statistics ({successful_fits} good fits, {n_samples} valid LW)")
            return None

        # Sample-size-dependent multiplier for bounds calculation
        if n_samples < RELIABLE_SAMPLES:
            alpha = 3.0  # Tighter bounds for small samples
            log_warning(f"Small sample size (n={n_samples}), using α=3×MAD")
        else:
            alpha = 5.0  # Standard multiplier for reliable statistics

        # Compute robust statistics using median (resistant to outliers)
        spectrum_lw_f1 = np.median(linewidth_stats['lw_f1'])
        spectrum_lw_f2 = np.median(linewidth_stats['lw_f2'])
        spectrum_lg_f1 = np.median(linewidth_stats['lg_f1']) if len(linewidth_stats['lg_f1']) > 0 else 1.0
        spectrum_lg_f2 = np.median(linewidth_stats['lg_f2']) if len(linewidth_stats['lg_f2']) > 0 else 1.0

        # Calculate MAD (Median Absolute Deviation) for uncertainty estimate
        lw_f1_array = np.array(linewidth_stats['lw_f1'])
        lw_f2_array = np.array(linewidth_stats['lw_f2'])
        mad_f1 = np.median(np.abs(lw_f1_array - spectrum_lw_f1))
        mad_f2 = np.median(np.abs(lw_f2_array - spectrum_lw_f2))

        # Also compute std, min, max for ML training data collection
        std_f1 = np.std(lw_f1_array) if len(lw_f1_array) > 1 else 0.0
        std_f2 = np.std(lw_f2_array) if len(lw_f2_array) > 1 else 0.0

        return {
            'lw_f1_median': spectrum_lw_f1,
            'lw_f2_median': spectrum_lw_f2,
            'lg_f1_median': spectrum_lg_f1,
            'lg_f2_median': spectrum_lg_f2,
            'n_samples': n_samples,  # Use actual valid samples, not successful_fits
            'n_good_peaks': n_samples,  # Alias for ML collector
            'lw_f1_mad': mad_f1,
            'lw_f2_mad': mad_f2,
            'lw_f1_std': std_f1,  # For ML training data
            'lw_f2_std': std_f2,
            'lw_f1_min': float(np.min(lw_f1_array)),
            'lw_f1_max': float(np.max(lw_f1_array)),
            'lw_f2_min': float(np.min(lw_f2_array)),
            'lw_f2_max': float(np.max(lw_f2_array)),
            'alpha': alpha  # Sample-size-dependent multiplier for PASS2 bounds
        }

    def _extract_good_peaks_for_optimization(self, fitting_results):
        """
        Extract good isolated peaks for adaptive optimization.

        Returns list of peak dicts with position, R², and residuals
        for use by AdaptiveOptimizer.

        Args:
            fitting_results: List of cluster results from PASS1 isolated fitting

        Returns:
            List of dicts with x_ppm, y_ppm, r_squared, residuals for good peaks
        """
        good_peaks = []
        for cluster_result in fitting_results:
            if not cluster_result or not cluster_result.get('success'):
                continue
            for peak_result in cluster_result.get('peak_results', []):
                if not peak_result or not peak_result.get('success'):
                    continue

                r_squared = peak_result.get('r_squared', 0)
                if r_squared >= 0.70:
                    # Get position - handle different key naming conventions
                    # Standard format uses peak_x/peak_y or center_x/center_y
                    x_ppm = (peak_result.get('peak_x') or peak_result.get('center_x') or
                             peak_result.get('pos_x') or peak_result.get('pos_f2') or peak_result.get('x_ppm'))
                    y_ppm = (peak_result.get('peak_y') or peak_result.get('center_y') or
                             peak_result.get('pos_y') or peak_result.get('pos_f1') or peak_result.get('y_ppm'))

                    # Skip peaks without valid positions
                    if x_ppm is None or y_ppm is None:
                        continue

                    good_peaks.append({
                        'x_ppm': x_ppm,
                        'y_ppm': y_ppm,
                        'r_squared': r_squared,
                        'residuals': peak_result.get('residuals'),
                        'fwhm_f1': peak_result.get('fwhm_f1'),
                        'fwhm_f2': peak_result.get('fwhm_f2')
                    })

        # Cap at 15 peaks maximum - random sample if more available
        max_peaks = 15
        n_total = len(good_peaks)
        if n_total > max_peaks:
            np.random.seed(42)  # Reproducibility
            indices = np.random.choice(n_total, max_peaks, replace=False)
            good_peaks = [good_peaks[i] for i in sorted(indices)]
            log_info(f"Sampled {max_peaks} peaks from {n_total} for optimization (R² >= 0.70)")
        else:
            log_info(f"Selected {n_total} peaks for optimization (R² >= 0.70)")

        return good_peaks

    def _prepare_shared_context(self, all_peaks_context=None, spectrum_statistics=None,
                                 reference_linewidths=None):
        """
        Create shared memory context containing all data needed
        by worker processes without breaking current logic.

        Args:
            all_peaks_context: List of all peaks for overlap detection
            spectrum_statistics: Optional dict with spectrum-wide linewidth stats from Pass 1
            reference_linewidths: Optional dict mapping assignment -> linewidth params for reuse
        """
        
        # 1. Create shared memory for spectral data
        spectrum_shape = self.original_fitter.nmr_data.shape
        spectrum_dtype = self.original_fitter.nmr_data.dtype
        
        shared_spectrum = shared_memory.SharedMemory(
            create=True, 
            size=self.original_fitter.nmr_data.nbytes
        )
        self.shared_memory_blocks.append(shared_spectrum)
        
        # Copy spectral data to shared memory
        shared_array = np.ndarray(spectrum_shape, dtype=spectrum_dtype, 
                                 buffer=shared_spectrum.buf)
        shared_array[:] = self.original_fitter.nmr_data[:]
        
        # 2. Serialize all essential parameters
        shared_context = {
            # Shared memory identifiers
            'spectrum_memory_name': shared_spectrum.name,
            'spectrum_shape': spectrum_shape,
            'spectrum_dtype': str(spectrum_dtype),
            
            # PPM axes (copy to each worker)  
            'ppm_x_axis': self.original_fitter.ppm_x_axis.copy(),
            'ppm_y_axis': self.original_fitter.ppm_y_axis.copy(),
            
            # Integrator parameters (from parent if available)
            'integrator_params': self._serialize_integrator_params(),
            
            # Preserve all current EnhancedVoigtFitter parameters
            'baseline_params': self._serialize_baseline_params(),
            'fitting_params': self._serialize_fitting_params(),
            'quality_params': self._serialize_quality_params(),
            'gui_params': self._get_gui_params(),
            
            # Advanced parameters
            'fitting_windows': self._get_fitting_windows(),
            'optimization_settings': self._get_optimization_settings(),
            'validation_thresholds': self._get_validation_thresholds(),

            # NEW: Overlap resolution parameters
            'overlap_resolution_params': self._serialize_overlap_resolution_params(),

            # CRITICAL: All peaks context for automatic 2D routing
            'all_peaks_context': all_peaks_context if all_peaks_context is not None else [],

            # NEW: PS2D configuration synchronization
            'ps2d_config': self._serialize_ps2d_config(),

            # NEW: Spectrum-wide linewidth statistics from Pass 1
            'spectrum_statistics': spectrum_statistics,

            # NEW: Per-peak reference linewidths for reuse
            'reference_linewidths': reference_linewidths,

            # Path information for worker imports
            'lunaNMR_path': os.path.dirname(os.path.dirname(__file__)),
        }

        return shared_context

    def _get_gui_params(self):
        """
        Get gui_params from fitter or parent integrator.

        CRITICAL: gui_params must be available for PS2D multi-peak activation.
        Priority: parent integrator > fitter > empty dict
        """
        # Try parent integrator first (most reliable source)
        parent = getattr(self.original_fitter, 'parent', None)
        if parent and hasattr(parent, 'gui_params') and parent.gui_params:
            gui_params = parent.gui_params
            use_ps2d = gui_params.get('use_ps2d_multi_peak', False)
            return gui_params

        # Fall back to fitter's gui_params
        if hasattr(self.original_fitter, 'gui_params') and self.original_fitter.gui_params:
            gui_params = self.original_fitter.gui_params
            use_ps2d = gui_params.get('use_ps2d_multi_peak', False)
            return gui_params

        # Last resort: empty dict
        log_warning("Parallel: No gui_params found, using empty dict (PS2D multi-peak will NOT activate)")
        return {}

    def _serialize_baseline_params(self):
        """Extract all baseline correction parameters"""
        return {
            'lam': getattr(self.original_fitter, 'baseline_lambda', 1e6),
            'p': getattr(self.original_fitter, 'baseline_p', 0.001),
            'max_iter': getattr(self.original_fitter, 'baseline_max_iter', 50),
            'tolerance': getattr(self.original_fitter, 'baseline_tolerance', 1e-6),
            'method': getattr(self.original_fitter, 'baseline_method', 'arpls')
        }

    def _serialize_fitting_params(self):
        """Extract all Voigt fitting parameters"""
        return {
            'max_iter': getattr(self.original_fitter, 'max_iterations', 1000),
            'tolerance': getattr(self.original_fitter, 'convergence_tolerance', 1e-8),
            'method': getattr(self.original_fitter, 'optimization_method', 'leastsq'),
            'bounds_method': getattr(self.original_fitter, 'bounds_method', 'soft'),
            'initial_guess_method': getattr(self.original_fitter, 'initial_guess_method', 'moments')
        }

    def _serialize_quality_params(self):
        """Extract all quality assessment parameters"""
        return {
            'min_r_squared': getattr(self.original_fitter, 'min_r_squared', 0.8),
            'max_residual': getattr(self.original_fitter, 'max_residual', 0.1),
            'snr_threshold': getattr(self.original_fitter, 'snr_threshold', 3.0),
            'quality_checks': getattr(self.original_fitter, 'quality_checks', True)
        }

    def _get_fitting_windows(self):
        """Extract fitting window parameters"""
        return {
            'window_x': getattr(self.original_fitter, 'fitting_window_x', 0.2),
            'window_y': getattr(self.original_fitter, 'fitting_window_y', 2.0),
            'auto_window': getattr(self.original_fitter, 'auto_window_sizing', True)
        }

    def _serialize_integrator_params(self):
        """Extract integrator-specific parameters from parent"""
        params = {}
        
        # Get parent integrator if available
        parent = getattr(self.original_fitter, 'parent', None)
        if parent:
            # Copy important integrator parameters
            for param_name in ['fitting_parameters', 'gui_params', 'processing_mode', 'noise_threshold',
                               'titration_pos_margin_f1', 'titration_pos_margin_f2']:
                if hasattr(parent, param_name):
                    params[param_name] = getattr(parent, param_name)
        
        return params

    def _get_optimization_settings(self):
        """Extract optimization algorithm settings"""
        return {
            'algorithm': getattr(self.original_fitter, 'optimization_algorithm', 'lm'),
            'step_size': getattr(self.original_fitter, 'step_size', 1e-8),
            'gradient_tolerance': getattr(self.original_fitter, 'gradient_tolerance', 1e-12)
        }

    def _serialize_overlap_resolution_params(self):
        """Extract overlap resolution parameters"""
        # Check if overlap_config exists (may be disabled in enhanced_voigt_fitter)
        if not hasattr(self.original_fitter, 'overlap_config'):
            return {
                'enabled': False,
                'threshold': 0.5,
                'config': None
            }

        overlap_config = getattr(self.original_fitter, 'overlap_config', None)
        # Convert OverlapResolutionConfig to dict for multiprocessing serialization
        if overlap_config is not None and hasattr(overlap_config, 'to_dict'):
            overlap_config = overlap_config.to_dict()

        return {
            'enabled': getattr(self.original_fitter, 'overlap_detection_enabled', False),
            'threshold': getattr(self.original_fitter, 'overlap_detection_threshold', 0.5),
            'config': overlap_config
        }

    def _serialize_ps2d_config(self):
        """
        Extract PS2D configuration for workers to ensure consistent behavior.

        Workers need access to PS2D config for:
        - Overlap detection thresholds
        - Position margins (F1, F2)
        - Max iterations for LM optimizer
        - Nucleus type assumptions
        """
        try:
            from lunaNMR.core.ps2d_config import get_ps2d_config
            config = get_ps2d_config()

            ps2d_dict = {
                'nucleus_type': getattr(config, 'nucleus_type', '15N'),
                'overlap_threshold_x': getattr(config, 'overlap_threshold_x', 0.08),
                'overlap_threshold_y': getattr(config, 'overlap_threshold_y', 0.8),
                'max_iterations': getattr(config, 'max_iterations', 100),
                'pos_margin_f1': getattr(config, 'pos_margin_f1', 0.05),
                'pos_margin_f2': getattr(config, 'pos_margin_f2', 0.02),
                'radF1': getattr(config, 'radF1', 0.6),
                'radF2': getattr(config, 'radF2', 0.06)
            }
            return ps2d_dict
        except ImportError:
            log_warning("Parallel: PS2D config not available, using defaults")
            return {}

    def _get_validation_thresholds(self):
        """Extract validation threshold parameters"""
        return {
            'peak_detection_threshold': getattr(self.original_fitter, 'peak_threshold', 0.1),
            'multipeak_threshold': getattr(self.original_fitter, 'multipeak_threshold', 0.3),
            'noise_level': getattr(self.original_fitter, 'noise_level', None),
            'intensity_snr_threshold': 3.0  # Skip fitting if detected_intensity < noise_level * this
        }
        
    def _create_cluster_tasks(self, clusters, peak_metadata, peak_list, shared_context):
        """
        Create cluster-based fitting tasks (IDENTICAL to sequential mode).

        Each task represents ONE cluster that will be fitted ONCE.
        This ensures no peak is fitted multiple times.

        Args:
            clusters: List of clusters from identify_overlap_clusters()
            peak_metadata: Dict mapping (x_ppm, y_ppm) → peak info
            peak_list: Original peak DataFrame
            shared_context: Shared memory context

        Returns:
            List of cluster tasks for parallel workers
        """
        cluster_tasks = []

        # Build all_peaks_context for 2D fitting (need full peak dictionaries)
        all_peaks_context = shared_context['all_peaks_context']

        for cluster_idx, cluster in enumerate(clusters):
            cluster_size = len(cluster)

            # Convert cluster positions to full peak dictionaries
            cluster_dicts = []
            cluster_assignments = []
            cluster_peak_numbers = []

            for peak_x, peak_y in cluster:
                # Find matching peak in all_peaks_context
                peak_dict = None
                for ctx_peak in all_peaks_context:
                    ctx_x = ctx_peak.get('x_ppm') or ctx_peak.get('pos_x')
                    ctx_y = ctx_peak.get('y_ppm') or ctx_peak.get('pos_y')
                    if abs(ctx_x - peak_x) < 0.001 and abs(ctx_y - peak_y) < 0.01:
                        peak_dict = ctx_peak
                        break

                # If not found, create minimal dict
                if peak_dict is None:
                    peak_dict = {
                        'x_ppm': peak_x,
                        'y_ppm': peak_y,
                        'pos_x': peak_x,
                        'pos_y': peak_y,
                        'intensity': None
                    }

                cluster_dicts.append(peak_dict)

                # Get metadata for this peak
                meta = peak_metadata.get((peak_x, peak_y))
                if meta:
                    cluster_assignments.append(meta['assignment'])
                    cluster_peak_numbers.append(meta['peak_number'])
                else:
                    cluster_assignments.append('Unknown')
                    cluster_peak_numbers.append(0)

            # Create cluster task
            cluster_task = {
                'task_id': cluster_idx,
                'cluster_idx': cluster_idx,
                'cluster_size': cluster_size,
                'cluster_positions': cluster,  # List of (x, y) tuples
                'cluster_dicts': cluster_dicts,  # Full peak dictionaries
                'cluster_assignments': cluster_assignments,
                'cluster_peak_numbers': cluster_peak_numbers,
                'shared_context': shared_context
            }
            cluster_tasks.append(cluster_task)

        log_info(f"Created {len(cluster_tasks)} cluster-based tasks")
        log_info(f"  Isolated peaks: {sum(1 for t in cluster_tasks if t['cluster_size'] == 1)}")
        log_info(f"  Overlap groups: {sum(1 for t in cluster_tasks if t['cluster_size'] > 1)}")
        return cluster_tasks
        
    def _execute_parallel_cluster_fitting(self, cluster_tasks, progress_callback, pass_name=None):
        """
        Execute cluster-based fitting in parallel (IDENTICAL logic to sequential mode).

        Each worker fits ONE cluster and returns results for ALL peaks in that cluster.

        Args:
            cluster_tasks: List of cluster task dictionaries
            progress_callback: Optional progress callback function
            pass_name: Optional pass indicator ('PASS 1', 'PASS 2') for progress reporting
        """

        cluster_results = []
        successful_clusters = 0
        failed_clusters = 0
        total_peaks_processed = 0

        pass_prefix = f"{pass_name}: " if pass_name else ""

        try:
            # Test multiprocessing capability first
            with Pool(processes=1) as test_pool:
                test_result = test_pool.apply(_test_parallel_worker, ("test",))
                if test_result != "test_success":
                    raise Exception(f"Multiprocessing test failed: {test_result}")

            # Execute parallel cluster fitting
            with Pool(processes=self.max_workers) as pool:
                # Submit all cluster tasks
                async_results = []
                for task in cluster_tasks:
                    async_result = pool.apply_async(_parallel_cluster_worker, (task,))
                    async_results.append((task['task_id'], task['cluster_size'], async_result))

                # Collect results with progress tracking
                for task_id, cluster_size, async_result in async_results:
                    try:
                        result = async_result.get(timeout=600)  # 10 minute timeout per cluster

                        if result['success']:
                            cluster_results.append(result)
                            successful_clusters += 1
                            total_peaks_processed += result['peaks_fitted']

                            # Progress reporting
                            r_squared = result.get('r_squared', 0)
                            log_info(f"Cluster {task_id + 1}: {result['peaks_fitted']} peaks, R²={r_squared:.3f}")

                            # Log individual peaks for progress dialog
                            if progress_callback:
                                progress = ((successful_clusters + failed_clusters) / len(async_results)) * 100
                                task_desc = f"{pass_prefix}Parallel cluster fitting:\n{successful_clusters + failed_clusters}/{len(async_results)} clusters\n{total_peaks_processed} peaks fitted"
                                for peak_result in result.get('peak_results', []):
                                    assignment = peak_result.get('assignment', 'Unknown')
                                    peak_r2 = peak_result.get('r_squared', peak_result.get('avg_r_squared', 0))
                                    progress_callback(progress, task_desc, f"{assignment} fitted (R²={peak_r2:.3f})")

                        else:
                            failed_clusters += 1
                            error_msg = result.get('error', 'Unknown error')
                            log_warning(f"Cluster {task_id + 1}: {error_msg}")

                            # Log failed cluster for progress dialog
                            if progress_callback:
                                progress = ((successful_clusters + failed_clusters) / len(async_results)) * 100
                                task_desc = f"{pass_prefix}Parallel cluster fitting:\n{successful_clusters + failed_clusters}/{len(async_results)} clusters\n{total_peaks_processed} peaks fitted"
                                progress_callback(progress, task_desc, f"Cluster {task_id + 1} failed: {error_msg}")

                    except mp.TimeoutError:
                        failed_clusters += 1
                        log_warning(f"Cluster {task_id + 1}: Timeout (>10 minutes)")
                    except Exception as e:
                        failed_clusters += 1
                        log_error(f"Cluster {task_id + 1}: Execution error - {str(e)}")

        except Exception as e:
            log_error(f"Parallel cluster execution failed: {e}")
            raise  # Re-raise to trigger fallback

        log_info("Parallel cluster execution summary:")
        log_info(f"  Successful clusters: {successful_clusters}")
        log_info(f"  Failed clusters: {failed_clusters}")
        log_info(f"  Total peaks fitted: {total_peaks_processed}")

        return cluster_results
        
    def _consolidate_cluster_results(self, cluster_results, peak_list, peak_metadata):
        """
        Consolidate cluster results and return peaks in ORIGINAL peak_list order.

        This is IDENTICAL to sequential mode's consolidation logic.

        Args:
            cluster_results: List of cluster fit results from workers
            peak_list: Original peak DataFrame
            peak_metadata: Dict mapping (x, y) → peak info

        Returns:
            List of peak results in original peak_list order
        """

        # Build results cache from cluster results (same as sequential mode)
        # Use peak_number as key to avoid floating point precision issues
        results_by_number = {}  # peak_number → fit_result
        results_by_assignment = {}  # assignment → fit_result

        for cluster_result in cluster_results:
            if not cluster_result['success']:
                continue

            # Each cluster result contains 'peak_results': list of dicts for each peak
            for peak_result in cluster_result.get('peak_results', []):
                peak_number = peak_result.get('peak_number')
                assignment = peak_result.get('assignment')

                if peak_number:
                    results_by_number[peak_number] = peak_result
                if assignment:
                    results_by_assignment[str(assignment)] = peak_result


        # Return results in original peak_list order using peak_number match
        # CRITICAL: Add placeholders for failed fits to maintain 1:1 index mapping (same as sequential mode)
        consolidated_results = []
        failed_count = 0

        for i, (peak_idx, peak_row) in enumerate(peak_list.iterrows()):
            peak_number = i + 1
            assignment = str(peak_row.get('Assignment', f'Peak_{peak_number}'))
            peak_x = float(peak_row['Position_X'])
            peak_y = float(peak_row['Position_Y'])

            # Try matching by peak_number first, then assignment
            result = results_by_number.get(peak_number)
            if not result:
                result = results_by_assignment.get(assignment)

            if result:
                # Ensure processing_mode is set to 'parallel'
                result['processing_mode'] = 'parallel'
                result['peak_number'] = peak_number  # Ensure correct peak_number
                consolidated_results.append(result)
            else:
                # Fitting failed - add placeholder to maintain index alignment (matches sequential mode)
                placeholder = {
                    'assignment': assignment,
                    'peak_number': peak_number,
                    'peak_position': (peak_x, peak_y),
                    'peak_x': peak_x,
                    'peak_y': peak_y,
                    'amplitude': 0.0,
                    'height': 0.0,
                    'volume': 0.0,
                    'r_squared': 0.0,
                    'avg_r_squared': 0.0,
                    'center_x': peak_x,
                    'center_y': peak_y,
                    'sigma_x': 0.0,
                    'gamma_x': 0.0,
                    'sigma_y': 0.0,
                    'gamma_y': 0.0,
                    'fitting_quality': 'Failed',
                    'quality': 'Failed',
                    'success': False,
                    'fitted': False,
                    'method': 'none',
                    'processing_mode': 'parallel',
                    'failure_reason': 'Fitting did not converge or failed acceptance criteria'
                }
                consolidated_results.append(placeholder)
                failed_count += 1

        successful_count = len(consolidated_results) - failed_count
        log_info(f"Consolidated {len(consolidated_results)} peaks in original order (from {len(peak_list)} total peaks)")
        log_info(f"  Successful fits: {successful_count}")
        log_info(f"  Failed fits (placeholders): {failed_count}")

        # Verify 1:1 mapping maintained
        if len(consolidated_results) != len(peak_list):
            log_warning(f"Peak count mismatch! {len(consolidated_results)} != {len(peak_list)}")
            log_warning("This suggests a logic error in consolidation")

        # Warn if success rate is very low
        if successful_count < len(peak_list) * 0.5:
            success_rate = (successful_count / len(peak_list) * 100) if len(peak_list) > 0 else 0
            log_warning(f"Low success rate: {successful_count}/{len(peak_list)} ({success_rate:.1f}%)")
            log_warning("Many peaks failed to fit - check parameters or data quality")

        return consolidated_results
        
    def _sequential_fallback(self, peak_list, progress_callback):
        """
        Fallback to existing sequential processing if parallel fails.
        """
        log_info("Using sequential fallback processing")
        # Call existing enhanced_peak_fitting method for each peak
        results = []
        for i, (peak_idx, peak_row) in enumerate(peak_list.iterrows()):
            peak_x = float(peak_row['Position_X'])
            peak_y = float(peak_row['Position_Y'])
            assignment = peak_row.get('Assignment', f'Peak_{i+1}')

            result = None
            try:
                result = self.original_fitter.enhanced_peak_fitting(peak_x, peak_y, assignment)
                if result:
                    result['processing_mode'] = 'sequential_fallback'
            except Exception as e:
                log_error(f"Sequential fallback failed for peak {i+1}: {e}")
                result = None

            if result:
                results.append(result)
            else:
                # Always append a placeholder so results stay 1:1 with peak_list —
                # a shorter list breaks cascade index alignment downstream.
                results.append({
                    'assignment': assignment, 'peak_number': i + 1,
                    'peak_x': peak_x, 'peak_y': peak_y,
                    'center_x': peak_x, 'center_y': peak_y,
                    'height': 0.0, 'volume': 0.0, 'amplitude': 0.0,
                    'r_squared': 0.0, 'quality': 'Failed', 'fitting_quality': 'Failed',
                    'success': False, 'fitted': False, 'method': 'none',
                    'processing_mode': 'sequential_fallback',
                    'failure_reason': 'Sequential fallback fit failed',
                })

            if progress_callback:
                progress = ((i + 1) / len(peak_list)) * 100
                progress_callback(progress, f"Sequential: {i+1}/{len(peak_list)}", assignment)

        return results


def _test_parallel_worker(_test_data):
    """Simple test function for multiprocessing validation"""
    return "test_success"


def _create_skipped_peak_result(peak_x, peak_y, assignment, peak_number, skip_reason, detected_intensity=0):
    """
    Create a standardized result for peaks that are skipped due to low/negative intensity.

    Args:
        peak_x, peak_y: Peak position
        assignment: Peak assignment string
        peak_number: Peak number in list
        skip_reason: Why the peak was skipped ('negative_intensity' or 'below_snr_threshold')
        detected_intensity: The detected intensity value (for logging)

    Returns:
        Dict with standardized skipped peak result
    """
    return {
        'assignment': assignment,
        'peak_number': peak_number,
        'peak_position': (peak_x, peak_y),
        'peak_x': peak_x,
        'peak_y': peak_y,
        'intensity': 0.0,
        'height': 0.0,
        'volume': 0.0,
        'detected_intensity': detected_intensity,
        'r_squared': 0.0,
        'avg_r_squared': 0.0,
        'center_x': peak_x,
        'center_y': peak_y,
        'sigma_x': 0.0,
        'gamma_x': 0.0,
        'sigma_y': 0.0,
        'gamma_y': 0.0,
        'lw_x': 0.0,
        'lw_y': 0.0,
        'fitting_quality': 'Skipped',
        'quality': 'Skipped',
        'success': False,
        'fitted': False,
        'method': 'skipped',
        'skip_reason': skip_reason,
        'processing_mode': 'parallel',
        'cluster_size': 1
    }


def _should_skip_peak(detected_intensity, noise_level, snr_threshold):
    """
    Check if a peak should be skipped due to low/negative intensity.

    Returns:
        tuple: (should_skip: bool, reason: str or None)
    """
    if detected_intensity is None:
        return False, None  # No intensity info, proceed with fitting

    if detected_intensity <= 0:
        return True, 'negative_intensity'

    if noise_level is not None and noise_level > 0:
        if detected_intensity < noise_level * snr_threshold:
            return True, 'below_snr_threshold'

    return False, None


def _parallel_cluster_worker(cluster_task):
    """
    Worker function that fits ONE cluster (IDENTICAL logic to sequential mode).

    This is the NEW cluster-based worker that ensures deterministic results:
    - Fits entire cluster ONCE (isolated or overlap group)
    - Returns results for ALL peaks in cluster
    - Uses same 2D fitting logic as sequential mode
    - Respects fix_positions and fix_linewidths from GUI

    Args:
        cluster_task: Dict containing cluster info and shared context

    Returns:
        Dict with success flag and peak_results list
    """
    try:
        # Step 1: Initialize worker environment
        worker_integrator = _initialize_worker_fitter(cluster_task['shared_context'])

        cluster_size = cluster_task['cluster_size']
        cluster_dicts = cluster_task['cluster_dicts']
        cluster_assignments = cluster_task['cluster_assignments']
        cluster_peak_numbers = cluster_task['cluster_peak_numbers']

        # Extract fix_positions and fix_linewidths from GUI parameters (CRITICAL!)
        fix_positions = worker_integrator.gui_params.get('fix_positions', False)
        fix_linewidths = worker_integrator.gui_params.get('fix_linewidths', False)

        peak_results = []

        # Extract intensity thresholds for skip check
        validation_thresholds = cluster_task['shared_context'].get('validation_thresholds', {})
        noise_level = validation_thresholds.get('noise_level')
        snr_threshold = validation_thresholds.get('intensity_snr_threshold', 3.0)

        if cluster_size == 1:
            # Isolated peak - standard 1D fitting (IDENTICAL to sequential)
            peak_x, peak_y = cluster_task['cluster_positions'][0]
            assignment = cluster_assignments[0]
            peak_number = cluster_peak_numbers[0]

            # Get detected intensity from original peak data
            detected_intensity = None
            if cluster_dicts:
                detected_intensity = cluster_dicts[0].get('intensity')

            # Pre-check: Skip fitting if intensity is negative or below SNR threshold
            should_skip, skip_reason = _should_skip_peak(detected_intensity, noise_level, snr_threshold)

            if should_skip:
                # Return skipped result - don't waste time fitting noise
                result = _create_skipped_peak_result(
                    peak_x, peak_y, assignment, peak_number, skip_reason, detected_intensity
                )
                peak_results.append(result)
            else:
                # Fit using enhanced_peak_fitting (routes through consensus or standard 1D)
                result = worker_integrator.enhanced_peak_fitting(
                    peak_x, peak_y, assignment,
                    all_peaks_context=cluster_task['shared_context']['all_peaks_context']
                )

                if result:
                    result['peak_number'] = peak_number
                    result['processing_mode'] = 'parallel'
                    result['cluster_size'] = 1

                    # Add detection fields from original peak data (cluster_dicts)
                    if cluster_dicts:
                        orig_peak = cluster_dicts[0]
                        result['detected'] = orig_peak.get('detected', True)
                        result['detection_quality'] = orig_peak.get('detection_quality', 'Matched')
                        result['distance_from_reference'] = orig_peak.get('distance_from_reference', 0)
                        result['distance_from_reference_x'] = orig_peak.get('distance_from_reference_x', 0)
                        result['distance_from_reference_y'] = orig_peak.get('distance_from_reference_y', 0)
                        result['distance_from_reference_elliptical'] = orig_peak.get('distance_from_reference_elliptical', 0)
                        result['reference_retained'] = orig_peak.get('reference_retained', False)
                        # Copy detected_intensity if not already set by enhanced_peak_fitting
                        if result.get('detected_intensity') is None:
                            result['detected_intensity'] = orig_peak.get('intensity')

                    peak_results.append(result)

        else:
            # Overlap group - 2D simultaneous fitting (IDENTICAL to sequential)
            target_assignment = cluster_assignments[0]

            # Build reference linewidths dict for peaks in this cluster
            cluster_reference_linewidths = None
            all_reference_linewidths = cluster_task['shared_context'].get('reference_linewidths')
            if all_reference_linewidths:
                cluster_reference_linewidths = {}
                for assignment in cluster_assignments:
                    if assignment in all_reference_linewidths:
                        cluster_reference_linewidths[assignment] = all_reference_linewidths[assignment]
                # Reference linewidths are used as initial guesses, NOT fixed constraints
                # The optimizer can still adjust linewidths to find the best fit
                # fix_linewidths remains as set by GUI (user choice)

            # Pre-check: If ALL peaks in cluster have negative/low SNR intensity, skip entire cluster
            all_skip = True
            skip_reasons = []
            for idx, peak_dict in enumerate(cluster_dicts):
                detected_intensity = peak_dict.get('intensity')
                should_skip, skip_reason = _should_skip_peak(detected_intensity, noise_level, snr_threshold)
                if not should_skip:
                    all_skip = False
                    break
                skip_reasons.append((idx, skip_reason, detected_intensity))

            if all_skip:
                # All peaks in cluster are noise - skip entire cluster
                cluster_positions = cluster_task['cluster_positions']
                for idx in range(len(cluster_dicts)):
                    peak_x, peak_y = cluster_positions[idx]
                    assignment = cluster_assignments[idx]
                    peak_number = cluster_peak_numbers[idx]
                    detected_intensity = cluster_dicts[idx].get('intensity', 0)
                    # Find the skip reason for this peak
                    skip_reason = 'negative_intensity'  # default
                    for sr_idx, sr_reason, sr_int in skip_reasons:
                        if sr_idx == idx:
                            skip_reason = sr_reason
                            break
                    result = _create_skipped_peak_result(
                        peak_x, peak_y, assignment, peak_number, skip_reason, detected_intensity
                    )
                    result['cluster_size'] = cluster_size
                    result['cluster_idx'] = cluster_task['cluster_idx']
                    peak_results.append(result)

                # Return early - don't waste time on 2D fitting
                return {
                    'success': True,
                    'task_id': cluster_task['task_id'],
                    'cluster_idx': cluster_task['cluster_idx'],
                    'cluster_size': cluster_size,
                    'peaks_fitted': 0,
                    'peaks_skipped': len(peak_results),
                    'r_squared': 0.0,
                    'peak_results': peak_results
                }

            # Call fit_overlap_group_2d() (SAME AS SEQUENTIAL!)
            group_result = worker_integrator.fit_overlap_group_2d(
                cluster_dicts,
                target_assignment,
                peak_assignments=cluster_assignments,
                fix_positions=fix_positions,  # Pass GUI checkbox state
                fix_linewidths=fix_linewidths,  # Pass GUI checkbox state
                reference_linewidths=cluster_reference_linewidths  # Pass reference linewidths for reuse
            )

            if group_result and group_result.get('success', False):
                # Extract 2D region for visualization
                region_2d = worker_integrator.extract_2d_region_for_overlap_group(cluster_dicts)

                if region_2d is None:
                    return {
                        'success': False,
                        'task_id': cluster_task['task_id'],
                        'error': 'Failed to extract 2D region for visualization'
                    }

                # Reconstruct 2D fitted surface
                fitted_2d_surface, individual_surfaces, baseline = worker_integrator._reconstruct_2d_surface(
                    region_2d, group_result['peaks']
                )

                # OPTIMAL MATCHING using Hungarian algorithm (IDENTICAL to sequential)
                # This ensures 1-to-1 matching and prevents peak swapping
                from scipy.optimize import linear_sum_assignment

                cluster_positions = cluster_task['cluster_positions']
                fitted_peaks = group_result['peaks']

                if len(cluster_positions) > 0 and len(fitted_peaks) > 0:
                    # Build cost matrix: distance from each original peak to each fitted peak
                    cost_matrix = np.zeros((len(cluster_positions), len(fitted_peaks)))
                    for i, (peak_x, peak_y) in enumerate(cluster_positions):
                        for j, peak_fit in enumerate(fitted_peaks):
                            fit_x = peak_fit['pos_f2']
                            fit_y = peak_fit['pos_f1']
                            cost_matrix[i, j] = np.sqrt((fit_x - peak_x)**2 + (fit_y - peak_y)**2)

                    # Hungarian algorithm: find optimal 1-to-1 assignment
                    row_indices, col_indices = linear_sum_assignment(cost_matrix)

                    # Create mapping: original peak index → fitted peak index
                    peak_to_fit_mapping = dict(zip(row_indices, col_indices))
                else:
                    peak_to_fit_mapping = {}

                # Extract result for EACH peak using optimal matching
                for i, (peak_x, peak_y) in enumerate(cluster_positions):
                    # Get optimally matched fitted peak
                    if i in peak_to_fit_mapping:
                        best_match = fitted_peaks[peak_to_fit_mapping[i]]
                    else:
                        best_match = None

                    if best_match:
                        # Get detected_intensity from original cluster_dicts (before fitting)
                        orig_detected_intensity = cluster_dicts[i].get('intensity') if i < len(cluster_dicts) else None

                        # Convert to standard format (IDENTICAL to sequential lines 338-387)
                        result = {
                            'assignment': cluster_assignments[i],
                            'peak_number': cluster_peak_numbers[i],
                            'peak_position': (best_match['pos_f2'], best_match['pos_f1']),
                            'peak_x': peak_x,
                            'peak_y': peak_y,
                            'intensity': best_match['intensity'],  # Fitted Voigt parameter (for ML training)
                            'amplitude': best_match.get('amplitude', best_match['intensity']),
                            'height': best_match.get('height', best_match['intensity']),
                            'volume': best_match.get('volume', best_match['intensity']),
                            'detected_intensity': orig_detected_intensity,  # From peak detection (before fitting)
                            'r_squared': group_result['r_squared'],
                            'avg_r_squared': group_result['r_squared'],
                            'center_x': best_match['pos_f2'],
                            'center_y': best_match['pos_f1'],
                            'sigma_x': best_match['lw_gau_f2'],
                            'gamma_x': best_match['lw_lor_f2'],
                            'sigma_y': best_match['lw_gau_f1'],
                            'gamma_y': best_match['lw_lor_f1'],
                            'fitting_quality': 'Excellent' if group_result['r_squared'] > 0.9 else 'Good',
                            'quality': 'Excellent' if group_result['r_squared'] > 0.9 else 'Good',
                            'success': True,
                            'fitted': True,
                            'method': '2d_simultaneous_multi_peak',
                            'processing_mode': 'parallel',
                            'cluster_idx': cluster_task['cluster_idx'],  # Cluster ID for ML training
                            'cluster_size': cluster_size,
                            'overlap_group_size': cluster_size,
                            # Visualization data
                            'x_fit': {
                                'center': best_match['pos_f2'],
                                'sigma': best_match['lw_gau_f2'],
                                'gamma': best_match['lw_lor_f2'],
                                'amplitude': best_match['intensity'],
                                'r_squared': group_result['r_squared'],
                                'success': True,
                                'method': '2d_simultaneous'
                            },
                            'y_fit': {
                                'center': best_match['pos_f1'],
                                'sigma': best_match['lw_gau_f1'],
                                'gamma': best_match['lw_lor_f1'],
                                'amplitude': best_match['intensity'],
                                'r_squared': group_result['r_squared'],
                                'success': True,
                                'method': '2d_simultaneous'
                            },
                            # Cluster-level fit statistics
                            'chi2': group_result.get('chi2', 0),
                            'iterations': group_result.get('iterations', 0),
                            # Convergence flags from PS2D
                            'formal_convergence': group_result.get('formal_convergence', False),
                            'pragmatic_acceptance': group_result.get('pragmatic_acceptance', False),
                            'chi2_reduction_success': group_result.get('chi2_reduction_success', False),
                            # 2D visualization data
                            'region_2d': region_2d,
                            'fitted_2d_surface': fitted_2d_surface,
                            'individual_surfaces': individual_surfaces,
                            'all_peaks': group_result['peaks'],
                            'baseline': baseline  # Baseline offset for visualization
                        }

                        # Add detection fields from original peak data (cluster_dicts)
                        if i < len(cluster_dicts):
                            orig_peak = cluster_dicts[i]
                            result['detected'] = orig_peak.get('detected', True)
                            result['detection_quality'] = orig_peak.get('detection_quality', 'Matched')
                            result['distance_from_reference'] = orig_peak.get('distance_from_reference', 0)
                            result['distance_from_reference_x'] = orig_peak.get('distance_from_reference_x', 0)
                            result['distance_from_reference_y'] = orig_peak.get('distance_from_reference_y', 0)
                            result['distance_from_reference_elliptical'] = orig_peak.get('distance_from_reference_elliptical', 0)
                            result['reference_retained'] = orig_peak.get('reference_retained', False)

                        peak_results.append(result)

        # Return success with all peaks fitted in this cluster
        return {
            'success': True,
            'task_id': cluster_task['task_id'],
            'cluster_idx': cluster_task['cluster_idx'],
            'cluster_size': cluster_size,
            'peaks_fitted': len(peak_results),
            'r_squared': peak_results[0].get('r_squared', 0) if peak_results else 0,
            'peak_results': peak_results
        }

    except Exception as e:
        return {
            'success': False,
            'task_id': cluster_task['task_id'],
            'cluster_idx': cluster_task.get('cluster_idx', '?'),
            'cluster_size': cluster_task.get('cluster_size', 0),
            'peaks_fitted': 0,
            'error': str(e),
            'traceback': traceback.format_exc()
        }


def _initialize_worker_fitter(shared_context):
    """
    Initialize a complete EnhancedVoigtFitter instance in worker process
    that maintains identical behavior to original.
    """
    
    # 1. Setup import path
    lunaNMR_path = shared_context['lunaNMR_path']
    if lunaNMR_path not in sys.path:
        sys.path.insert(0, lunaNMR_path)

    # 2. CRITICAL: Import PS2D modules in worker process (MUST happen before integrator creation)
    try:
        from lunaNMR.core.ps2d_2d_fitter import Ps2dMultiPeakFitter2D
        from lunaNMR.core.ps2d_data_selector import select_data_2d_for_overlap_group
    except ImportError:
        # Workers will fall back to 1D fitting for overlapping peaks
        pass

    # 3. Reconstruct shared spectral data
    shared_spectrum = shared_memory.SharedMemory(
        name=shared_context['spectrum_memory_name']
    )
    
    spectrum_array = np.ndarray(
        shared_context['spectrum_shape'],
        dtype=np.dtype(shared_context['spectrum_dtype']),
        buffer=shared_spectrum.buf
    )
    
    # 3. Create complete integrator instance (needed for real enhanced_peak_fitting)
    from lunaNMR.core.core_integrator import EnhancedVoigtIntegrator
    worker_integrator = EnhancedVoigtIntegrator()
    
    # 4. Restore complete integrator state
    worker_integrator.nmr_data = spectrum_array.copy()
    worker_integrator.ppm_x_axis = shared_context['ppm_x_axis']
    worker_integrator.ppm_y_axis = shared_context['ppm_y_axis']
    
    # 5. Set integrator-specific parameters
    integrator_params = shared_context.get('integrator_params', {})
    for param_name, param_value in integrator_params.items():
        if hasattr(worker_integrator, param_name):
            setattr(worker_integrator, param_name, param_value)

    # CRITICAL: Set gui_params on integrator (needed for PS2D multi-peak activation)
    worker_integrator.gui_params = shared_context['gui_params']

    # NEW: Set spectrum-wide linewidth statistics (Pass 2 only)
    spectrum_statistics = shared_context.get('spectrum_statistics', None)
    if spectrum_statistics:
        worker_integrator.spectrum_statistics = spectrum_statistics

    # CRITICAL: Synchronize PS2D configuration
    ps2d_config = shared_context.get('ps2d_config', {})
    if ps2d_config:
        try:
            from lunaNMR.core.ps2d_config import get_ps2d_config
            config = get_ps2d_config()

            # Directly set attributes on config object (same as GUI does)
            for key, value in ps2d_config.items():
                if hasattr(config, key):
                    setattr(config, key, value)

        except ImportError:
            pass  # Could not synchronize PS2D config

    # DEBUG: Verify gui_params were set
    #if worker_integrator.gui_params:
    #    print(f"✅ Worker: gui_params set on integrator (use_ps2d_multi_peak={worker_integrator.gui_params.get('use_ps2d_multi_peak', 'N/A')})")
    #else:
    #    print(f"❌ Worker: gui_params is None or empty!")

    # 6. Restore enhanced fitter parameters
    if worker_integrator.enhanced_fitter:
        worker_fitter = worker_integrator.enhanced_fitter
        _restore_baseline_params(worker_fitter, shared_context['baseline_params'])
        _restore_fitting_params(worker_fitter, shared_context['fitting_params'])
        _restore_quality_params(worker_fitter, shared_context['quality_params'])
        worker_fitter.gui_params = shared_context['gui_params']

        # CRITICAL: Update enhanced_fitter's parent reference to point to worker_integrator
        # This ensures enhanced_fitter calls adaptive_fit_1d on the correct integrator with gui_params
        worker_fitter.parent = worker_integrator
        
        # Advanced settings
        _restore_fitting_windows(worker_fitter, shared_context['fitting_windows'])
        _restore_optimization_settings(worker_fitter, shared_context['optimization_settings'])
        _restore_validation_thresholds(worker_fitter, shared_context['validation_thresholds'])

        # NEW: Overlap resolution settings
        _restore_overlap_resolution_params(worker_fitter, shared_context['overlap_resolution_params'])

        # Set noise_level on worker_integrator for fit_overlap_group_2d() weak peak detection
        worker_integrator.noise_level = shared_context['validation_thresholds'].get('noise_level')

    # CRITICAL: Verify worker integrator has all required methods for PS2D 2D routing
    required_methods = [
        'check_if_peaks_need_2d_fitting',
        'fit_overlap_group_2d',
        'extract_2d_region_for_overlap_group',
        'identify_overlap_clusters'
    ]
    missing_methods = [m for m in required_methods if not hasattr(worker_integrator, m)]

    # Cleanup shared memory reference in worker
    try:
        shared_spectrum.close()
    except Exception:
        pass

    return worker_integrator


def _restore_baseline_params(worker_fitter, baseline_params):
    """Restore baseline correction parameters"""
    worker_fitter.baseline_lambda = baseline_params['lam']
    worker_fitter.baseline_p = baseline_params['p']
    worker_fitter.baseline_max_iter = baseline_params['max_iter']
    worker_fitter.baseline_tolerance = baseline_params['tolerance']
    worker_fitter.baseline_method = baseline_params['method']

def _restore_fitting_params(worker_fitter, fitting_params):
    """Restore Voigt fitting parameters"""
    worker_fitter.max_iterations = fitting_params['max_iter']
    worker_fitter.convergence_tolerance = fitting_params['tolerance']
    worker_fitter.optimization_method = fitting_params['method']
    worker_fitter.bounds_method = fitting_params['bounds_method']
    worker_fitter.initial_guess_method = fitting_params['initial_guess_method']

def _restore_quality_params(worker_fitter, quality_params):
    """Restore quality assessment parameters"""
    worker_fitter.min_r_squared = quality_params['min_r_squared']
    worker_fitter.max_residual = quality_params['max_residual']
    worker_fitter.snr_threshold = quality_params['snr_threshold']
    worker_fitter.quality_checks = quality_params['quality_checks']

def _restore_fitting_windows(worker_fitter, fitting_windows):
    """Restore fitting window parameters"""
    worker_fitter.fitting_window_x = fitting_windows['window_x']
    worker_fitter.fitting_window_y = fitting_windows['window_y']
    worker_fitter.auto_window_sizing = fitting_windows['auto_window']

def _restore_optimization_settings(worker_fitter, optimization_settings):
    """Restore optimization algorithm settings"""
    worker_fitter.optimization_algorithm = optimization_settings['algorithm']
    worker_fitter.step_size = optimization_settings['step_size']
    worker_fitter.gradient_tolerance = optimization_settings['gradient_tolerance']

def _restore_validation_thresholds(worker_fitter, validation_thresholds):
    """Restore validation threshold parameters"""
    worker_fitter.peak_threshold = validation_thresholds['peak_detection_threshold']
    worker_fitter.multipeak_threshold = validation_thresholds['multipeak_threshold']
    worker_fitter.noise_level = validation_thresholds['noise_level']

def _restore_overlap_resolution_params(worker_fitter, overlap_params):
    """Restore overlap resolution settings"""
    # CRITICAL: Set the overlap_detection_enabled flag
    worker_fitter.overlap_detection_enabled = overlap_params.get('enabled', False)
    worker_fitter.overlap_detection_threshold = overlap_params.get('threshold', 0.5)

    # Restore configuration if present
    config = overlap_params.get('config')
    if config is not None:
        # Convert dict back to OverlapResolutionConfig if needed
        if isinstance(config, dict):
            try:
                from lunaNMR.utils.overlap_config import OverlapResolutionConfig
                worker_fitter.overlap_config = OverlapResolutionConfig(user_config=config)
            except ImportError:
                # Fallback: store as dict
                worker_fitter.overlap_config = config
        else:
            worker_fitter.overlap_config = config

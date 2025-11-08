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
import pandas as pd
from typing import List, Dict, Any, Optional, Tuple
import time
import traceback
import sys
import os
from functools import partial

# CRITICAL FIX for Linux deadlock: Set multiprocessing start method to 'spawn'
# On Linux, default 'fork' causes pipe buffer deadlock when workers print verbose output
# On Mac, 'spawn' is already the default, so this has no effect there
try:
    mp.set_start_method('spawn', force=False)
    print(f"✅ Multiprocessing start method set to 'spawn' (Linux deadlock fix)")
except RuntimeError:
    # Already set (can only be called once)
    current_method = mp.get_start_method()
    print(f"ℹ️  Multiprocessing start method already set to '{current_method}'")

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
    
    def __init__(self, enhanced_voigt_fitter, max_workers=None):
        """
        Initialize parallel Voigt processor
        
        Args:
            enhanced_voigt_fitter: Original EnhancedVoigtFitter instance
            max_workers: Maximum number of worker processes (default: 75% of cores)
        """
        self.original_fitter = enhanced_voigt_fitter
        self.shared_memory_blocks = []  # Track for cleanup
        
        # Configure worker processes
        if max_workers is None:
            self.max_workers = max(1, int(cpu_count() * 0.75))
        else:
            self.max_workers = max_workers
            
        print(f"🚀 ParallelVoigtProcessor initialized with {self.max_workers} workers")
        
    def fit_all_peaks_parallel(self, peak_list, progress_callback=None):
        """
        Complete parallel workflow using IDENTICAL clustering algorithm as sequential mode.

        CRITICAL: This now uses identify_overlap_clusters() ONCE before distributing tasks,
        ensuring identical results to sequential mode (just faster with multiple cores).

        Args:
            peak_list: DataFrame with peak information
            progress_callback: Optional progress callback function

        Returns:
            List of fitting results (same format as sequential)
        """
        if len(peak_list) == 0:
            return []

        print(f"🔬 Starting cluster-based parallel fitting of {len(peak_list)} peaks")
        start_time = time.time()

        try:
            # STEP 1: Build all_peaks_context (identical to sequential mode)
            all_peaks_context = []
            for _, row in peak_list.iterrows():
                intensity = row.get('Height', row.get('Intensity', None))
                all_peaks_context.append({
                    'assignment': str(row.get('Assignment', 'Unknown')),
                    'x_ppm': float(row['Position_X']),
                    'y_ppm': float(row['Position_Y']),
                    'pos_x': float(row['Position_X']),
                    'pos_y': float(row['Position_Y']),
                    'intensity': intensity
                })
            print(f"   🎯 2D overlap detection enabled with {len(all_peaks_context)} peaks context")

            # STEP 2: Call identify_overlap_clusters() ONCE (CRITICAL for deterministic clustering)
            print(f"   🔍 Identifying overlap clusters using hierarchical algorithm...")
            clusters = self.original_fitter.parent.identify_overlap_clusters(all_peaks_context)
            print(f"   ✅ Found {len(clusters)} clusters from {len(peak_list)} peaks")

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

            # Phase 1: Data preparation and sharing
            shared_context = self._prepare_shared_context(all_peaks_context)

            # Phase 2: Create cluster-based tasks (not peak-based!)
            cluster_tasks = self._create_cluster_tasks(clusters, peak_metadata, peak_list, shared_context)

            # Phase 3: Parallel execution (process clusters in parallel)
            cluster_results = self._execute_parallel_cluster_fitting(cluster_tasks, progress_callback)

            # Phase 4: Consolidate and return in original peak order
            consolidated_results = self._consolidate_cluster_results(cluster_results, peak_list, peak_metadata)

            elapsed_time = time.time() - start_time
            print(f"✅ Parallel cluster-based fitting completed in {elapsed_time:.1f}s")
            print(f"   Results: {len(consolidated_results)} successful fits from {len(clusters)} clusters")

            return consolidated_results

        except Exception as e:
            print(f"❌ Parallel Voigt fitting failed: {e}")
            print("🔄 Falling back to sequential processing")
            import traceback
            traceback.print_exc()
            return self._sequential_fallback(peak_list, progress_callback)

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
    
    def _prepare_shared_context(self, all_peaks_context=None):
        """
        Create shared memory context containing all data needed
        by worker processes without breaking current logic.

        Args:
            all_peaks_context: List of all peaks for overlap detection
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
            print(f"📋 Parallel: gui_params from parent integrator (use_ps2d_multi_peak={use_ps2d})")
            return gui_params

        # Fall back to fitter's gui_params
        if hasattr(self.original_fitter, 'gui_params') and self.original_fitter.gui_params:
            gui_params = self.original_fitter.gui_params
            use_ps2d = gui_params.get('use_ps2d_multi_peak', False)
            print(f"📋 Parallel: gui_params from fitter (use_ps2d_multi_peak={use_ps2d})")
            return gui_params

        # Last resort: empty dict
        print(f"⚠️ Parallel: No gui_params found, using empty dict (PS2D multi-peak will NOT activate)")
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
            for param_name in ['fitting_parameters', 'gui_params', 'processing_mode', 'noise_threshold']:
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
            print(f"📋 Parallel: PS2D config serialized:")
            print(f"   nucleus={ps2d_dict['nucleus_type']}, radF1={ps2d_dict['radF1']:.3f}, radF2={ps2d_dict['radF2']:.4f}")
            print(f"   max_iterations={ps2d_dict['max_iterations']}")
            return ps2d_dict
        except ImportError:
            print(f"⚠️ Parallel: PS2D config not available, using defaults")
            return {}

    def _get_validation_thresholds(self):
        """Extract validation threshold parameters"""
        return {
            'peak_detection_threshold': getattr(self.original_fitter, 'peak_threshold', 0.1),
            'multipeak_threshold': getattr(self.original_fitter, 'multipeak_threshold', 0.3),
            'noise_level': getattr(self.original_fitter, 'noise_level', None)
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

        print(f"📋 Created {len(cluster_tasks)} cluster-based tasks")
        print(f"   Isolated peaks: {sum(1 for t in cluster_tasks if t['cluster_size'] == 1)}")
        print(f"   Overlap groups: {sum(1 for t in cluster_tasks if t['cluster_size'] > 1)}")
        return cluster_tasks
        
    def _execute_parallel_cluster_fitting(self, cluster_tasks, progress_callback):
        """
        Execute cluster-based fitting in parallel (IDENTICAL logic to sequential mode).

        Each worker fits ONE cluster and returns results for ALL peaks in that cluster.
        """

        cluster_results = []
        successful_clusters = 0
        failed_clusters = 0
        total_peaks_processed = 0

        print(f"⚡ Starting parallel cluster execution with {self.max_workers} workers")

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
                            print(f"✅ Cluster {task_id + 1}: {result['peaks_fitted']} peaks, R²={r_squared:.3f}")

                        else:
                            failed_clusters += 1
                            error_msg = result.get('error', 'Unknown error')
                            print(f"❌ Cluster {task_id + 1}: {error_msg}")

                            if 'traceback' in result:
                                print(f"   Traceback: {result['traceback'][:200]}...")

                        # Progress callback
                        if progress_callback:
                            progress = ((successful_clusters + failed_clusters) / len(async_results)) * 100
                            progress_callback(
                                progress,
                                f"Parallel cluster fitting: {successful_clusters + failed_clusters}/{len(async_results)} clusters",
                                f"{total_peaks_processed} peaks fitted"
                            )

                    except mp.TimeoutError:
                        failed_clusters += 1
                        print(f"⏰ Cluster {task_id + 1}: Timeout (>10 minutes)")
                    except Exception as e:
                        failed_clusters += 1
                        print(f"❌ Cluster {task_id + 1}: Execution error - {str(e)}")

        except Exception as e:
            print(f"❌ Parallel cluster execution failed: {e}")
            raise  # Re-raise to trigger fallback

        print(f"📊 Parallel cluster execution summary:")
        print(f"   ✅ Successful clusters: {successful_clusters}")
        print(f"   ❌ Failed clusters: {failed_clusters}")
        print(f"   📈 Total peaks fitted: {total_peaks_processed}")

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

        # DEBUG: Print what we have
        print(f"📋 Consolidation: {len(results_by_number)} results by peak_number, {len(results_by_assignment)} by assignment")
        print(f"   Available peak_numbers: {sorted(results_by_number.keys())[:10]}... (showing first 10)")

        # Return results in original peak_list order using peak_number match
        consolidated_results = []
        for i, (peak_idx, peak_row) in enumerate(peak_list.iterrows()):
            peak_number = i + 1
            assignment = str(peak_row.get('Assignment', f'Peak_{peak_number}'))

            # Try matching by peak_number first, then assignment
            result = results_by_number.get(peak_number)
            if not result:
                result = results_by_assignment.get(assignment)

            if result:
                # Ensure processing_mode is set to 'parallel'
                result['processing_mode'] = 'parallel'
                result['peak_number'] = peak_number  # Ensure correct peak_number
                consolidated_results.append(result)

        print(f"📋 Consolidated {len(consolidated_results)} peaks in original order (from {len(peak_list)} total peaks)")

        if len(consolidated_results) < len(peak_list) * 0.5:
            print(f"⚠️ WARNING: Only {len(consolidated_results)}/{len(peak_list)} peaks consolidated")
            print(f"   This suggests a matching problem between cluster results and peak_list")

        return consolidated_results
        
    def _sequential_fallback(self, peak_list, progress_callback):
        """
        Fallback to existing sequential processing if parallel fails.
        """
        print("🔄 Using sequential fallback processing")
        # Call existing enhanced_peak_fitting method for each peak
        results = []
        for i, (peak_idx, peak_row) in enumerate(peak_list.iterrows()):
            peak_x = float(peak_row['Position_X'])
            peak_y = float(peak_row['Position_Y'])
            assignment = peak_row.get('Assignment', f'Peak_{i+1}')
            
            try:
                result = self.original_fitter.enhanced_peak_fitting(peak_x, peak_y, assignment)
                if result:
                    result['processing_mode'] = 'sequential_fallback'
                    results.append(result)
                    
                if progress_callback:
                    progress = ((i + 1) / len(peak_list)) * 100
                    progress_callback(progress, f"Sequential: {i+1}/{len(peak_list)}", assignment)
                    
            except Exception as e:
                print(f"❌ Sequential fallback failed for peak {i+1}: {e}")
                
        return results
        
    def _cleanup_shared_memory(self):
        """Clean up all shared memory blocks"""
        for shm in self.shared_memory_blocks:
            try:
                shm.close()
                shm.unlink()
            except Exception as e:
                print(f"⚠️ Warning: Could not cleanup shared memory: {e}")
        self.shared_memory_blocks.clear()


def _test_parallel_worker(test_data):
    """Simple test function for multiprocessing validation"""
    return "test_success"


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

        if cluster_size == 1:
            # Isolated peak - standard 1D fitting (IDENTICAL to sequential)
            peak_x, peak_y = cluster_task['cluster_positions'][0]
            assignment = cluster_assignments[0]
            peak_number = cluster_peak_numbers[0]

            # Fit using enhanced_peak_fitting (routes through consensus or standard 1D)
            result = worker_integrator.enhanced_peak_fitting(
                peak_x, peak_y, assignment,
                all_peaks_context=cluster_task['shared_context']['all_peaks_context']
            )

            if result:
                result['peak_number'] = peak_number
                result['processing_mode'] = 'parallel'
                result['cluster_size'] = 1
                peak_results.append(result)

        else:
            # Overlap group - 2D simultaneous fitting (IDENTICAL to sequential)
            target_assignment = cluster_assignments[0]

            # Call fit_overlap_group_2d() (SAME AS SEQUENTIAL!)
            group_result = worker_integrator.fit_overlap_group_2d(
                cluster_dicts,
                target_assignment,
                peak_assignments=cluster_assignments,
                fix_positions=fix_positions,  # Pass GUI checkbox state
                fix_linewidths=fix_linewidths  # Pass GUI checkbox state
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
                fitted_2d_surface, individual_surfaces = worker_integrator._reconstruct_2d_surface(
                    region_2d, group_result['peaks']
                )

                # Extract result for EACH peak in cluster (IDENTICAL to sequential)
                for i, (peak_x, peak_y) in enumerate(cluster_task['cluster_positions']):
                    # Find matching peak in group_result
                    best_match = None
                    min_dist = float('inf')
                    for peak_fit in group_result['peaks']:
                        fit_x = peak_fit['pos_f2']
                        fit_y = peak_fit['pos_f1']
                        dist = np.sqrt((fit_x - peak_x)**2 + (fit_y - peak_y)**2)
                        if dist < min_dist:
                            min_dist = dist
                            best_match = peak_fit

                    if best_match:
                        # Convert to standard format (IDENTICAL to sequential lines 338-387)
                        result = {
                            'assignment': cluster_assignments[i],
                            'peak_number': cluster_peak_numbers[i],
                            'peak_position': (best_match['pos_f2'], best_match['pos_f1']),
                            'peak_x': peak_x,
                            'peak_y': peak_y,
                            'amplitude': best_match['intensity'],
                            'height': best_match['intensity'],
                            'volume': best_match['intensity'],
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
                            'region_2d': region_2d,
                            'fitted_2d_surface': fitted_2d_surface,
                            'individual_surfaces': individual_surfaces,
                            'all_peaks': group_result['peaks']
                        }
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
        print(f"✅ Worker: PS2D 2D fitter modules imported successfully")
    except ImportError as e:
        print(f"❌ Worker: PS2D 2D fitter imports failed - {e}")
        print(f"   Workers will fall back to 1D fitting for overlapping peaks")

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

            nucleus = ps2d_config.get('nucleus_type', 'Unknown')
            radF1 = ps2d_config.get('radF1', 'N/A')
            radF2 = ps2d_config.get('radF2', 'N/A')
            max_iter = ps2d_config.get('max_iterations', 'N/A')
            print(f"✅ Worker: PS2D config synchronized (nucleus={nucleus})")
            print(f"   radF1={radF1}, radF2={radF2}, max_iterations={max_iter}")
        except ImportError:
            print(f"⚠️ Worker: Could not synchronize PS2D config (module not available)")

    # DEBUG: Verify gui_params were set
    if worker_integrator.gui_params:
        print(f"✅ Worker: gui_params set on integrator (use_ps2d_multi_peak={worker_integrator.gui_params.get('use_ps2d_multi_peak', 'N/A')})")
    else:
        print(f"❌ Worker: gui_params is None or empty!")

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
        print(f"✅ Worker: Enhanced fitter parent updated to worker_integrator")
        
        # Advanced settings
        _restore_fitting_windows(worker_fitter, shared_context['fitting_windows'])
        _restore_optimization_settings(worker_fitter, shared_context['optimization_settings'])
        _restore_validation_thresholds(worker_fitter, shared_context['validation_thresholds'])

        # NEW: Overlap resolution settings
        _restore_overlap_resolution_params(worker_fitter, shared_context['overlap_resolution_params'])

    # CRITICAL: Verify worker integrator has all required methods for PS2D 2D routing
    required_methods = [
        'check_if_peaks_need_2d_fitting',
        'fit_overlap_group_2d',
        'extract_2d_region_for_overlap_group',
        'identify_overlap_clusters'
    ]
    missing_methods = [m for m in required_methods if not hasattr(worker_integrator, m)]

    if missing_methods:
        print(f"⚠️ Worker missing PS2D methods: {missing_methods}")
        print(f"   2D multi-peak fitting may fail, will fall back to 1D")
    else:
        print(f"✅ Worker: All PS2D 2D routing methods present")

    # Cleanup shared memory reference in worker
    try:
        shared_spectrum.close()
    except:
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
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
        Complete parallel workflow that replicates all current 
        EnhancedVoigtFitter behavior across multiple processes.
        
        Args:
            peak_list: DataFrame with peak information
            progress_callback: Optional progress callback function
            
        Returns:
            List of fitting results (same format as sequential)
        """
        if len(peak_list) == 0:
            return []
            
        print(f"🔬 Starting complete parallel Voigt fitting of {len(peak_list)} peaks")
        start_time = time.time()
        
        try:
            # Phase 1: Data preparation and sharing
            shared_context = self._prepare_shared_context()
            
            # Phase 2: Task distribution  
            peak_tasks = self._create_peak_tasks(peak_list, shared_context)
            
            # Phase 3: Parallel execution
            results = self._execute_parallel_fitting(peak_tasks, progress_callback)
            
            # Phase 4: Result consolidation
            consolidated_results = self._consolidate_results(results)
            
            elapsed_time = time.time() - start_time
            print(f"✅ Parallel Voigt fitting completed in {elapsed_time:.1f}s")
            print(f"   Results: {len(consolidated_results)} successful fits")
            
            return consolidated_results
            
        except Exception as e:
            print(f"❌ Parallel Voigt fitting failed: {e}")
            print("🔄 Falling back to sequential processing")
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
    
    def _prepare_shared_context(self):
        """
        Create shared memory context containing all data needed
        by worker processes without breaking current logic.
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
            'gui_params': getattr(self.original_fitter, 'gui_params', {}),
            
            # Advanced parameters
            'fitting_windows': self._get_fitting_windows(),
            'optimization_settings': self._get_optimization_settings(),
            'validation_thresholds': self._get_validation_thresholds(),
            
            # Path information for worker imports
            'lunaNMR_path': os.path.dirname(os.path.dirname(__file__)),
        }
        
        return shared_context

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

    def _get_validation_thresholds(self):
        """Extract validation threshold parameters"""  
        return {
            'peak_detection_threshold': getattr(self.original_fitter, 'peak_threshold', 0.1),
            'multipeak_threshold': getattr(self.original_fitter, 'multipeak_threshold', 0.3),
            'noise_level': getattr(self.original_fitter, 'noise_level', None)
        }
        
    def _create_peak_tasks(self, peak_list, shared_context):
        """
        Create individual peak fitting tasks that maintain full compatibility
        with current EnhancedVoigtFitter workflow.
        """
        peak_tasks = []
        
        for i, (peak_idx, peak_row) in enumerate(peak_list.iterrows()):
            # Extract all peak information (preserving current logic)
            peak_task = {
                'task_id': i,
                'peak_number': i + 1,
                'assignment': peak_row.get('Assignment', f'Peak_{i+1}'),
                'peak_x': float(peak_row['Position_X']),
                'peak_y': float(peak_row['Position_Y']),
                
                # Include any additional peak-specific parameters
                'height_hint': peak_row.get('Height', None),
                'linewidth_hint': peak_row.get('Linewidth', None), 
                'quality_threshold': peak_row.get('Quality_Threshold', None),
                
                # Peak-specific fitting parameters
                'peak_specific_params': {
                    'custom_window_x': peak_row.get('Window_X', None),
                    'custom_window_y': peak_row.get('Window_Y', None),
                    'fitting_priority': peak_row.get('Priority', 'normal'),
                    'expected_multiplicity': peak_row.get('Multiplicity', 1)
                },
                
                # Reference to shared context
                'shared_context': shared_context
            }
            peak_tasks.append(peak_task)
        
        print(f"📋 Created {len(peak_tasks)} parallel peak fitting tasks")
        return peak_tasks
        
    def _execute_parallel_fitting(self, peak_tasks, progress_callback):
        """
        Execute all peak fitting tasks in parallel while maintaining
        complete compatibility with current error handling and progress tracking.
        """
        
        results = []
        successful_fits = 0
        failed_fits = 0
        
        print(f"⚡ Starting parallel execution with {self.max_workers} workers")
        
        try:
            # Test multiprocessing capability first
            with Pool(processes=1) as test_pool:
                test_result = test_pool.apply(_test_parallel_worker, ("test",))
                if test_result != "test_success":
                    raise Exception(f"Multiprocessing test failed: {test_result}")
            
            # Execute parallel fitting
            with Pool(processes=self.max_workers) as pool:
                # Submit all tasks
                async_results = []
                for task in peak_tasks:
                    async_result = pool.apply_async(_parallel_voigt_worker, (task,))
                    async_results.append((task['task_id'], task['assignment'], async_result))
                
                # Collect results with progress tracking
                for task_id, assignment, async_result in async_results:
                    try:
                        result = async_result.get(timeout=300)  # 5 minute timeout per peak
                        
                        if result['success']:
                            results.append(result)
                            successful_fits += 1
                            
                            # Progress reporting (preserving current format)
                            r_squared = result['result'].get('avg_r_squared', 0)
                            print(f"✅ Peak {result['peak_number']} ({assignment}): R²={r_squared:.3f}")
                            
                        else:
                            failed_fits += 1
                            error_msg = result.get('error', 'Unknown error')
                            print(f"❌ Peak {result.get('peak_number', '?')} ({assignment}): {error_msg}")
                            
                            if 'traceback' in result:
                                print(f"   Traceback: {result['traceback'][:200]}...")
                        
                        # Progress callback (maintaining current interface)
                        if progress_callback:
                            progress = ((successful_fits + failed_fits) / len(async_results)) * 100
                            progress_callback(
                                progress,
                                f"Parallel fitting: {successful_fits + failed_fits}/{len(async_results)} completed",
                                f"Processing {assignment}"
                            )
                    
                    except mp.TimeoutError:
                        failed_fits += 1
                        print(f"⏰ Peak {task_id + 1} ({assignment}): Timeout (>5 minutes)")
                    except Exception as e:
                        failed_fits += 1
                        print(f"❌ Peak {task_id + 1} ({assignment}): Execution error - {str(e)}")
        
        except Exception as e:
            print(f"❌ Parallel pool execution failed: {e}")
            raise  # Re-raise to trigger fallback
        
        print(f"📊 Parallel execution summary:")
        print(f"   ✅ Successful: {successful_fits}")
        print(f"   ❌ Failed: {failed_fits}")
        
        return results
        
    def _consolidate_results(self, parallel_results):
        """
        Consolidate parallel results into format identical to current
        EnhancedVoigtFitter output structure.
        """
        
        # Sort results by task_id to maintain original peak order
        parallel_results.sort(key=lambda x: x['task_id'])
        
        # Extract actual fitting results (preserving all current fields)
        consolidated_results = []
        for result in parallel_results:
            if result['success']:
                # Get exact result from integrator (no modifications)
                fitting_result = result['result']
                
                # Only add minimal parallel metadata if result is valid
                if fitting_result and isinstance(fitting_result, dict):
                    # Create copy to avoid modifying original
                    fitting_result = fitting_result.copy()
                    fitting_result['processing_mode'] = 'parallel'
                    fitting_result['peak_number'] = result['peak_number']
                
                consolidated_results.append(fitting_result)
        
        print(f"📋 Consolidated {len(consolidated_results)} successful parallel results")
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


def _parallel_voigt_worker(peak_task):
    """
    Worker function that replicates COMPLETE EnhancedVoigtFitter
    workflow for a single peak in an isolated process.
    
    This is the core worker that maintains 100% compatibility with
    the existing EnhancedVoigtFitter.enhanced_peak_fitting method.
    """
    try:
        # Step 1: Initialize worker environment
        worker_integrator = _initialize_worker_fitter(peak_task['shared_context'])
        
        # Step 2: Execute REAL enhanced_peak_fitting (same as GUI)
        result = worker_integrator.enhanced_peak_fitting(
            peak_task['peak_x'], 
            peak_task['peak_y'], 
            str(peak_task['assignment'])  # Ensure string type
        )
        
        if result:
            # Add parallel processing metadata
            result['peak_number'] = peak_task['peak_number']
            result['processing_mode'] = 'parallel'
            result['task_id'] = peak_task['task_id']
            
            return {
                'success': True,
                'task_id': peak_task['task_id'],
                'peak_number': peak_task['peak_number'],
                'assignment': peak_task['assignment'],
                'result': result
            }
        else:
            return {
                'success': False,
                'task_id': peak_task['task_id'],
                'peak_number': peak_task.get('peak_number', '?'),
                'assignment': peak_task.get('assignment', 'Unknown'),
                'error': 'Enhanced peak fitting returned None',
            }
        
    except Exception as e:
        return {
            'success': False,
            'task_id': peak_task['task_id'],
            'peak_number': peak_task.get('peak_number', '?'),
            'assignment': peak_task.get('assignment', 'Unknown'),
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
    
    # 2. Reconstruct shared spectral data
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
    
    # 6. Restore enhanced fitter parameters
    if worker_integrator.enhanced_fitter:
        worker_fitter = worker_integrator.enhanced_fitter
        _restore_baseline_params(worker_fitter, shared_context['baseline_params'])
        _restore_fitting_params(worker_fitter, shared_context['fitting_params'])
        _restore_quality_params(worker_fitter, shared_context['quality_params'])
        worker_fitter.gui_params = shared_context['gui_params']
        
        # Advanced settings
        _restore_fitting_windows(worker_fitter, shared_context['fitting_windows'])
        _restore_optimization_settings(worker_fitter, shared_context['optimization_settings'])
        _restore_validation_thresholds(worker_fitter, shared_context['validation_thresholds'])
    
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
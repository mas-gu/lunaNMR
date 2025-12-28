#!/usr/bin/env python3
"""
Parallel Peak Fitting Module

This module provides multiprocessing capabilities for NMR peak fitting to improve
performance on multi-core systems. It handles the entire peak fitting pipeline
with parallel processing.

Key Features:
- Parallel peak fitting using multiprocessing.Pool
- Progress tracking with limited output
- Uses 75% of available CPU cores
- Fault tolerance and error handling

Author: Guillaume Mas
Date: 2025
"""

import multiprocessing as mp
from multiprocessing import Pool, cpu_count
import numpy as np
import time
import traceback


class ParallelPeakFitter:
    """Parallel peak fitting manager"""

    def __init__(self, integrator, max_workers=None):
        """
        Initialize parallel peak fitter

        Args:
            integrator: The VoigtIntegrator instance
            max_workers: Maximum number of worker processes (default: 75% of cores)
        """
        self.integrator = integrator
        if max_workers is None:
            # Use 75% of available cores
            self.max_workers = max(1, int(cpu_count() * 0.75))
        else:
            self.max_workers = max_workers

    def fit_peaks_parallel(self, peak_list, progress_callback=None, all_peaks_context=None):
        """
        Fit multiple peaks in parallel

        Args:
            peak_list: DataFrame with peak information
            progress_callback: Optional callback for progress updates
            all_peaks_context: Optional list of all peaks for 2D overlap detection

        Returns:
            List of fitting results
        """
        if len(peak_list) == 0:
            return []

        # Build all_peaks_context if not provided (for 2D overlap detection)
        if all_peaks_context is None:
            all_peaks_context = []
            for _, row in peak_list.iterrows():
                all_peaks_context.append({
                    'assignment': str(row.get('Assignment', 'Unknown')),
                    'x_ppm': float(row['Position_X']),
                    'y_ppm': float(row['Position_Y']),
                    'pos_x': float(row['Position_X']),
                    'pos_y': float(row['Position_Y'])
                })

        # Serialize integrator data ONCE (not per-peak) for efficiency
        # This avoids redundant serialization of large nmr_data arrays
        shared_integrator_data = self._serialize_integrator_data()

        # Prepare peak data for parallel processing
        peak_tasks = []
        for i, (peak_idx, peak_row) in enumerate(peak_list.iterrows()):
            peak_number = i + 1
            assignment = peak_row.get('Assignment', f'Peak_{peak_number}')
            peak_x = float(peak_row['Position_X'])
            peak_y = float(peak_row['Position_Y'])

            peak_tasks.append({
                'peak_number': peak_number,
                'peak_x': peak_x,
                'peak_y': peak_y,
                'assignment': assignment,
                'integrator_data': shared_integrator_data,
                'all_peaks_context': all_peaks_context
            })

        start_time = time.time()

        # Create multiprocessing pool and fit peaks
        results = []
        successful_fits = 0
        failed_fits = 0

        try:
            # Test if multiprocessing works with a simple task first
            with Pool(processes=self.max_workers) as test_pool:
                test_result = test_pool.apply(_test_worker, ("test",))
                if test_result != "test_ok":
                    raise Exception(f"Multiprocessing test failed: {test_result}")

            with Pool(processes=self.max_workers) as pool:
                # Submit all tasks and collect results
                async_results = []
                for task in peak_tasks:
                    async_result = pool.apply_async(_fit_single_peak_worker, (task,))
                    async_results.append((task['peak_number'], task['assignment'], async_result))

                # Collect results as they complete
                for peak_number, assignment, async_result in async_results:
                    try:
                        result = async_result.get(timeout=120)  # Increased timeout to 120s

                        if result and result.get('success', True):  # Default to True if key missing
                            # Check if it's an actual successful result or error report
                            if 'error' in result:
                                failed_fits += 1
                            else:
                                results.append(result)
                                successful_fits += 1
                        else:
                            failed_fits += 1

                        # Progress callback
                        if progress_callback:
                            progress = (successful_fits + failed_fits) / len(peak_tasks) * 100
                            progress_callback(successful_fits + failed_fits, len(peak_tasks), f"Peak_{peak_number}")

                    except mp.TimeoutError:
                        failed_fits += 1
                    except Exception:
                        failed_fits += 1

        except Exception:
            # Fallback to sequential processing
            return self._fit_peaks_sequential_fallback(peak_list, progress_callback)

        return results

    def _serialize_integrator_data(self):
        """Serialize integrator data for multiprocessing"""
        integrator_data = {
            'nmr_data': self.integrator.nmr_data,
            'ppm_x_axis': self.integrator.ppm_x_axis,
            'ppm_y_axis': self.integrator.ppm_y_axis,
            'gui_params': getattr(self.integrator, 'gui_params', None)
        }

        # Add optional attributes if they exist
        if hasattr(self.integrator, 'nmr_dict'):
            integrator_data['nmr_dict'] = self.integrator.nmr_dict
        if hasattr(self.integrator, 'peak_list'):
            integrator_data['peak_list'] = self.integrator.peak_list
        if hasattr(self.integrator, 'nmr_file_path'):
            integrator_data['nmr_file_path'] = self.integrator.nmr_file_path
        if hasattr(self.integrator, 'peak_list_path'):
            integrator_data['peak_list_path'] = self.integrator.peak_list_path

        # Try to serialize enhanced_fitter (might fail, that's OK)
        try:
            if hasattr(self.integrator, 'enhanced_fitter') and self.integrator.enhanced_fitter is not None:
                # Don't serialize the entire enhanced_fitter (too complex)
                # Instead, let workers create their own
                integrator_data['has_enhanced_fitter'] = True
            else:
                integrator_data['has_enhanced_fitter'] = False
        except Exception:
            integrator_data['has_enhanced_fitter'] = False

        return integrator_data

    def _fit_peaks_sequential_fallback(self, peak_list, progress_callback=None):
        """Fallback to sequential processing if parallel fails"""
        results = []

        for i, (peak_idx, peak_row) in enumerate(peak_list.iterrows()):
            peak_number = i + 1
            assignment = peak_row.get('Assignment', f'Peak_{peak_number}')
            peak_x = float(peak_row['Position_X'])
            peak_y = float(peak_row['Position_Y'])

            try:
                result = self.integrator.enhanced_peak_fitting(peak_x, peak_y, assignment)
                if result:
                    result['peak_number'] = peak_number
                    results.append(result)

                if progress_callback:
                    progress_callback(i + 1, len(peak_list), assignment)

            except Exception:
                pass  # Silent failure - results list tracks success

        return results


def _fit_single_peak_worker(task_data):
    """
    Worker function for parallel peak fitting
    This function runs in a separate process
    """
    try:
        import sys
        import os

        # Add current directory to Python path for imports
        current_dir = os.path.dirname(os.path.abspath(__file__))
        if current_dir not in sys.path:
            sys.path.insert(0, current_dir)

        # Import the integrator class
        from lunaNMR.core.core_integrator import EnhancedVoigtIntegrator

        # Create integrator instance
        integrator = EnhancedVoigtIntegrator()

        # Restore integrator state
        integrator_data = task_data['integrator_data']
        integrator.nmr_data = integrator_data['nmr_data']
        integrator.ppm_x_axis = integrator_data['ppm_x_axis']
        integrator.ppm_y_axis = integrator_data['ppm_y_axis']
        integrator.gui_params = integrator_data['gui_params']

        # Restore other essential attributes
        if 'nmr_dict' in integrator_data:
            integrator.nmr_dict = integrator_data['nmr_dict']
        if 'peak_list' in integrator_data:
            integrator.peak_list = integrator_data['peak_list']
        if 'enhanced_fitter' in integrator_data:
            integrator.enhanced_fitter = integrator_data['enhanced_fitter']

        # Initialize enhanced fitter if not available
        if not hasattr(integrator, 'enhanced_fitter') or integrator.enhanced_fitter is None:
            try:
                from lunaNMR.core.enhanced_voigt_fitter import EnhancedVoigtFitter
                integrator.enhanced_fitter = EnhancedVoigtFitter()
            except ImportError:
                integrator.enhanced_fitter = None

        # Perform the fitting with all_peaks_context for 2D overlap detection
        all_peaks_context = task_data.get('all_peaks_context', None)
        result = integrator.enhanced_peak_fitting(
            task_data['peak_x'],
            task_data['peak_y'],
            task_data['assignment'],
            all_peaks_context=all_peaks_context
        )

        if result:
            result['peak_number'] = task_data['peak_number']
            return result
        else:
            return {
                'success': False,
                'error': 'Fitting returned None',
                'peak_number': task_data['peak_number']
            }

    except ImportError as e:
        return {
            'success': False,
            'error': f'Import error: {str(e)}',
            'traceback': traceback.format_exc(),
            'peak_number': task_data['peak_number']
        }
    except Exception as e:
        return {
            'success': False,
            'error': str(e),
            'traceback': traceback.format_exc(),
            'peak_number': task_data['peak_number']
        }


def _test_worker(_test_data):
    """Simple test function for multiprocessing"""
    return "test_ok"

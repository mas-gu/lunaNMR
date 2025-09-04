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
import pandas as pd
from typing import List, Dict, Any, Tuple
import time
import traceback
from functools import partial

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

        print(f" Parallel fitting initialized with {self.max_workers} workers (75% of {cpu_count()} cores)")

    def fit_peaks_parallel(self, peak_list, progress_callback=None):
        """
        Fit multiple peaks in parallel

        Args:
            peak_list: DataFrame with peak information
            progress_callback: Optional callback for progress updates

        Returns:
            List of fitting results
        """
        if len(peak_list) == 0:
            return []

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
                'integrator_data': self._serialize_integrator_data()
            })

        print(f" Starting parallel fitting of {len(peak_tasks)} peaks...")
        start_time = time.time()

        # Create multiprocessing pool and fit peaks
        results = []
        successful_fits = 0
        failed_fits = 0

        try:
            # Test if multiprocessing works with a simple task first
            print(f" Testing multiprocessing with {self.max_workers} workers...")
            with Pool(processes=self.max_workers) as test_pool:
                test_result = test_pool.apply(_test_worker, ("test",))
                if test_result != "test_ok":
                    raise Exception(f"Multiprocessing test failed: {test_result}")
            print(f" Multiprocessing test successful")

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
                                print(f"❌ Peak {peak_number} ({assignment}): {result['error']}")
                                if 'traceback' in result:
                                    print(f"   Traceback: {result['traceback']}")
                            else:
                                results.append(result)
                                successful_fits += 1
                                r_squared = result.get('avg_r_squared', 0)
                                print(f"✅ Peak {peak_number} ({assignment}): R²={r_squared:.3f}")
                        else:
                            failed_fits += 1
                            error_msg = result.get('error', 'Unknown error') if result else 'No result returned'
                            print(f"❌ Peak {peak_number} ({assignment}): {error_msg}")

                        # Progress callback
                        if progress_callback:
                            progress = (successful_fits + failed_fits) / len(peak_tasks) * 100
                            progress_callback(successful_fits + failed_fits, len(peak_tasks), f"Peak_{peak_number}")

                    except mp.TimeoutError:
                        failed_fits += 1
                        print(f" Peak {peak_number} ({assignment}): Timeout (>120s)")
                    except Exception as e:
                        failed_fits += 1
                        print(f"Peak {peak_number} ({assignment}): Error - {str(e)}")

        except Exception as e:
            print(f"❌ Parallel processing failed: {e}")
            # Fallback to sequential processing
            return self._fit_peaks_sequential_fallback(peak_list, progress_callback)

        elapsed_time = time.time() - start_time
        print(f" Parallel fitting completed in {elapsed_time:.1f}s")
        print(f" Results: {successful_fits} successful, {failed_fits} failed")

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
                    print(f"✅ Peak {peak_number} ({assignment}): Sequential fit successful")
                else:
                    print(f"❌ Peak {peak_number} ({assignment}): Sequential fit failed")

                if progress_callback:
                    progress_callback(i + 1, len(peak_list), assignment)

            except Exception as e:
                print(f" Peak {peak_number} ({assignment}): Sequential error - {str(e)}")

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

        # Perform the fitting
        result = integrator.enhanced_peak_fitting(
            task_data['peak_x'],
            task_data['peak_y'],
            task_data['assignment']
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


def _test_worker(test_data):
    """Simple test function for multiprocessing"""
    return "test_ok"


def _init_worker():
    """Initialize worker processes"""
    # Set up worker process (if needed)
    pass


# Gaussian Mixture Model implementation for very close peaks
class GaussianMixtureModel:
    """Gaussian Mixture Model for overlapping peak detection and fitting"""

    def __init__(self, n_components=2, max_iter=100):
        """
        Initialize GMM for peak fitting

        Args:
            n_components: Number of Gaussian components
            max_iter: Maximum iterations for EM algorithm
        """
        self.n_components = n_components
        self.max_iter = max_iter

    def fit_overlapping_peaks(self, x_data, y_data, n_peaks=None):
        """
        Fit overlapping peaks using Gaussian Mixture Model

        Args:
            x_data: X coordinates
            y_data: Y intensities
            n_peaks: Number of peaks to fit (auto-detect if None)

        Returns:
            Dictionary with fitted parameters
        """
        try:
            from sklearn.mixture import GaussianMixture
            import numpy as np

            # Auto-detect number of peaks if not provided
            if n_peaks is None:
                n_peaks = self._estimate_n_peaks(x_data, y_data)

            # Prepare data for GMM (weight by intensity)
            weights = y_data / np.sum(y_data)
            X = x_data.reshape(-1, 1)

            # Fit Gaussian Mixture Model
            gmm = GaussianMixture(n_components=n_peaks, max_iter=self.max_iter)
            gmm.fit(X, sample_weight=weights)

            # Extract peak parameters
            peaks = []
            for i in range(n_peaks):
                center = gmm.means_[i, 0]
                sigma = np.sqrt(gmm.covariances_[i, 0, 0])
                amplitude = gmm.weights_[i] * np.sum(y_data)

                peaks.append({
                    'center': center,
                    'sigma': sigma,
                    'amplitude': amplitude,
                    'gamma': sigma * 0.5  # Approximate Lorentzian component
                })

            # Calculate combined fit
            y_fitted = self._calculate_combined_voigt(x_data, peaks)
            r_squared = self._calculate_r_squared(y_data, y_fitted)

            return {
                'success': True,
                'n_peaks': n_peaks,
                'peaks': peaks,
                'fitted_curve': y_fitted,
                'r_squared': r_squared,
                'method': 'gaussian_mixture_model'
            }

        except ImportError:
            print("⚠️ sklearn not available, GMM fitting disabled")
            return None
        except Exception as e:
            print(f"GMM fitting failed: {e}")
            return None

    def _estimate_n_peaks(self, x_data, y_data, max_peaks=6):
        """Estimate optimal number of peaks using BIC/AIC"""
        try:
            from sklearn.mixture import GaussianMixture

            best_n_peaks = 1
            best_bic = float('inf')

            weights = y_data / np.sum(y_data)
            X = x_data.reshape(-1, 1)

            for n in range(1, max_peaks + 1):
                try:
                    gmm = GaussianMixture(n_components=n, max_iter=50)
                    gmm.fit(X, sample_weight=weights)
                    bic = gmm.bic(X)

                    if bic < best_bic:
                        best_bic = bic
                        best_n_peaks = n

                except Exception:
                    break  # Stop if fitting fails

            return min(best_n_peaks, 4)  # Limit to reasonable number

        except Exception:
            return 2  # Default fallback

    def _calculate_combined_voigt(self, x_data, peaks):
        """Calculate combined Voigt profile from multiple peaks"""
        from scipy.special import wofz

        result = np.zeros_like(x_data)
        for peak in peaks:
            # Voigt profile calculation
            sigma = peak['sigma']
            gamma = peak['gamma']
            center = peak['center']
            amplitude = peak['amplitude']

            z = ((x_data - center) + 1j * gamma) / (sigma * np.sqrt(2))
            voigt = amplitude * np.real(wofz(z)) / (sigma * np.sqrt(2 * np.pi))
            result += voigt

        return result

    def _calculate_r_squared(self, y_true, y_pred):
        """Calculate R² coefficient"""
        ss_res = np.sum((y_true - y_pred) ** 2)
        ss_tot = np.sum((y_true - np.mean(y_true)) ** 2)
        return 1 - (ss_res / ss_tot) if ss_tot > 0 else 0

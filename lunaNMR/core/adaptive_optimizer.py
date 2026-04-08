# ABOUTME: Adaptive parameter optimization for PS2D fitting via grid search on isolated peaks.
# ABOUTME: Learns optimal radF1/radF2 from reference spectrum and locks params for series integration.
"""
Adaptive PS2D Parameter Optimizer
=================================

This module provides adaptive parameter optimization for PS2D 2D fitting.
Instead of using hardcoded parameters, it learns optimal values from the
actual spectrum via grid search on isolated peaks.

Key Features:
- Grid search over radF1/radF2 multipliers (9 combinations)
- Composite score objective (R² + residual SNR + position + linewidth)
- Train/validation split (70-30) to prevent overfitting
- Linked parameters: radF1 = overlap_threshold_y, radF2 = overlap_threshold_x
- Series parameter locking for consistency across spectra

Future Improvement:
- For larger parameter spaces (>20 combinations), consider replacing grid search
  with Bayesian Optimization (e.g., scikit-optimize) for more efficient exploration.
  Current grid search is optimal for 9 combinations (~1s total).

Usage:
    optimizer = AdaptiveOptimizer(spectrum_data, noise_level)
    optimal_params = optimizer.optimize(isolated_peaks, median_lw_f1, median_lw_f2)
"""

import numpy as np
from typing import Dict, List, Tuple, Optional, Any
from datetime import datetime
import multiprocessing as mp
import time


def _fit_peak_worker(args):
    """
    Worker function for parallel peak evaluation in adaptive optimization.

    Evaluates one peak against all parameter combinations in the grid search.
    Uses the same shared memory pattern as _parallel_cluster_worker.

    Parameters
    ----------
    args : tuple
        (peak_idx, peak, is_train, combinations, shared_context)

        peak_idx : int
            Index of this peak in the sample set
        peak : dict
            Peak data with 'x_ppm', 'y_ppm'
        is_train : bool
            True if training peak, False if validation peak
        combinations : list of tuple
            List of (mult_f1, mult_f2) combinations to evaluate
        shared_context : dict
            Shared context from ParallelVoigtProcessor (spectrum in shared memory)

    Returns
    -------
    tuple : (peak_idx, is_train, scores_by_combo)
        scores_by_combo is list of (combo_idx, score, success) for each combination
    """
    peak_idx, peak, is_train, combinations, shared_context = args

    median_lw_f1 = shared_context['adaptive_context']['median_lw_f1']
    median_lw_f2 = shared_context['adaptive_context']['median_lw_f2']
    noise_level = shared_context['adaptive_context']['noise_level']
    weights = shared_context['adaptive_context']['weights']
    all_peaks_context = shared_context['adaptive_context']['all_peaks_context']

    # Initialize worker fitter (reconstructs integrator with shared memory spectrum)
    from lunaNMR.core.parallel_voigt_processor import _initialize_worker_fitter
    worker_integrator = _initialize_worker_fitter(shared_context)

    # Import PS2D config
    from lunaNMR.core.ps2d_config import get_ps2d_config
    config = get_ps2d_config()

    scores_by_combo = []

    # Evaluate this peak against all combinations
    for combo_idx, (mult_f1, mult_f2) in enumerate(combinations):
        radF1 = mult_f1 * median_lw_f1
        radF2 = mult_f2 * median_lw_f2

        # Save original config values
        old_radF1 = config.radF1
        old_radF2 = config.radF2

        try:
            # Temporarily update config for this combination
            config.radF1 = radF1
            config.radF2 = radF2

            # Fit peak using full PS2D pipeline (same as sequential mode)
            result = worker_integrator.enhanced_peak_fitting(
                peak.get('x_ppm'),
                peak.get('y_ppm'),
                all_peaks_context=all_peaks_context
            )

            if result is None:
                scores_by_combo.append((combo_idx, 0.0, False))
                continue

            # Extract r_squared
            r_squared = _extract_r_squared_helper(result)
            if r_squared is None:
                scores_by_combo.append((combo_idx, 0.0, False))
                continue

            # Compute composite score
            score = _compute_composite_score_helper(
                result,
                initial_pos=(peak.get('x_ppm'), peak.get('y_ppm')),
                noise_level=noise_level,
                config=config,
                weights=weights
            )
            scores_by_combo.append((combo_idx, score, True))

        except Exception:
            scores_by_combo.append((combo_idx, 0.0, False))

        finally:
            # Always restore config
            config.radF1 = old_radF1
            config.radF2 = old_radF2

    return (peak_idx, is_train, scores_by_combo)


def _extract_r_squared_helper(result: Dict) -> Optional[float]:
    """Helper function for extracting R² in worker processes."""
    if result is None:
        return None

    if 'r_squared' in result:
        return result['r_squared']

    if 'fit_result' in result and isinstance(result['fit_result'], dict):
        if 'r_squared' in result['fit_result']:
            return result['fit_result']['r_squared']

    if 'peak_results' in result and isinstance(result['peak_results'], list):
        for pr in result['peak_results']:
            if pr and 'r_squared' in pr:
                return pr['r_squared']

    if isinstance(result, list) and len(result) > 0:
        if isinstance(result[0], dict) and 'r_squared' in result[0]:
            return result[0]['r_squared']

    if 'results' in result and isinstance(result['results'], dict):
        if 'r_squared' in result['results']:
            return result['results']['r_squared']

    return None


def _compute_composite_score_helper(result: Dict,
                                    initial_pos: Tuple[float, float],
                                    noise_level: float,
                                    config,
                                    weights: Dict[str, float]) -> float:
    """Helper function for computing composite score in worker processes."""
    # Extract weights
    weight_r2 = weights['r2']
    weight_snr = weights['snr']
    weight_position = weights['position']
    weight_linewidth = weights['linewidth']

    # 1. R² component
    r2 = _extract_r_squared_helper(result) or 0
    r2_score = max(0, min(1, r2))

    # 2. Residual SNR component
    residuals = result.get('residuals', None)
    if residuals is not None and len(residuals) > 0 and noise_level > 0:
        res_snr = np.std(residuals) / noise_level
        snr_score = 1.0 / (1.0 + abs(res_snr - 1.0))
    else:
        snr_score = 0.5

    # 3. Position stability component
    fitted_x = result.get('pos_x') or result.get('pos_f2', initial_pos[0])
    fitted_y = result.get('pos_y') or result.get('pos_f1', initial_pos[1])
    pos_shift = np.sqrt((fitted_x - initial_pos[0])**2 +
                       (fitted_y - initial_pos[1])**2)
    pos_score = 1.0 / (1.0 + 100 * pos_shift)

    # 4. Linewidth physicality component
    lw_f1 = result.get('fwhm_f1') or result.get('lw_y', None)
    if lw_f1 is not None and config is not None:
        min_lw = getattr(config, 'min_linewidth_f1', 0.01)
        max_lw = getattr(config, 'max_linewidth_f1', 1.0)
        if min_lw < lw_f1 < max_lw * 0.9:
            lw_score = 1.0
        else:
            lw_score = 0.5
    else:
        lw_score = 0.75

    # Weighted combination
    composite = (weight_r2 * r2_score +
                weight_snr * snr_score +
                weight_position * pos_score +
                weight_linewidth * lw_score)

    return composite


class AdaptiveOptimizer:
    """
    Adaptive parameter optimizer for PS2D fitting.

    Finds optimal radF1/radF2 by grid search on isolated peaks using
    a composite score that prevents R²-only pathologies.
    """

    # Grid search multipliers
    MULTIPLIERS = [1.0, 1.5, 2.0]  # scaled back from [1.0, 1.5, 2.0, 2.5, 3.0]

    # Train/validation split ratio
    TRAIN_RATIO = 0.7

    # Minimum peaks required for optimization
    MIN_PEAKS_FOR_OPTIMIZATION = 10

    # Composite score weights
    WEIGHT_R2 = 0.2  # was 0.4
    WEIGHT_SNR = 0.6  # was 0.3
    WEIGHT_POSITION = 0.1  # was 0.2
    WEIGHT_LINEWIDTH = 0.1  # was 0.1

    def __init__(self, noise_level: float, config=None):
        """
        Initialize the adaptive optimizer.

        Parameters
        ----------
        noise_level : float
            Estimated noise level of the spectrum (for composite score)
        config : PS2DConfig, optional
            PS2D configuration for fallback values and constraints
        """
        self.noise_level = noise_level
        self.config = config

        # Results storage
        self.optimization_results = None
        self.best_params = None
        self.search_history = []

    def optimize(self,
                 isolated_peaks: List[Dict],
                 median_lw_f1: float,
                 median_lw_f2: float,
                 fit_function: callable,
                 progress_callback: callable = None,
                 use_parallel: bool = False,
                 shared_context: Dict = None,
                 all_peaks_context: List[Dict] = None) -> Dict[str, Any]:
        """
        Run grid search optimization to find optimal radF1/radF2.

        Re-fits isolated peaks with different radF1/radF2 combinations
        and selects the parameters that give the best composite score.

        Parameters
        ----------
        isolated_peaks : list of dict
            Isolated peaks with fit results from PASS1.
            Each peak should have: 'x_ppm', 'y_ppm', 'r_squared', etc.
        median_lw_f1 : float
            Median linewidth in F1 dimension (from PASS1)
        median_lw_f2 : float
            Median linewidth in F2 dimension (from PASS1)
        fit_function : callable
            Function to fit a single peak: fit_function(peak, radF1, radF2) -> result_dict
            Result must contain 'r_squared' key.
        progress_callback : callable, optional
            Function to report progress: progress_callback(percent, task_desc, detail)
        use_parallel : bool, optional
            Enable parallel optimization (default: False)
        shared_context : dict, optional
            Shared context for parallel mode (spectrum in shared memory, etc.)
        all_peaks_context : list of dict, optional
            List of all peaks for overlap detection (required for parallel mode)

        Returns
        -------
        dict : Optimal parameters and optimization metadata

        Note: Future improvement could use Bayesian Optimization for larger
        parameter spaces (>20 combinations). Current grid search is appropriate
        for 9 combinations.
        """
        # Check minimum peaks requirement
        if len(isolated_peaks) < self.MIN_PEAKS_FOR_OPTIMIZATION:
            return self._get_fallback_params(median_lw_f1, median_lw_f2, reason="insufficient_peaks")

        # Route to parallel or sequential mode
        if use_parallel and shared_context is not None:
            return self._optimize_parallel(
                isolated_peaks, median_lw_f1, median_lw_f2,
                shared_context, all_peaks_context, progress_callback
            )

        # Sequential mode (original implementation)
        return self._optimize_sequential(
            isolated_peaks, median_lw_f1, median_lw_f2,
            fit_function, progress_callback
        )

    def _optimize_sequential(self,
                            isolated_peaks: List[Dict],
                            median_lw_f1: float,
                            median_lw_f2: float,
                            fit_function: callable,
                            progress_callback: callable = None) -> Dict[str, Any]:
        """
        Sequential grid search optimization (original implementation).

        This is the original optimize() logic extracted into a separate method.
        """
        # Split into train/validation sets
        train_peaks, val_peaks = self._split_train_validation(isolated_peaks)

        opt_start_time = time.time()

        n_combinations = len(self.MULTIPLIERS) ** 2
        n_fits_per_combo = len(train_peaks) + len(val_peaks)

        # Grid search
        best_score = -np.inf
        best_params = None
        self.search_history = []

        for combo_idx, (mult_f1, mult_f2) in enumerate([(m1, m2) for m1 in self.MULTIPLIERS for m2 in self.MULTIPLIERS]):
            combo_start = time.time()
            radF1 = mult_f1 * median_lw_f1
            radF2 = mult_f2 * median_lw_f2

            # Evaluate on training set
            train_scores = []
            train_failures = 0
            total_fits_in_combo = len(train_peaks) + len(val_peaks)
            for peak_idx, peak in enumerate(train_peaks):
                # Update GUI progress
                if progress_callback:
                    # Calculate overall progress across all combinations
                    fits_done = combo_idx * n_fits_per_combo + peak_idx
                    total_fits = n_combinations * n_fits_per_combo
                    progress = int(100 * fits_done / total_fits)
                    task_desc = f"OPTIMIZE [{combo_idx+1}/{n_combinations}]: {len(train_scores)}/{peak_idx+1} peaks"
                    detail = f"F1={mult_f1:.1f}×, F2={mult_f2:.1f}× (train {peak_idx+1}/{len(train_peaks)})"
                    progress_callback(progress, task_desc, detail)

                try:
                    result = fit_function(peak, radF1, radF2)
                    if result is None:
                        train_failures += 1
                        continue

                    # Extract r_squared - handle nested result structures
                    r_squared = self._extract_r_squared(result)
                    if r_squared is not None:
                        score = self.compute_composite_score(
                            result,
                            initial_pos=(peak.get('x_ppm'), peak.get('y_ppm'))
                        )
                        train_scores.append(score)
                    else:
                        train_failures += 1
                except Exception:
                    train_failures += 1
                    continue

            # Check if we got enough training scores
            if len(train_scores) == 0:
                continue

            train_score = np.mean(train_scores)

            # Evaluate on validation set
            val_scores = []
            val_failures = 0
            for peak_idx, peak in enumerate(val_peaks):
                # Update GUI progress
                if progress_callback:
                    fits_done = combo_idx * n_fits_per_combo + len(train_peaks) + peak_idx
                    total_fits = n_combinations * n_fits_per_combo
                    progress = int(100 * fits_done / total_fits)
                    task_desc = f"OPTIMIZE [{combo_idx+1}/{n_combinations}]: validating"
                    detail = f"F1={mult_f1:.1f}×, F2={mult_f2:.1f}× (val {peak_idx+1}/{len(val_peaks)})"
                    progress_callback(progress, task_desc, detail)

                try:
                    result = fit_function(peak, radF1, radF2)
                    if result is None:
                        val_failures += 1
                        continue

                    r_squared = self._extract_r_squared(result)
                    if r_squared is not None:
                        score = self.compute_composite_score(
                            result,
                            initial_pos=(peak.get('x_ppm'), peak.get('y_ppm'))
                        )
                        val_scores.append(score)
                    else:
                        val_failures += 1
                except Exception:
                    val_failures += 1
                    continue

            if len(val_scores) == 0:
                continue

            val_score = np.mean(val_scores)
            combo_time = time.time() - combo_start

            # Record in history
            self.search_history.append({
                'multiplier_f1': mult_f1,
                'multiplier_f2': mult_f2,
                'radF1': radF1,
                'radF2': radF2,
                'train_score': train_score,
                'val_score': val_score,
                'n_train_fits': len(train_scores),
                'n_val_fits': len(val_scores)
            })

            # Select based on VALIDATION score (not training)
            if val_score > best_score:
                best_score = val_score
                best_params = {
                    'multiplier_f1': mult_f1,
                    'multiplier_f2': mult_f2,
                    'radF1': radF1,
                    'radF2': radF2,
                    'train_score': train_score,
                    'val_score': val_score
                }

        # Check if optimization succeeded
        if best_params is None or best_score < 0.6:
            return self._get_fallback_params(median_lw_f1, median_lw_f2, reason="poor_score")

        # Build result
        self.best_params = best_params
        self.optimization_results = {
            'radF1': best_params['radF1'],
            'radF2': best_params['radF2'],
            'overlap_threshold_x': best_params['radF2'],  # Linked
            'overlap_threshold_y': best_params['radF1'],  # Linked
            'multiplier_f1': best_params['multiplier_f1'],
            'multiplier_f2': best_params['multiplier_f2'],
            'validation_score': best_params['val_score'],
            'train_score': best_params['train_score'],
            'generalization_gap': best_params['train_score'] - best_params['val_score'],
            'n_peaks_used': len(isolated_peaks),
            'n_train': len(train_peaks),
            'n_val': len(val_peaks),
            'median_lw_f1': median_lw_f1,
            'median_lw_f2': median_lw_f2,
            'search_history': self.search_history,
            'success': True
        }

        return self.optimization_results

    def _optimize_parallel(self,
                          isolated_peaks: List[Dict],
                          median_lw_f1: float,
                          median_lw_f2: float,
                          shared_context: Dict,
                          all_peaks_context: List[Dict],
                          progress_callback: callable = None) -> Dict[str, Any]:
        """
        Parallel grid search optimization using multiprocessing.

        Distributes peaks across workers. Each worker evaluates all 9 combinations
        for its assigned peak(s), then the main process aggregates results.

        Parameters
        ----------
        isolated_peaks : list of dict
            Isolated peaks for optimization
        median_lw_f1 : float
            Median F1 linewidth
        median_lw_f2 : float
            Median F2 linewidth
        shared_context : dict
            Shared context with spectrum in shared memory
        all_peaks_context : list of dict
            All peaks for overlap detection
        progress_callback : callable, optional
            Progress reporting function

        Returns
        -------
        dict : Optimal parameters and optimization metadata
        """
        # Split into train/validation sets
        train_peaks, val_peaks = self._split_train_validation(isolated_peaks)

        opt_start_time = time.time()

        # Build list of all parameter combinations
        combinations = [(m1, m2) for m1 in self.MULTIPLIERS for m2 in self.MULTIPLIERS]
        n_combinations = len(combinations)

        # Build extended shared context for adaptive optimization
        adaptive_context = {
            'median_lw_f1': median_lw_f1,
            'median_lw_f2': median_lw_f2,
            'noise_level': self.noise_level,
            'weights': {
                'r2': self.WEIGHT_R2,
                'snr': self.WEIGHT_SNR,
                'position': self.WEIGHT_POSITION,
                'linewidth': self.WEIGHT_LINEWIDTH
            },
            'all_peaks_context': all_peaks_context
        }

        # Add adaptive context to shared_context
        extended_context = shared_context.copy()
        extended_context['adaptive_context'] = adaptive_context

        # Prepare worker tasks (one task per peak)
        worker_tasks = []

        # Training peaks
        for peak_idx, peak in enumerate(train_peaks):
            worker_tasks.append((peak_idx, peak, True, combinations, extended_context))

        # Validation peaks
        for peak_idx, peak in enumerate(val_peaks):
            # Use offset index to distinguish from train peaks
            worker_tasks.append((len(train_peaks) + peak_idx, peak, False, combinations, extended_context))

        # Execute parallel optimization
        from multiprocessing import Pool
        max_workers = mp.cpu_count() - 2  # Leave 2 cores free
        max_workers = max(1, max_workers)

        # Storage for aggregated results
        # combo_results[combo_idx] = {'train_scores': [], 'val_scores': []}
        combo_results = [{'train_scores': [], 'val_scores': []} for _ in range(n_combinations)]

        try:
            with Pool(processes=max_workers) as pool:
                # Submit all peak tasks
                async_results = []
                for task in worker_tasks:
                    async_result = pool.apply_async(_fit_peak_worker, (task,))
                    async_results.append(async_result)

                # Collect results with progress tracking
                completed = 0
                for async_result in async_results:
                    try:
                        peak_idx, is_train, scores_by_combo = async_result.get(timeout=600)

                        # Aggregate scores by combination
                        for combo_idx, score, success in scores_by_combo:
                            if success:
                                if is_train:
                                    combo_results[combo_idx]['train_scores'].append(score)
                                else:
                                    combo_results[combo_idx]['val_scores'].append(score)

                        # Progress reporting
                        completed += 1
                        if progress_callback:
                            progress = int(100 * completed / len(async_results))
                            task_desc = f"OPTIMIZE [parallel]: {completed}/{len(async_results)} peaks"
                            detail = f"Completed peak {peak_idx + 1}"
                            progress_callback(progress, task_desc, detail)

                    except mp.TimeoutError:
                        # Skip this peak if it times out
                        completed += 1
                        continue
                    except Exception:
                        # Skip this peak on error
                        completed += 1
                        continue

        except Exception as e:
            # Fallback to sequential mode on error
            from lunaNMR.utils.output_manager import log_warning
            log_warning(f"Parallel optimization failed: {e}. Falling back to sequential.")
            return self._get_fallback_params(median_lw_f1, median_lw_f2, reason="parallel_failed")

        # Aggregate results across combinations
        best_score = -np.inf
        best_params = None
        self.search_history = []

        for combo_idx, (mult_f1, mult_f2) in enumerate(combinations):
            train_scores = combo_results[combo_idx]['train_scores']
            val_scores = combo_results[combo_idx]['val_scores']

            if len(train_scores) == 0 or len(val_scores) == 0:
                continue

            train_score = np.mean(train_scores)
            val_score = np.mean(val_scores)

            radF1 = mult_f1 * median_lw_f1
            radF2 = mult_f2 * median_lw_f2

            # Record in history
            self.search_history.append({
                'multiplier_f1': mult_f1,
                'multiplier_f2': mult_f2,
                'radF1': radF1,
                'radF2': radF2,
                'train_score': train_score,
                'val_score': val_score,
                'n_train_fits': len(train_scores),
                'n_val_fits': len(val_scores)
            })

            # Select based on VALIDATION score
            if val_score > best_score:
                best_score = val_score
                best_params = {
                    'multiplier_f1': mult_f1,
                    'multiplier_f2': mult_f2,
                    'radF1': radF1,
                    'radF2': radF2,
                    'train_score': train_score,
                    'val_score': val_score
                }

        # Check if optimization succeeded
        if best_params is None or best_score < 0.6:
            return self._get_fallback_params(median_lw_f1, median_lw_f2, reason="poor_score")

        # Build result
        self.best_params = best_params
        self.optimization_results = {
            'radF1': best_params['radF1'],
            'radF2': best_params['radF2'],
            'overlap_threshold_x': best_params['radF2'],  # Linked
            'overlap_threshold_y': best_params['radF1'],  # Linked
            'multiplier_f1': best_params['multiplier_f1'],
            'multiplier_f2': best_params['multiplier_f2'],
            'validation_score': best_params['val_score'],
            'train_score': best_params['train_score'],
            'generalization_gap': best_params['train_score'] - best_params['val_score'],
            'n_peaks_used': len(isolated_peaks),
            'n_train': len(train_peaks),
            'n_val': len(val_peaks),
            'median_lw_f1': median_lw_f1,
            'median_lw_f2': median_lw_f2,
            'search_history': self.search_history,
            'success': True,
            'parallel': True  # Flag to indicate parallel mode was used
        }

        return self.optimization_results

    def _extract_r_squared(self, result: Dict) -> Optional[float]:
        """
        Extract R² value from various result structures.

        enhanced_peak_fitting can return results in different formats:
        - Direct: {'r_squared': 0.95, ...}
        - Nested: {'success': True, 'fit_result': {'r_squared': 0.95, ...}}
        - List: [{'r_squared': 0.95, ...}]  (for multi-peak)

        Parameters
        ----------
        result : dict
            Fit result from enhanced_peak_fitting

        Returns
        -------
        float or None : R² value if found, None otherwise
        """
        if result is None:
            return None

        # Direct access
        if 'r_squared' in result:
            return result['r_squared']

        # Check nested 'fit_result'
        if 'fit_result' in result and isinstance(result['fit_result'], dict):
            if 'r_squared' in result['fit_result']:
                return result['fit_result']['r_squared']

        # Check 'peak_results' (for cluster results)
        if 'peak_results' in result and isinstance(result['peak_results'], list):
            for pr in result['peak_results']:
                if pr and 'r_squared' in pr:
                    return pr['r_squared']

        # Check if result is a list (multi-peak result)
        if isinstance(result, list) and len(result) > 0:
            if isinstance(result[0], dict) and 'r_squared' in result[0]:
                return result[0]['r_squared']

        # Check 'results' key
        if 'results' in result and isinstance(result['results'], dict):
            if 'r_squared' in result['results']:
                return result['results']['r_squared']

        return None

    def compute_composite_score(self,
                                result: Dict,
                                initial_pos: Tuple[float, float]) -> float:
        """
        Compute composite fit score for parameter optimization.

        Combines multiple metrics to prevent R²-only pathologies:
        - R² (40%): goodness of fit
        - Residual SNR (30%): should be ~1.0 (residuals = noise)
        - Position stability (20%): small shift from initial
        - Linewidth physicality (10%): not hitting bounds

        Parameters
        ----------
        result : dict
            Fit result with keys: 'r_squared', 'residuals', 'pos_x'/'pos_f2',
            'pos_y'/'pos_f1', 'fwhm_f1', etc.
        initial_pos : tuple
            Initial (x, y) position in ppm

        Returns
        -------
        float : Composite score in [0, 1]
        """
        # 1. R² component (40%)
        r2 = self._extract_r_squared(result) or 0
        r2_score = max(0, min(1, r2))

        # 2. Residual SNR component (30%)
        # Target: residual_snr ≈ 1.0 (residuals should look like pure noise)
        residuals = result.get('residuals', None)
        if residuals is not None and len(residuals) > 0 and self.noise_level > 0:
            res_snr = np.std(residuals) / self.noise_level
            snr_score = 1.0 / (1.0 + abs(res_snr - 1.0))
        else:
            snr_score = 0.5  # Neutral if can't compute

        # 3. Position stability component (20%)
        # Penalize large shifts from initial position
        fitted_x = result.get('pos_x') or result.get('pos_f2', initial_pos[0])
        fitted_y = result.get('pos_y') or result.get('pos_f1', initial_pos[1])
        pos_shift = np.sqrt((fitted_x - initial_pos[0])**2 +
                           (fitted_y - initial_pos[1])**2)
        pos_score = 1.0 / (1.0 + 100 * pos_shift)

        # 4. Linewidth physicality component (10%)
        # Penalize if linewidth hits bounds
        lw_f1 = result.get('fwhm_f1') or result.get('lw_y', None)
        if lw_f1 is not None and self.config is not None:
            min_lw = getattr(self.config, 'min_linewidth_f1', 0.01)
            max_lw = getattr(self.config, 'max_linewidth_f1', 1.0)
            if min_lw < lw_f1 < max_lw * 0.9:
                lw_score = 1.0
            else:
                lw_score = 0.5  # Penalize boundary hits
        else:
            lw_score = 0.75  # Neutral if can't check

        # Weighted combination
        composite = (self.WEIGHT_R2 * r2_score +
                    self.WEIGHT_SNR * snr_score +
                    self.WEIGHT_POSITION * pos_score +
                    self.WEIGHT_LINEWIDTH * lw_score)

        return composite

    def _split_train_validation(self, peaks: List[Dict]) -> Tuple[List[Dict], List[Dict]]:
        """
        Split peaks into training and validation sets (70-30).

        Uses random shuffle with fixed seed for reproducibility.
        """
        n_peaks = len(peaks)
        n_train = int(n_peaks * self.TRAIN_RATIO)

        # Shuffle with fixed seed for reproducibility
        indices = np.arange(n_peaks)
        np.random.seed(42)
        np.random.shuffle(indices)

        train_indices = indices[:n_train]
        val_indices = indices[n_train:]

        train_peaks = [peaks[i] for i in train_indices]
        val_peaks = [peaks[i] for i in val_indices]

        return train_peaks, val_peaks

    def _get_fallback_params(self,
                             median_lw_f1: float,
                             median_lw_f2: float,
                             reason: str) -> Dict[str, Any]:
        """
        Return fallback parameters when optimization fails.

        Uses 1.5× multiplier as default (middle of search space).
        """
        default_mult = 1.5

        return {
            'radF1': default_mult * median_lw_f1,
            'radF2': default_mult * median_lw_f2,
            'overlap_threshold_x': default_mult * median_lw_f2,
            'overlap_threshold_y': default_mult * median_lw_f1,
            'multiplier_f1': default_mult,
            'multiplier_f2': default_mult,
            'validation_score': None,
            'train_score': None,
            'generalization_gap': None,
            'n_peaks_used': 0,
            'median_lw_f1': median_lw_f1,
            'median_lw_f2': median_lw_f2,
            'search_history': [],
            'success': False,
            'fallback_reason': reason
        }

def estimate_noise_level(spectrum: np.ndarray, method: str = 'corner') -> float:
    """
    Estimate noise level from 2D spectrum.

    Parameters
    ----------
    spectrum : np.ndarray
        2D spectrum data
    method : str
        'corner': Use corner regions (typically baseline)
        'mad': Use MAD-based robust estimate

    Returns
    -------
    float : Estimated noise level (standard deviation)
    """
    if method == 'corner':
        ny, nx = spectrum.shape
        # Use 10% corners
        size_y = max(1, ny // 10)
        size_x = max(1, nx // 10)

        corners = [
            spectrum[:size_y, :size_x],           # Top-left
            spectrum[:size_y, -size_x:],          # Top-right
            spectrum[-size_y:, :size_x],          # Bottom-left
            spectrum[-size_y:, -size_x:]          # Bottom-right
        ]

        # Use median of corner stds (robust to one corner having signal)
        corner_stds = [np.std(c) for c in corners]
        return np.median(corner_stds)

    elif method == 'mad':
        # Median Absolute Deviation (robust to outliers/signals)
        # MAD * 1.4826 ≈ std for normal distribution
        median_val = np.median(spectrum)
        mad = np.median(np.abs(spectrum - median_val))
        return 1.4826 * mad

    else:
        raise ValueError(f"Unknown noise estimation method: {method}")


def create_series_params(optimization_result: Dict[str, Any],
                         cluster_assignments: Dict[str, Any],
                         reference_spectrum: str = None) -> Dict[str, Any]:
    """
    Create series parameters structure from optimization results.

    Parameters
    ----------
    optimization_result : dict
        Result from AdaptiveOptimizer.optimize()
    cluster_assignments : dict
        Peak ID to cluster assignment mapping
    reference_spectrum : str, optional
        Name/path of reference spectrum

    Returns
    -------
    dict : Complete series parameters for locking
    """
    return {
        # Optimized parameters
        'radF1': optimization_result['radF1'],
        'radF2': optimization_result['radF2'],
        'overlap_threshold_x': optimization_result['overlap_threshold_x'],
        'overlap_threshold_y': optimization_result['overlap_threshold_y'],

        # Multipliers (for documentation)
        'multiplier_f1': optimization_result['multiplier_f1'],
        'multiplier_f2': optimization_result['multiplier_f2'],

        # Cluster assignments (locked for series)
        'cluster_assignments': cluster_assignments,

        # Reference statistics (for deviation detection)
        'reference_LW_f1': optimization_result['median_lw_f1'],
        'reference_LW_f2': optimization_result['median_lw_f2'],

        # Metadata
        'reference_spectrum': reference_spectrum,
        'optimization_score': optimization_result['validation_score'],
        'success': optimization_result['success'],
        'timestamp': datetime.now().isoformat()
    }


def check_lw_deviation(current_lw_f1: float,
                       current_lw_f2: float,
                       series_params: Dict[str, Any],
                       threshold: float = 0.3) -> Tuple[bool, str]:
    """
    Check if current spectrum's linewidths deviate from reference.

    Parameters
    ----------
    current_lw_f1 : float
        Current spectrum's median F1 linewidth
    current_lw_f2 : float
        Current spectrum's median F2 linewidth
    series_params : dict
        Locked series parameters with reference_LW_f1/f2
    threshold : float
        Maximum allowed relative deviation (default: 0.3 = 30%)

    Returns
    -------
    tuple : (has_deviation, warning_message)
    """
    ref_lw_f1 = series_params.get('reference_LW_f1', current_lw_f1)
    ref_lw_f2 = series_params.get('reference_LW_f2', current_lw_f2)

    # Compute relative deviations
    dev_f1 = abs(current_lw_f1 - ref_lw_f1) / ref_lw_f1 if ref_lw_f1 > 0 else 0
    dev_f2 = abs(current_lw_f2 - ref_lw_f2) / ref_lw_f2 if ref_lw_f2 > 0 else 0

    warnings = []
    has_deviation = False

    if dev_f1 > threshold:
        has_deviation = True
        direction = "broader" if current_lw_f1 > ref_lw_f1 else "narrower"
        warnings.append(f"F1 linewidth {direction} by {dev_f1*100:.0f}% "
                       f"({current_lw_f1:.4f} vs {ref_lw_f1:.4f} ppm)")

    if dev_f2 > threshold:
        has_deviation = True
        direction = "broader" if current_lw_f2 > ref_lw_f2 else "narrower"
        warnings.append(f"F2 linewidth {direction} by {dev_f2*100:.0f}% "
                       f"({current_lw_f2:.5f} vs {ref_lw_f2:.5f} ppm)")

    if has_deviation:
        warning_msg = "Linewidth deviation detected:\n" + "\n".join(f"  - {w}" for w in warnings)
        warning_msg += "\nLocked parameters may not be optimal for this spectrum."
    else:
        warning_msg = ""

    return has_deviation, warning_msg

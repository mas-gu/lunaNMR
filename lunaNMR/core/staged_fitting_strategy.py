#!/usr/bin/env python3
"""
Staged Fitting Strategy Module for LunaNMR
===========================================

Implements a 3-stage parameter release strategy inspired by  multi-stage
optimization approach. This progressive fitting strategy improves convergence and
robustness for complex multi-peak Voigt profile fitting.

The strategy uses warm starts where each stage uses the previous stage's result
as the initial guess (p0), progressively releasing constraints:

Stage 1: Lock positions, fit amplitudes & widths
Stage 2: Release positions with tight bounds (±5-10%)
Stage 3: Full optimization with wide bounds

This approach prevents the optimizer from getting stuck in local minima by
first optimizing the easier parameters (amplitudes and widths) before allowing
the more challenging position parameters to vary.

Author: Guillaume Mas
Date: 2025
"""

import numpy as np
from scipy.optimize import curve_fit, OptimizeWarning
from scipy.special import wofz
from typing import Dict, List, Tuple, Callable, Optional, Any
import warnings
import time
from dataclasses import dataclass, field


@dataclass
class StageResult:
    """Data class for individual stage fitting results"""
    stage: int
    parameters: Optional[np.ndarray] = None
    covariance: Optional[np.ndarray] = None
    r_squared: float = 0.0
    reduced_chi_squared: float = np.inf
    residual_std: float = np.inf
    success: bool = False
    message: str = ""
    iterations: int = 0
    function_evals: int = 0
    fit_time: float = 0.0
    convergence_quality: str = "failed"

    def to_dict(self) -> Dict:
        """Convert to dictionary for logging"""
        return {
            'stage': self.stage,
            'success': self.success,
            'r_squared': self.r_squared,
            'reduced_chi_squared': self.reduced_chi_squared,
            'residual_std': self.residual_std,
            'fit_time': self.fit_time,
            'convergence_quality': self.convergence_quality,
            'message': self.message
        }


@dataclass
class StagedFittingResult:
    """Data class for complete staged fitting results"""
    final_parameters: Optional[np.ndarray] = None
    final_covariance: Optional[np.ndarray] = None
    final_r_squared: float = 0.0
    final_reduced_chi_squared: float = np.inf
    success: bool = False
    stage_history: List[StageResult] = field(default_factory=list)
    total_time: float = 0.0
    best_stage: int = 0
    convergence_trajectory: Dict[str, List] = field(default_factory=dict)
    quality_assessment: Dict = field(default_factory=dict)

    def to_dict(self) -> Dict:
        """Convert to dictionary for logging"""
        return {
            'success': self.success,
            'final_r_squared': self.final_r_squared,
            'final_reduced_chi_squared': self.final_reduced_chi_squared,
            'total_time': self.total_time,
            'best_stage': self.best_stage,
            'stage_history': [s.to_dict() for s in self.stage_history],
            'quality_assessment': self.quality_assessment
        }


class BoundsGenerator:
    """
    Generate parameter bounds for staged fitting.

    Handles bounds generation for multi-peak Voigt models with parameter ordering:
    [amp1, center1, sigma1, gamma1, amp2, center2, sigma2, gamma2, ..., baseline]
    """

    # Default parameter ranges for NMR Voigt fitting
    DEFAULT_RANGES = {
        '1H': {'min': 5.5, 'max': 12.0, 'typical_width': 0.01},
        '15N': {'min': 100.0, 'max': 140.0, 'typical_width': 0.5},
        '13C': {'min': 5.0, 'max': 50.0, 'typical_width': 0.6},
        'default': {'min': 0.0, 'max': 200.0, 'typical_width': 1.0}
    }

    @staticmethod
    def fixed_positions(peak_centers: List[float],
                       n_peaks: int,
                       amplitude_range: Tuple[float, float] = (0, np.inf),
                       width_ranges: Dict[str, Tuple[float, float]] = None,
                       baseline_range: Tuple[float, float] = (-np.inf, np.inf),
                       nucleus_type: str = '1H') -> Tuple[List, List]:
        """
        Stage 1 bounds: positions exactly locked at initial values.

        Parameters:
        -----------
        peak_centers : List[float]
            Initial peak center positions
        n_peaks : int
            Number of peaks
        amplitude_range : Tuple[float, float]
            (min, max) bounds for amplitudes
        width_ranges : Dict[str, Tuple[float, float]]
            Bounds for 'sigma' and 'gamma' widths
        baseline_range : Tuple[float, float]
            (min, max) bounds for baseline
        nucleus_type : str
            Nucleus type for width estimation

        Returns:
        --------
        Tuple[List, List] : (lower_bounds, upper_bounds)
        """
        if width_ranges is None:
            typical_width = BoundsGenerator.DEFAULT_RANGES.get(
                nucleus_type, BoundsGenerator.DEFAULT_RANGES['default']
            )['typical_width']
            width_ranges = {
                'sigma': (typical_width * 0.1, typical_width * 10),
                'gamma': (typical_width * 0.1, typical_width * 10)
            }

        lower_bounds = []
        upper_bounds = []

        for i in range(n_peaks):
            # Amplitude bounds
            lower_bounds.append(amplitude_range[0])
            upper_bounds.append(amplitude_range[1])

            # Position bounds - LOCKED at initial value
            center = peak_centers[i]
            epsilon = abs(center) * 1e-10 if center != 0 else 1e-10
            lower_bounds.append(center - epsilon)
            upper_bounds.append(center + epsilon)

            # Sigma bounds
            lower_bounds.append(width_ranges['sigma'][0])
            upper_bounds.append(width_ranges['sigma'][1])

            # Gamma bounds
            lower_bounds.append(width_ranges['gamma'][0])
            upper_bounds.append(width_ranges['gamma'][1])

        # Baseline bounds
        lower_bounds.append(baseline_range[0])
        upper_bounds.append(baseline_range[1])

        return lower_bounds, upper_bounds

    @staticmethod
    def constrained_positions(peak_centers: List[float],
                            tolerance: float,
                            n_peaks: int,
                            amplitude_range: Tuple[float, float] = (0, np.inf),
                            width_ranges: Dict[str, Tuple[float, float]] = None,
                            baseline_range: Tuple[float, float] = (-np.inf, np.inf),
                            nucleus_type: str = '1H',
                            x_range: Optional[Tuple[float, float]] = None) -> Tuple[List, List]:
        """
        Stage 2 bounds: positions allowed within ±tolerance of initial values.

        Parameters:
        -----------
        peak_centers : List[float]
            Initial peak center positions
        tolerance : float
            Relative tolerance for position variation (e.g., 0.05 = ±5%)
        n_peaks : int
            Number of peaks
        amplitude_range : Tuple[float, float]
            (min, max) bounds for amplitudes
        width_ranges : Dict[str, Tuple[float, float]]
            Bounds for 'sigma' and 'gamma' widths
        baseline_range : Tuple[float, float]
            (min, max) bounds for baseline
        nucleus_type : str
            Nucleus type for width estimation
        x_range : Optional[Tuple[float, float]]
            Data range to enforce bounds within

        Returns:
        --------
        Tuple[List, List] : (lower_bounds, upper_bounds)
        """
        if width_ranges is None:
            typical_width = BoundsGenerator.DEFAULT_RANGES.get(
                nucleus_type, BoundsGenerator.DEFAULT_RANGES['default']
            )['typical_width']
            width_ranges = {
                'sigma': (typical_width * 0.1, typical_width * 10),
                'gamma': (typical_width * 0.1, typical_width * 10)
            }

        lower_bounds = []
        upper_bounds = []

        for i in range(n_peaks):
            # Amplitude bounds
            lower_bounds.append(amplitude_range[0])
            upper_bounds.append(amplitude_range[1])

            # Position bounds - constrained based on window size
            center = peak_centers[i]

            # Calculate window-based position tolerance (not center-based)
            if x_range is not None:
                window_width = x_range[1] - x_range[0]

                # Infer dimension and apply appropriate tolerance
                # X dimension (1H, narrow): want ±0.03 ppm for ~0.15 ppm window → 0.2 factor
                # Y dimension (15N/13C, wide): want ±0.4 ppm for ~4.0 ppm window → 0.1 factor
                if window_width < 1.0:
                    position_delta_fraction = 0.2  # X dimension: tighter absolute constraint
                else:
                    position_delta_fraction = 0.1  # Y dimension: moderate absolute constraint

                position_delta = window_width * position_delta_fraction
            else:
                # Fallback if no x_range provided (should not happen)
                position_delta = abs(center) * tolerance

            # Ensure minimum movement is possible
            if position_delta < 1e-6:
                position_delta = 1e-6

            pos_lower = center - position_delta
            pos_upper = center + position_delta

            # Enforce data range if provided
            if x_range is not None:
                pos_lower = max(pos_lower, x_range[0])
                pos_upper = min(pos_upper, x_range[1])

            # Safety check
            if pos_lower >= pos_upper:
                pos_lower = center - 1e-6
                pos_upper = center + 1e-6

            lower_bounds.append(pos_lower)
            upper_bounds.append(pos_upper)

            # Sigma bounds
            lower_bounds.append(width_ranges['sigma'][0])
            upper_bounds.append(width_ranges['sigma'][1])

            # Gamma bounds
            lower_bounds.append(width_ranges['gamma'][0])
            upper_bounds.append(width_ranges['gamma'][1])

        # Baseline bounds
        lower_bounds.append(baseline_range[0])
        upper_bounds.append(baseline_range[1])

        return lower_bounds, upper_bounds

    @staticmethod
    def full_freedom(peak_centers: List[float],
                    x_range: Tuple[float, float],
                    n_peaks: int,
                    amplitude_range: Tuple[float, float] = (0, np.inf),
                    width_ranges: Dict[str, Tuple[float, float]] = None,
                    baseline_range: Tuple[float, float] = (-np.inf, np.inf),
                    nucleus_type: str = '1H',
                    min_separation: Optional[float] = None) -> Tuple[List, List]:
        """
        Stage 3 bounds: wide bounds allowing full exploration within data range.

        Parameters:
        -----------
        peak_centers : List[float]
            Initial peak center positions (for reference)
        x_range : Tuple[float, float]
            Data range (min, max) for position bounds
        n_peaks : int
            Number of peaks
        amplitude_range : Tuple[float, float]
            (min, max) bounds for amplitudes
        width_ranges : Dict[str, Tuple[float, float]]
            Bounds for 'sigma' and 'gamma' widths
        baseline_range : Tuple[float, float]
            (min, max) bounds for baseline
        nucleus_type : str
            Nucleus type for width estimation
        min_separation : Optional[float]
            Minimum required separation between peaks

        Returns:
        --------
        Tuple[List, List] : (lower_bounds, upper_bounds)
        """
        if width_ranges is None:
            typical_width = BoundsGenerator.DEFAULT_RANGES.get(
                nucleus_type, BoundsGenerator.DEFAULT_RANGES['default']
            )['typical_width']
            width_ranges = {
                'sigma': (typical_width * 0.05, typical_width * 20),
                'gamma': (typical_width * 0.05, typical_width * 20)
            }

        lower_bounds = []
        upper_bounds = []

        # Calculate window-based position tolerance (same as Stage 2)
        x_min, x_max = x_range
        window_width = x_max - x_min

        # Infer dimension and apply appropriate tolerance
        # X dimension (1H, narrow): want ±0.03 ppm for ~0.15 ppm window → 0.2 factor
        # Y dimension (15N/13C, wide): want ±0.4 ppm for ~4.0 ppm window → 0.1 factor
        if window_width < 1.0:
            position_delta_fraction = 0.2  # X dimension: tighter absolute constraint
        else:
            position_delta_fraction = 0.1  # Y dimension: moderate absolute constraint

        position_tolerance = window_width * position_delta_fraction

        for i in range(n_peaks):
            # Amplitude bounds
            lower_bounds.append(amplitude_range[0])
            upper_bounds.append(amplitude_range[1])

            # Position bounds - constrained to ±tolerance around peak center
            center = peak_centers[i]
            pos_lower = max(x_min, center - position_tolerance)
            pos_upper = min(x_max, center + position_tolerance)

            lower_bounds.append(pos_lower)
            upper_bounds.append(pos_upper)

            # Sigma bounds - wider for full exploration
            lower_bounds.append(width_ranges['sigma'][0])
            upper_bounds.append(width_ranges['sigma'][1])

            # Gamma bounds - wider for full exploration
            lower_bounds.append(width_ranges['gamma'][0])
            upper_bounds.append(width_ranges['gamma'][1])

        # Baseline bounds
        lower_bounds.append(baseline_range[0])
        upper_bounds.append(baseline_range[1])

        return lower_bounds, upper_bounds


class StagedFittingStrategy:
    """
    3-stage fitting strategy with progressive parameter release.

    Implements warm-start optimization where each stage uses the previous
    stage's result as the initial guess, progressively releasing constraints
    on parameter values.
    """

    def __init__(self, config: Optional[Dict] = None):
        """
        Initialize staged fitting strategy.

        Parameters:
        -----------
        config : Optional[Dict]
            Configuration dictionary with stage-specific settings
        """
        self.config = config or {}

        # Default stage configurations
        self.stage_configs = {
            'stage1': {
                'name': 'Fixed Positions',
                'position_lock': True,
                'max_iterations': 1000,
                'ftol': 1e-6,
                'xtol': 1e-6,
                'gtol': 1e-6
            },
            'stage2': {
                'name': 'Constrained Positions',
                'position_tolerance': 0.08,  # ±8% default
                'max_iterations': 1500,
                'ftol': 1e-7,
                'xtol': 1e-7,
                'gtol': 1e-7
            },
            'stage3': {
                'name': 'Full Freedom',
                'full_freedom': True,
                'max_iterations': 2000,
                'ftol': 1e-8,
                'xtol': 1e-8,
                'gtol': 1e-8
            }
        }

        # Update with user configuration
        for stage_key in ['stage1', 'stage2', 'stage3']:
            if stage_key in self.config:
                self.stage_configs[stage_key].update(self.config[stage_key])

        # Quality thresholds
        self.quality_thresholds = {
            'excellent': 0.95,
            'good': 0.85,
            'fair': 0.70,
            'poor': 0.0
        }

        # Convergence tracking
        self.convergence_history = []

    def fit_staged(self,
                   x_data: np.ndarray,
                   y_data: np.ndarray,
                   peak_centers: List[float],
                   model_func: Callable,
                   initial_params: Optional[np.ndarray] = None,
                   amplitude_range: Tuple[float, float] = (0, np.inf),
                   width_ranges: Optional[Dict[str, Tuple[float, float]]] = None,
                   baseline_range: Tuple[float, float] = (-np.inf, np.inf),
                   nucleus_type: str = '1H',
                   verbose: bool = False) -> StagedFittingResult:
        """
        Execute 3-stage fitting with warm starts.

        Parameters:
        -----------
        x_data : np.ndarray
            Independent variable data
        y_data : np.ndarray
            Dependent variable data
        peak_centers : List[float]
            Initial peak center positions
        model_func : Callable
            Model function with signature model_func(x, *params)
        initial_params : Optional[np.ndarray]
            Initial parameter guess. If None, will be estimated.
        amplitude_range : Tuple[float, float]
            Bounds for peak amplitudes
        width_ranges : Optional[Dict[str, Tuple[float, float]]]
            Bounds for sigma and gamma widths
        baseline_range : Tuple[float, float]
            Bounds for baseline parameter
        nucleus_type : str
            Nucleus type for parameter estimation
        verbose : bool
            Print detailed progress information

        Returns:
        --------
        StagedFittingResult : Complete fitting results with stage history
        """
        start_time = time.time()
        result = StagedFittingResult()
        result.convergence_trajectory = {
            'r_squared': [],
            'residual_std': [],
            'chi_squared': []
        }

        n_peaks = len(peak_centers)
        n_data = len(x_data)

        if verbose:
            print(f"\n{'='*60}")
            print(f"Staged Fitting Strategy: {n_peaks} peaks, {n_data} data points")
            print(f"{'='*60}")

        # Prepare initial parameters if not provided
        if initial_params is None:
            initial_params = self._estimate_initial_parameters(
                x_data, y_data, peak_centers, nucleus_type
            )

        # Validate initial parameters
        expected_n_params = n_peaks * 4 + 1  # amp, center, sigma, gamma per peak + baseline
        if len(initial_params) != expected_n_params:
            result.success = False
            result.stage_history.append(StageResult(
                stage=0,
                success=False,
                message=f"Parameter count mismatch: expected {expected_n_params}, got {len(initial_params)}"
            ))
            return result

        # Stage 1: Fixed positions
        stage1_result = self._fit_stage_1(
            x_data, y_data, initial_params, peak_centers, model_func,
            amplitude_range, width_ranges, baseline_range, nucleus_type, verbose
        )
        result.stage_history.append(stage1_result)
        result.convergence_trajectory['r_squared'].append(stage1_result.r_squared)
        result.convergence_trajectory['residual_std'].append(stage1_result.residual_std)
        result.convergence_trajectory['chi_squared'].append(stage1_result.reduced_chi_squared)

        if not stage1_result.success:
            if verbose:
                print(f"  Stage 1 failed: {stage1_result.message}")
            result.success = False
            result.total_time = time.time() - start_time
            return result

        # Use Stage 1 result as warm start for Stage 2
        p0_stage2 = stage1_result.parameters

        # Stage 2: Constrained positions
        stage2_result = self._fit_stage_2(
            x_data, y_data, p0_stage2, peak_centers, model_func,
            amplitude_range, width_ranges, baseline_range, nucleus_type, verbose
        )
        result.stage_history.append(stage2_result)
        result.convergence_trajectory['r_squared'].append(stage2_result.r_squared)
        result.convergence_trajectory['residual_std'].append(stage2_result.residual_std)
        result.convergence_trajectory['chi_squared'].append(stage2_result.reduced_chi_squared)

        if not stage2_result.success:
            if verbose:
                print(f"  Stage 2 failed, using Stage 1 result")
            # Fall back to Stage 1 result
            result.final_parameters = stage1_result.parameters
            result.final_covariance = stage1_result.covariance
            result.final_r_squared = stage1_result.r_squared
            result.final_reduced_chi_squared = stage1_result.reduced_chi_squared
            result.success = True
            result.best_stage = 1
        else:
            # Use Stage 2 result as warm start for Stage 3
            p0_stage3 = stage2_result.parameters

            # Stage 3: Full freedom
            stage3_result = self._fit_stage_3(
                x_data, y_data, p0_stage3, peak_centers, model_func,
                amplitude_range, width_ranges, baseline_range, nucleus_type, verbose
            )
            result.stage_history.append(stage3_result)
            result.convergence_trajectory['r_squared'].append(stage3_result.r_squared)
            result.convergence_trajectory['residual_std'].append(stage3_result.residual_std)
            result.convergence_trajectory['chi_squared'].append(stage3_result.reduced_chi_squared)

            if not stage3_result.success:
                if verbose:
                    print(f"  Stage 3 failed, using Stage 2 result")
                # Fall back to Stage 2 result
                result.final_parameters = stage2_result.parameters
                result.final_covariance = stage2_result.covariance
                result.final_r_squared = stage2_result.r_squared
                result.final_reduced_chi_squared = stage2_result.reduced_chi_squared
                result.success = True
                result.best_stage = 2
            else:
                # Stage 3 succeeded - select best result
                best_stage = self._select_best_stage(result.stage_history)
                result.best_stage = best_stage
                best_result = result.stage_history[best_stage - 1]

                result.final_parameters = best_result.parameters
                result.final_covariance = best_result.covariance
                result.final_r_squared = best_result.r_squared
                result.final_reduced_chi_squared = best_result.reduced_chi_squared
                result.success = True

        # Quality assessment
        result.quality_assessment = self._assess_quality(result, n_peaks, n_data)

        result.total_time = time.time() - start_time

        if verbose:
            print(f"\n{'='*60}")
            print(f"Staged Fitting Complete:")
            print(f"  Best stage: {result.best_stage}")
            print(f"  Final R²: {result.final_r_squared:.6f}")
            print(f"  Final χ²_red: {result.final_reduced_chi_squared:.6f}")
            print(f"  Quality: {result.quality_assessment.get('overall_quality', 'unknown')}")
            print(f"  Total time: {result.total_time:.3f}s")
            print(f"{'='*60}\n")

        return result

    def _fit_stage_1(self,
                     x_data: np.ndarray,
                     y_data: np.ndarray,
                     p0: np.ndarray,
                     peak_centers: List[float],
                     model_func: Callable,
                     amplitude_range: Tuple[float, float],
                     width_ranges: Optional[Dict],
                     baseline_range: Tuple[float, float],
                     nucleus_type: str,
                     verbose: bool) -> StageResult:
        """
        Stage 1: Fixed positions, optimize amplitudes and widths.

        This stage locks peak positions at their initial values and optimizes
        only the amplitude, sigma, gamma, and baseline parameters.
        """
        stage = 1
        config = self.stage_configs['stage1']

        if verbose:
            print(f"\nStage {stage}: {config['name']}")
            print(f"  Locking positions, optimizing amplitudes & widths")

        result = StageResult(stage=stage)
        start_time = time.time()

        try:
            n_peaks = len(peak_centers)

            # Generate bounds with locked positions
            lower_bounds, upper_bounds = BoundsGenerator.fixed_positions(
                peak_centers=peak_centers,
                n_peaks=n_peaks,
                amplitude_range=amplitude_range,
                width_ranges=width_ranges,
                baseline_range=baseline_range,
                nucleus_type=nucleus_type
            )

            bounds = (lower_bounds, upper_bounds)

            # Perform curve fitting
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", OptimizeWarning)
                warnings.simplefilter("ignore", RuntimeWarning)

                popt, pcov = curve_fit(
                    model_func,
                    x_data,
                    y_data,
                    p0=p0,
                    bounds=bounds,
                    maxfev=config['max_iterations'],
                    ftol=config['ftol'],
                    xtol=config['xtol'],
                    gtol=config['gtol'],
                    method='trf'
                )

            result.parameters = popt
            result.covariance = pcov
            result.success = True

            # Validate result
            validation = self._validate_stage_result(
                popt, pcov, x_data, y_data, model_func, stage, n_peaks
            )

            result.r_squared = validation['r_squared']
            result.reduced_chi_squared = validation['reduced_chi_squared']
            result.residual_std = validation['residual_std']
            result.convergence_quality = validation['quality']
            result.message = "Stage 1 completed successfully"

            if verbose:
                print(f"  ✓ R² = {result.r_squared:.6f}, χ²_red = {result.reduced_chi_squared:.6f}")

        except Exception as e:
            result.success = False
            result.message = f"Stage 1 error: {str(e)}"
            if verbose:
                print(f"  ✗ {result.message}")

        result.fit_time = time.time() - start_time
        return result

    def _fit_stage_2(self,
                     x_data: np.ndarray,
                     y_data: np.ndarray,
                     p0_warm: np.ndarray,
                     peak_centers: List[float],
                     model_func: Callable,
                     amplitude_range: Tuple[float, float],
                     width_ranges: Optional[Dict],
                     baseline_range: Tuple[float, float],
                     nucleus_type: str,
                     verbose: bool) -> StageResult:
        """
        Stage 2: Constrained positions (warm start from Stage 1).

        Allows peak positions to vary within a tight tolerance (±5-10%)
        while continuing to optimize all other parameters.
        """
        stage = 2
        config = self.stage_configs['stage2']

        if verbose:
            tolerance_pct = config['position_tolerance'] * 100
            print(f"\nStage {stage}: {config['name']}")
            print(f"  Positions ±{tolerance_pct:.1f}%, warm start from Stage 1")

        result = StageResult(stage=stage)
        start_time = time.time()

        try:
            n_peaks = len(peak_centers)
            x_range = (np.min(x_data), np.max(x_data))

            # Extract current positions from warm start parameters
            current_centers = [p0_warm[i*4 + 1] for i in range(n_peaks)]

            # Generate constrained bounds
            lower_bounds, upper_bounds = BoundsGenerator.constrained_positions(
                peak_centers=current_centers,
                tolerance=config['position_tolerance'],
                n_peaks=n_peaks,
                amplitude_range=amplitude_range,
                width_ranges=width_ranges,
                baseline_range=baseline_range,
                nucleus_type=nucleus_type,
                x_range=x_range
            )

            bounds = (lower_bounds, upper_bounds)

            # Perform curve fitting with warm start
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", OptimizeWarning)
                warnings.simplefilter("ignore", RuntimeWarning)

                popt, pcov = curve_fit(
                    model_func,
                    x_data,
                    y_data,
                    p0=p0_warm,  # Warm start!
                    bounds=bounds,
                    maxfev=config['max_iterations'],
                    ftol=config['ftol'],
                    xtol=config['xtol'],
                    gtol=config['gtol'],
                    method='trf'
                )

            result.parameters = popt
            result.covariance = pcov
            result.success = True

            # Validate result
            validation = self._validate_stage_result(
                popt, pcov, x_data, y_data, model_func, stage, n_peaks
            )

            result.r_squared = validation['r_squared']
            result.reduced_chi_squared = validation['reduced_chi_squared']
            result.residual_std = validation['residual_std']
            result.convergence_quality = validation['quality']
            result.message = "Stage 2 completed successfully"

            if verbose:
                print(f"  ✓ R² = {result.r_squared:.6f}, χ²_red = {result.reduced_chi_squared:.6f}")

        except Exception as e:
            result.success = False
            result.message = f"Stage 2 error: {str(e)}"
            if verbose:
                print(f"  ✗ {result.message}")

        result.fit_time = time.time() - start_time
        return result

    def _fit_stage_3(self,
                     x_data: np.ndarray,
                     y_data: np.ndarray,
                     p0_warm: np.ndarray,
                     peak_centers: List[float],
                     model_func: Callable,
                     amplitude_range: Tuple[float, float],
                     width_ranges: Optional[Dict],
                     baseline_range: Tuple[float, float],
                     nucleus_type: str,
                     verbose: bool) -> StageResult:
        """
        Stage 3: Full optimization (warm start from Stage 2).

        Allows all parameters to vary freely within wide bounds,
        enabling full exploration of parameter space.
        """
        stage = 3
        config = self.stage_configs['stage3']

        if verbose:
            print(f"\nStage {stage}: {config['name']}")
            print(f"  Full parameter freedom, warm start from Stage 2")

        result = StageResult(stage=stage)
        start_time = time.time()

        try:
            n_peaks = len(peak_centers)
            x_range = (np.min(x_data), np.max(x_data))

            # Extract current positions from warm start parameters
            current_centers = [p0_warm[i*4 + 1] for i in range(n_peaks)]

            # Generate wide bounds for full freedom
            lower_bounds, upper_bounds = BoundsGenerator.full_freedom(
                peak_centers=current_centers,
                x_range=x_range,
                n_peaks=n_peaks,
                amplitude_range=amplitude_range,
                width_ranges=width_ranges,
                baseline_range=baseline_range,
                nucleus_type=nucleus_type
            )

            bounds = (lower_bounds, upper_bounds)

            # Perform curve fitting with warm start
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", OptimizeWarning)
                warnings.simplefilter("ignore", RuntimeWarning)

                popt, pcov = curve_fit(
                    model_func,
                    x_data,
                    y_data,
                    p0=p0_warm,  # Warm start!
                    bounds=bounds,
                    maxfev=config['max_iterations'],
                    ftol=config['ftol'],
                    xtol=config['xtol'],
                    gtol=config['gtol'],
                    method='trf'
                )

            result.parameters = popt
            result.covariance = pcov
            result.success = True

            # Validate result
            validation = self._validate_stage_result(
                popt, pcov, x_data, y_data, model_func, stage, n_peaks
            )

            result.r_squared = validation['r_squared']
            result.reduced_chi_squared = validation['reduced_chi_squared']
            result.residual_std = validation['residual_std']
            result.convergence_quality = validation['quality']
            result.message = "Stage 3 completed successfully"

            if verbose:
                print(f"  ✓ R² = {result.r_squared:.6f}, χ²_red = {result.reduced_chi_squared:.6f}")

        except Exception as e:
            result.success = False
            result.message = f"Stage 3 error: {str(e)}"
            if verbose:
                print(f"  ✗ {result.message}")

        result.fit_time = time.time() - start_time
        return result

    def _validate_stage_result(self,
                               params: np.ndarray,
                               covariance: np.ndarray,
                               x_data: np.ndarray,
                               y_data: np.ndarray,
                               model_func: Callable,
                               stage: int,
                               n_peaks: int) -> Dict:
        """
        Validate stage convergence and quality.

        Parameters:
        -----------
        params : np.ndarray
            Fitted parameters
        covariance : np.ndarray
            Covariance matrix
        x_data : np.ndarray
            Independent variable data
        y_data : np.ndarray
            Dependent variable data
        model_func : Callable
            Model function
        stage : int
            Stage number
        n_peaks : int
            Number of peaks

        Returns:
        --------
        Dict : Validation results with quality metrics
        """
        validation = {
            'r_squared': 0.0,
            'reduced_chi_squared': np.inf,
            'residual_std': np.inf,
            'quality': 'poor',
            'warnings': []
        }

        try:
            # Calculate residuals
            y_pred = model_func(x_data, *params)
            residuals = y_data - y_pred

            # R-squared
            ss_res = np.sum(residuals**2)
            ss_tot = np.sum((y_data - np.mean(y_data))**2)
            r_squared = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0.0

            # Reduced chi-squared
            n_data = len(x_data)
            n_params = len(params)
            dof = max(n_data - n_params, 1)
            reduced_chi_squared = ss_res / dof

            # Residual standard deviation
            residual_std = np.std(residuals)

            validation['r_squared'] = r_squared
            validation['reduced_chi_squared'] = reduced_chi_squared
            validation['residual_std'] = residual_std

            # Quality assessment
            if r_squared >= self.quality_thresholds['excellent']:
                validation['quality'] = 'excellent'
            elif r_squared >= self.quality_thresholds['good']:
                validation['quality'] = 'good'
            elif r_squared >= self.quality_thresholds['fair']:
                validation['quality'] = 'fair'
            else:
                validation['quality'] = 'poor'

            # Check for potential issues
            if covariance is not None:
                # Check for infinite uncertainties
                diag = np.diag(covariance)
                if np.any(np.isinf(diag)) or np.any(diag < 0):
                    validation['warnings'].append('Infinite or negative parameter uncertainties')

            # Check for unrealistic parameter values
            for i in range(n_peaks):
                amp = params[i*4]
                sigma = params[i*4 + 2]
                gamma = params[i*4 + 3]

                if amp <= 0:
                    validation['warnings'].append(f'Peak {i+1}: non-positive amplitude')
                if sigma <= 0:
                    validation['warnings'].append(f'Peak {i+1}: non-positive sigma')
                if gamma < 0:
                    validation['warnings'].append(f'Peak {i+1}: negative gamma')

        except Exception as e:
            validation['warnings'].append(f'Validation error: {str(e)}')

        return validation

    def _estimate_initial_parameters(self,
                                     x_data: np.ndarray,
                                     y_data: np.ndarray,
                                     peak_centers: List[float],
                                     nucleus_type: str) -> np.ndarray:
        """
        Estimate initial parameters for multi-peak Voigt fitting.

        Parameters:
        -----------
        x_data : np.ndarray
            Independent variable data
        y_data : np.ndarray
            Dependent variable data
        peak_centers : List[float]
            Initial peak center positions
        nucleus_type : str
            Nucleus type for parameter estimation

        Returns:
        --------
        np.ndarray : Initial parameter estimates
        """
        n_peaks = len(peak_centers)
        params = []

        # Get typical width for nucleus type
        typical_width = BoundsGenerator.DEFAULT_RANGES.get(
            nucleus_type, BoundsGenerator.DEFAULT_RANGES['default']
        )['typical_width']

        # Estimate baseline as median of lowest 20% of data
        sorted_y = np.sort(y_data)
        n_baseline = max(int(len(y_data) * 0.2), 1)
        baseline = np.median(sorted_y[:n_baseline])

        # Estimate parameters for each peak
        for center in peak_centers:
            # Find closest data point to center
            idx = np.argmin(np.abs(x_data - center))

            # Amplitude: peak height above baseline
            amplitude = max(y_data[idx] - baseline, 1.0)

            # Width estimates
            sigma = typical_width
            gamma = typical_width * 0.5

            params.extend([amplitude, center, sigma, gamma])

        # Add baseline
        params.append(baseline)

        return np.array(params)

    def _select_best_stage(self, stage_history: List[StageResult]) -> int:
        """
        Select the best stage based on quality metrics.

        Prioritizes later stages if they provide significant improvement,
        but falls back to earlier stages if later stages don't converge well.

        Parameters:
        -----------
        stage_history : List[StageResult]
            Results from all stages

        Returns:
        --------
        int : Best stage number (1, 2, or 3)
        """
        successful_stages = [s for s in stage_history if s.success]

        if not successful_stages:
            return 1  # Fallback

        # Find stage with highest R²
        best_stage = max(successful_stages, key=lambda s: s.r_squared)

        # Check if later stages provide significant improvement
        improvement_threshold = 0.01  # 1% improvement threshold

        for i, stage in enumerate(successful_stages[:-1]):
            next_stage = successful_stages[i + 1]
            improvement = next_stage.r_squared - stage.r_squared

            # If later stage doesn't improve significantly, use current
            if improvement < improvement_threshold:
                return stage.stage

        return best_stage.stage

    def _assess_quality(self,
                       result: StagedFittingResult,
                       n_peaks: int,
                       n_data: int) -> Dict:
        """
        Comprehensive quality assessment of staged fitting.

        Parameters:
        -----------
        result : StagedFittingResult
            Complete fitting result
        n_peaks : int
            Number of peaks
        n_data : int
            Number of data points

        Returns:
        --------
        Dict : Quality assessment metrics
        """
        assessment = {
            'overall_quality': 'poor',
            'stage_progression': [],
            'convergence_stable': False,
            'parameter_quality': 'unknown',
            'recommendations': []
        }

        # Assess stage progression
        r2_values = [s.r_squared for s in result.stage_history if s.success]
        if len(r2_values) > 1:
            # Check if R² improved through stages
            improvements = [r2_values[i+1] - r2_values[i] for i in range(len(r2_values)-1)]
            assessment['stage_progression'] = improvements

            # Check convergence stability
            if all(imp >= -0.01 for imp in improvements):  # Allow small degradation
                assessment['convergence_stable'] = True

        # Overall quality based on final R²
        final_r2 = result.final_r_squared
        if final_r2 >= self.quality_thresholds['excellent']:
            assessment['overall_quality'] = 'excellent'
        elif final_r2 >= self.quality_thresholds['good']:
            assessment['overall_quality'] = 'good'
        elif final_r2 >= self.quality_thresholds['fair']:
            assessment['overall_quality'] = 'fair'
        else:
            assessment['overall_quality'] = 'poor'
            assessment['recommendations'].append('Consider reviewing peak detection')

        # Check parameter quality
        if result.final_covariance is not None:
            diag = np.diag(result.final_covariance)
            if np.all(np.isfinite(diag)) and np.all(diag >= 0):
                assessment['parameter_quality'] = 'good'
            else:
                assessment['parameter_quality'] = 'poor'
                assessment['recommendations'].append('Parameter uncertainties may be unreliable')

        # Check degrees of freedom
        n_params = n_peaks * 4 + 1
        dof = n_data - n_params
        if dof < 10:
            assessment['recommendations'].append(f'Low degrees of freedom ({dof}), results may be unstable')

        return assessment


# Convenience function for quick usage
def fit_with_staged_strategy(x_data: np.ndarray,
                             y_data: np.ndarray,
                             peak_centers: List[float],
                             model_func: Callable,
                             nucleus_type: str = '1H',
                             config: Optional[Dict] = None,
                             verbose: bool = False) -> StagedFittingResult:
    """
    Convenience function for staged fitting.

    Parameters:
    -----------
    x_data : np.ndarray
        Independent variable data
    y_data : np.ndarray
        Dependent variable data
    peak_centers : List[float]
        Initial peak center positions
    model_func : Callable
        Model function with signature model_func(x, *params)
    nucleus_type : str
        Nucleus type for parameter estimation
    config : Optional[Dict]
        Configuration dictionary
    verbose : bool
        Print detailed progress

    Returns:
    --------
    StagedFittingResult : Complete fitting results
    """
    strategy = StagedFittingStrategy(config=config)

    # Estimate amplitude range from data
    y_min = np.min(y_data)
    y_max = np.max(y_data)
    amplitude_range = (0, (y_max - y_min) * 2)

    # Estimate baseline range
    baseline_range = (y_min - abs(y_min) * 0.5, y_max * 0.5)

    return strategy.fit_staged(
        x_data=x_data,
        y_data=y_data,
        peak_centers=peak_centers,
        model_func=model_func,
        amplitude_range=amplitude_range,
        baseline_range=baseline_range,
        nucleus_type=nucleus_type,
        verbose=verbose
    )


if __name__ == '__main__':
    """
    Example usage and testing
    """
    from scipy.special import wofz

    def voigt_profile(x, amplitude, center, sigma, gamma, baseline=0):
        """Single Voigt profile"""
        sigma = max(sigma, 1e-6)
        z = ((x - center) + 1j*gamma) / (sigma * np.sqrt(2))
        voigt = amplitude * np.real(wofz(z)) / (sigma * np.sqrt(2*np.pi))
        return voigt + baseline

    def multi_voigt_model(x, *params):
        """Multi-peak Voigt model"""
        *peak_params, baseline = params
        n_peaks = len(peak_params) // 4
        y_model = np.full(len(x), baseline)

        for i in range(n_peaks):
            p_start = i * 4
            amplitude, center, sigma, gamma = peak_params[p_start:p_start+4]
            y_model += voigt_profile(x, amplitude, center, sigma, gamma, 0)

        return y_model

    # Generate synthetic data with 2 peaks
    np.random.seed(42)
    x_data = np.linspace(7.0, 9.0, 300)

    # True parameters
    true_params = [
        100.0, 7.5, 0.02, 0.01,  # Peak 1
        80.0, 8.5, 0.015, 0.012,  # Peak 2
        5.0  # Baseline
    ]

    y_true = multi_voigt_model(x_data, *true_params)
    noise = np.random.normal(0, 2, len(x_data))
    y_data = y_true + noise

    # Test staged fitting
    peak_centers = [7.5, 8.5]

    print("Testing Staged Fitting Strategy")
    print("=" * 60)

    result = fit_with_staged_strategy(
        x_data=x_data,
        y_data=y_data,
        peak_centers=peak_centers,
        model_func=multi_voigt_model,
        nucleus_type='1H',
        verbose=True
    )

    if result.success:
        print("\nFitted Parameters:")
        print(f"  Peak 1: amp={result.final_parameters[0]:.2f}, "
              f"center={result.final_parameters[1]:.4f}, "
              f"sigma={result.final_parameters[2]:.4f}, "
              f"gamma={result.final_parameters[3]:.4f}")
        print(f"  Peak 2: amp={result.final_parameters[4]:.2f}, "
              f"center={result.final_parameters[5]:.4f}, "
              f"sigma={result.final_parameters[6]:.4f}, "
              f"gamma={result.final_parameters[7]:.4f}")
        print(f"  Baseline: {result.final_parameters[8]:.2f}")

        print("\nTrue Parameters:")
        print(f"  Peak 1: amp={true_params[0]:.2f}, "
              f"center={true_params[1]:.4f}, "
              f"sigma={true_params[2]:.4f}, "
              f"gamma={true_params[3]:.4f}")
        print(f"  Peak 2: amp={true_params[4]:.2f}, "
              f"center={true_params[5]:.4f}, "
              f"sigma={true_params[6]:.4f}, "
              f"gamma={true_params[7]:.4f}")
        print(f"  Baseline: {true_params[8]:.2f}")
    else:
        print("\nFitting failed!")

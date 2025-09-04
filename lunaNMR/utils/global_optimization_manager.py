"""
Global NMR Peak Optimization Manager

This module implements a sophisticated two-phase optimization system for NMR peak fitting:
1. Phase 1 (Survey): Complete first-pass fitting to extract linewidth statistics
2. Phase 2 (Refinement): Iterative improvement of failed peaks using constraints

Author: Guillaume Mas
Date: 2025
"""

import numpy as np
from typing import Dict, List, Tuple, Optional, Any
from dataclasses import dataclass
from enum import Enum
import logging
from collections import defaultdict

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class PeakQuality(Enum):
    """Peak fitting quality classification"""
    EXCELLENT = "excellent"  # R² > 0.95
    GOOD = "good"           # 0.85 ≤ R² ≤ 0.95
    POOR = "poor"           # 0.5 ≤ R² < 0.85
    FAILED = "failed"       # R² < 0.5


@dataclass
class PeakResult:
    """Container for peak fitting results with quality assessment"""
    peak_id: str
    x_position: float
    y_position: float
    x_r_squared: float
    y_r_squared: float
    x_linewidth: Optional[float] = None
    y_linewidth: Optional[float] = None
    quality: Optional[PeakQuality] = None
    fitting_method: str = "standard"
    optimization_round: int = 0
    error_message: Optional[str] = None

    def __post_init__(self):
        """Automatically classify peak quality based on R² values"""
        if self.quality is None:
            # Use the worse of X and Y dimensions for overall quality
            min_r_squared = min(self.x_r_squared, self.y_r_squared)

            if min_r_squared >= 0.95:
                self.quality = PeakQuality.EXCELLENT
            elif min_r_squared >= 0.85:
                self.quality = PeakQuality.GOOD
            elif min_r_squared >= 0.5:
                self.quality = PeakQuality.POOR
            else:
                self.quality = PeakQuality.FAILED

    @property
    def is_successful(self) -> bool:
        """Check if peak fitting was successful (Good or Excellent)"""
        return self.quality in [PeakQuality.EXCELLENT, PeakQuality.GOOD]

    @property
    def needs_improvement(self) -> bool:
        """Check if peak needs improvement (Poor or Failed)"""
        return self.quality in [PeakQuality.POOR, PeakQuality.FAILED]


class LinewidthStatistics:
    """
    Robust statistics collector for linewidth data from successful peak fits

    Maintains separate statistics for X and Y dimensions with robust methods
    to handle outliers and provide reliable constraints for failed peaks.
    """

    def __init__(self):
        self.x_linewidths: List[float] = []
        self.y_linewidths: List[float] = []
        self.peak_count = 0
        self.last_updated_round = 0

    def add_successful_peak(self, peak_result: PeakResult) -> None:
        """Add linewidth data from a successful peak fit"""
        if peak_result.is_successful:
            if peak_result.x_linewidth is not None:
                self.x_linewidths.append(peak_result.x_linewidth)
            if peak_result.y_linewidth is not None:
                self.y_linewidths.append(peak_result.y_linewidth)
            self.peak_count += 1

    def get_x_constraints(self) -> Tuple[float, float]:
        """Get robust X-dimension linewidth constraints (median ± 2×IQR)"""
        if len(self.x_linewidths) < 3:
            # Fallback to typical 1H linewidths
            return (0.005, 0.100)

        x_array = np.array(self.x_linewidths)
        median = np.median(x_array)
        q75, q25 = np.percentile(x_array, [75, 25])
        iqr = q75 - q25

        # Robust bounds: median ± 2×IQR
        lower = max(0.001, median - 2 * iqr)  # Minimum physical limit
        upper = median + 2 * iqr

        return (lower, upper)

    def get_y_constraints(self) -> Tuple[float, float]:
        """Get robust Y-dimension linewidth constraints (median ± 2×IQR)"""
        if len(self.y_linewidths) < 3:
            # Fallback to typical 13C/15N linewidths
            return (0.5, 8.0)

        y_array = np.array(self.y_linewidths)
        median = np.median(y_array)
        q75, q25 = np.percentile(y_array, [75, 25])
        iqr = q75 - q25

        # Robust bounds: median ± 2×IQR
        lower = max(0.1, median - 2 * iqr)  # Minimum physical limit
        upper = median + 2 * iqr

        return (lower, upper)

    def get_statistics_summary(self) -> Dict[str, Any]:
        """Get comprehensive statistics summary"""
        x_constraints = self.get_x_constraints()
        y_constraints = self.get_y_constraints()

        return {
            'peak_count': self.peak_count,
            'x_linewidths': {
                'count': len(self.x_linewidths),
                'median': np.median(self.x_linewidths) if self.x_linewidths else None,
                'iqr': np.percentile(self.x_linewidths, [75, 25]) if len(self.x_linewidths) >= 2 else None,
                'constraints': x_constraints
            },
            'y_linewidths': {
                'count': len(self.y_linewidths),
                'median': np.median(self.y_linewidths) if self.y_linewidths else None,
                'iqr': np.percentile(self.y_linewidths, [75, 25]) if len(self.y_linewidths) >= 2 else None,
                'constraints': y_constraints
            }
        }


@dataclass
class OptimizationRound:
    """Track results and progress for each optimization round"""
    round_number: int
    peaks_processed: int = 0
    peaks_improved: int = 0
    quality_distribution: Dict[PeakQuality, int] = None
    average_r_squared_improvement: float = 0.0
    processing_time: float = 0.0

    def __post_init__(self):
        if self.quality_distribution is None:
            self.quality_distribution = {
                PeakQuality.EXCELLENT: 0,
                PeakQuality.GOOD: 0,
                PeakQuality.POOR: 0,
                PeakQuality.FAILED: 0
            }

    @property
    def improvement_percentage(self) -> float:
        """Calculate percentage of peaks that improved this round"""
        if self.peaks_processed == 0:
            return 0.0
        return (self.peaks_improved / self.peaks_processed) * 100.0

    @property
    def success_rate(self) -> float:
        """Calculate percentage of successful peaks (Excellent + Good)"""
        total_successful = (self.quality_distribution[PeakQuality.EXCELLENT] +
                          self.quality_distribution[PeakQuality.GOOD])
        if self.peaks_processed == 0:
            return 0.0
        return (total_successful / self.peaks_processed) * 100.0


class GlobalOptimizationManager:
    """
    Two-Phase Global Peak Optimization Manager

    This class orchestrates the complete optimization workflow:

    Phase 1 (Survey):
    - Fits all peaks with standard methods
    - Collects linewidth statistics from successful peaks
    - Identifies peaks needing improvement

    Phase 2 (Iterative Refinement):
    - Applies linewidth constraints from successful peaks
    - Uses multi-peak deconvolution for complex cases
    - Continues until convergence (<5% improvement per round)
    """

    def __init__(self, convergence_threshold: float = 0.05, max_rounds: int = 5):
        """
        Initialize the Global Optimization Manager

        Args:
            convergence_threshold: Stop when <5% of peaks improve per round
            max_rounds: Maximum number of refinement rounds
        """
        self.convergence_threshold = convergence_threshold
        self.max_rounds = max_rounds

        # Results storage
        self.peak_results: Dict[str, PeakResult] = {}
        self.optimization_rounds: List[OptimizationRound] = []

        # Statistics and constraints
        self.linewidth_stats = LinewidthStatistics()

        # Progress tracking
        self.total_peaks = 0
        self.phase = "not_started"
        self.current_round = 0

        logger.info("🚀 GlobalOptimizationManager initialized")
        logger.info(f"   Convergence threshold: {convergence_threshold*100:.1f}%")
        logger.info(f"   Maximum rounds: {max_rounds}")

    def optimize_peak_list(self, peak_list: List[Tuple[float, float, str]],
                          integrator_instance) -> Dict[str, Any]:
        """
        Main optimization workflow - orchestrates both phases

        Args:
            peak_list: List of (x_pos, y_pos, peak_id) tuples
            integrator_instance: VoigtIntegrator instance for fitting

        Returns:
            Complete optimization report with all results and statistics
        """
        logger.info("🎯 Starting Two-Phase Global Optimization")
        logger.info(f"   Processing {len(peak_list)} peaks")

        self.total_peaks = len(peak_list)

        # Phase 1: Survey - fit all peaks and collect statistics
        logger.info("\n" + "="*60)
        logger.info("📊 PHASE 1: COMPLETE SURVEY")
        logger.info("="*60)

        self.phase = "survey"
        phase1_results = self._phase1_survey(peak_list, integrator_instance)

        # Phase 2: Iterative refinement of failed peaks
        logger.info("\n" + "="*60)
        logger.info("🔄 PHASE 2: ITERATIVE REFINEMENT")
        logger.info("="*60)

        self.phase = "refinement"
        phase2_results = self._phase2_iterative_refinement(integrator_instance)

        # Generate comprehensive report
        final_report = self._generate_final_report()

        logger.info("\n🎉 Global Optimization Complete!")
        self._print_final_summary()

        return final_report

    def _phase1_survey(self, peak_list: List[Tuple[float, float, str]],
                      integrator_instance) -> Dict[str, Any]:
        """
        Phase 1: Complete survey of all peaks with standard fitting
        """
        import time
        start_time = time.time()

        logger.info(f"Fitting {len(peak_list)} peaks with standard methods...")

        # Process each peak with standard fitting
        for i, (x_pos, y_pos, peak_id) in enumerate(peak_list, 1):
            logger.info(f"  Peak {i}/{len(peak_list)}: {peak_id} at ({x_pos:.3f}, {y_pos:.1f})")

            try:
                # Standard fitting (no dynamic optimization)
                result = integrator_instance.fit_peak_voigt_2d(
                    x_pos, y_pos, peak_id,
                    use_dynamic_optimization=False  # Phase 1 uses standard fitting
                )

                if result and result.get('success', False):
                    # Extract linewidth information from fitting results
                    x_linewidth = self._extract_linewidth_from_result(result, 'x')
                    y_linewidth = self._extract_linewidth_from_result(result, 'y')

                    # Extract fitting results
                    peak_result = PeakResult(
                        peak_id=peak_id,
                        x_position=x_pos,
                        y_position=y_pos,
                        x_r_squared=result.get('x_r_squared', 0.0),
                        y_r_squared=result.get('y_r_squared', 0.0),
                        x_linewidth=x_linewidth,
                        y_linewidth=y_linewidth,
                        fitting_method="standard",
                        optimization_round=1
                    )

                    self.peak_results[peak_id] = peak_result

                    # Collect statistics from successful peaks
                    if peak_result.is_successful:
                        self.linewidth_stats.add_successful_peak(peak_result)
                        logger.info(f"    ✅ Success: R²=({peak_result.x_r_squared:.3f}, {peak_result.y_r_squared:.3f})")
                    else:
                        logger.info(f"    ⚠️  Needs improvement: R²=({peak_result.x_r_squared:.3f}, {peak_result.y_r_squared:.3f})")

                else:
                    # Fitting failed completely
                    peak_result = PeakResult(
                        peak_id=peak_id,
                        x_position=x_pos,
                        y_position=y_pos,
                        x_r_squared=0.0,
                        y_r_squared=0.0,
                        fitting_method="standard",
                        optimization_round=1,
                        error_message="Standard fitting failed"
                    )
                    self.peak_results[peak_id] = peak_result
                    logger.info(f"    ❌ Failed: Standard fitting unsuccessful")

            except Exception as e:
                logger.error(f"    💥 Exception during fitting: {e}")
                peak_result = PeakResult(
                    peak_id=peak_id,
                    x_position=x_pos,
                    y_position=y_pos,
                    x_r_squared=0.0,
                    y_r_squared=0.0,
                    fitting_method="standard",
                    optimization_round=1,
                    error_message=str(e)
                )
                self.peak_results[peak_id] = peak_result

        # Create Phase 1 summary
        processing_time = time.time() - start_time
        round1 = OptimizationRound(round_number=1, processing_time=processing_time)

        # Count quality distribution
        for peak_result in self.peak_results.values():
            round1.quality_distribution[peak_result.quality] += 1
            round1.peaks_processed += 1

        self.optimization_rounds.append(round1)

        # Print Phase 1 results
        self._print_phase1_summary()

        return {"round": round1, "linewidth_stats": self.linewidth_stats.get_statistics_summary()}

    def _phase2_iterative_refinement(self, integrator_instance) -> Dict[str, Any]:
        """
        Phase 2: Iterative refinement of failed peaks using constraints
        """
        # Identify peaks that need improvement
        peaks_to_improve = [
            (peak_id, result) for peak_id, result in self.peak_results.items()
            if result.needs_improvement
        ]

        logger.info(f"🔄 Starting iterative refinement for {len(peaks_to_improve)} peaks")

        if not peaks_to_improve:
            logger.info("🎉 No peaks need improvement - all successful!")
            return {"message": "no_improvement_needed"}

        # Check if we have enough statistics for constraints
        stats_summary = self.linewidth_stats.get_statistics_summary()
        if stats_summary['peak_count'] < 3:
            logger.warning("⚠️  Insufficient successful peaks for reliable constraints")
            logger.warning("   Proceeding with fallback constraints")

        # Iterative refinement rounds
        for round_num in range(2, self.max_rounds + 2):  # Rounds 2, 3, 4, 5, 6
            logger.info(f"\n🔄 Refinement Round {round_num}")
            logger.info("-" * 40)

            round_results = self._refinement_round(round_num, peaks_to_improve, integrator_instance)

            # Check convergence
            if round_results.improvement_percentage < self.convergence_threshold * 100:
                logger.info(f"🏁 Convergence achieved: {round_results.improvement_percentage:.1f}% < {self.convergence_threshold*100:.1f}%")
                break

            # Update list of peaks that still need improvement
            peaks_to_improve = [
                (peak_id, self.peak_results[peak_id]) for peak_id, result in peaks_to_improve
                if self.peak_results[peak_id].needs_improvement
            ]

            if not peaks_to_improve:
                logger.info("🎉 All peaks now successful - refinement complete!")
                break

        return {"refinement_rounds": len(self.optimization_rounds) - 1}

    def _refinement_round(self, round_num: int, peaks_to_improve: List[Tuple[str, PeakResult]],
                         integrator_instance) -> OptimizationRound:
        """
        Execute a single refinement round with linewidth constraints
        """
        import time
        start_time = time.time()

        current_round = OptimizationRound(round_number=round_num)

        # Get current linewidth constraints
        x_constraints = self.linewidth_stats.get_x_constraints()
        y_constraints = self.linewidth_stats.get_y_constraints()

        logger.info(f"   Using X constraints: {x_constraints[0]:.4f} - {x_constraints[1]:.4f}")
        logger.info(f"   Using Y constraints: {y_constraints[0]:.1f} - {y_constraints[1]:.1f}")

        # Process peaks needing improvement
        for peak_id, old_result in peaks_to_improve:
            logger.info(f"   Refining {peak_id}...")

            try:
                # Enhanced fitting with constraints
                result = self._fit_with_constraints(
                    old_result.x_position,
                    old_result.y_position,
                    peak_id,
                    x_constraints,
                    y_constraints,
                    integrator_instance
                )

                if result and result.get('success', False):
                    # Create new peak result
                    new_peak_result = PeakResult(
                        peak_id=peak_id,
                        x_position=old_result.x_position,
                        y_position=old_result.y_position,
                        x_r_squared=result.get('x_r_squared', 0.0),
                        y_r_squared=result.get('y_r_squared', 0.0),
                        x_linewidth=result.get('x_linewidth'),
                        y_linewidth=result.get('y_linewidth'),
                        fitting_method=f"constrained_round_{round_num}",
                        optimization_round=round_num
                    )

                    # Check for improvement
                    old_min_r2 = min(old_result.x_r_squared, old_result.y_r_squared)
                    new_min_r2 = min(new_peak_result.x_r_squared, new_peak_result.y_r_squared)

                    if new_min_r2 > old_min_r2:
                        self.peak_results[peak_id] = new_peak_result
                        current_round.peaks_improved += 1

                        # Add to statistics if now successful
                        if new_peak_result.is_successful:
                            self.linewidth_stats.add_successful_peak(new_peak_result)

                        logger.info(f"     ✅ Improved: R²=({old_min_r2:.3f}→{new_min_r2:.3f})")
                    else:
                        logger.info(f"     ➖ No improvement: R²=({new_min_r2:.3f})")
                else:
                    logger.info(f"     ❌ Constrained fitting failed")

            except Exception as e:
                logger.error(f"     💥 Exception: {e}")

            current_round.peaks_processed += 1

        # Update round statistics
        for peak_result in self.peak_results.values():
            current_round.quality_distribution[peak_result.quality] += 1

        current_round.processing_time = time.time() - start_time
        self.optimization_rounds.append(current_round)

        logger.info(f"   Round {round_num}: {current_round.peaks_improved}/{current_round.peaks_processed} improved ({current_round.improvement_percentage:.1f}%)")

        return current_round

    def _fit_with_constraints(self, x_pos: float, y_pos: float, peak_id: str,
                            x_constraints: Tuple[float, float],
                            y_constraints: Tuple[float, float],
                            integrator_instance) -> Optional[Dict[str, Any]]:
        """
        Fit peak with linewidth constraints from successful peaks

        Creates linewidth constraints and applies them during fitting
        """
        try:
            # Create linewidth constraints for enhanced_voigt_fitter
            # For X-dimension (1H typically): sigma and gamma constraints
            x_linewidth_constraints = {
                'sigma_bounds': (x_constraints[0] * 0.6, x_constraints[1] * 0.6),  # Gaussian component
                'gamma_bounds': (x_constraints[0] * 0.4, x_constraints[1] * 0.4)   # Lorentzian component
            }

            # For Y-dimension (13C/15N typically): sigma and gamma constraints
            y_linewidth_constraints = {
                'sigma_bounds': (y_constraints[0] * 0.6, y_constraints[1] * 0.6),  # Gaussian component
                'gamma_bounds': (y_constraints[0] * 0.4, y_constraints[1] * 0.4)   # Lorentzian component
            }

            logger.info(f"     X constraints: σ=({x_linewidth_constraints['sigma_bounds'][0]:.4f}, {x_linewidth_constraints['sigma_bounds'][1]:.4f})")
            logger.info(f"     Y constraints: σ=({y_linewidth_constraints['sigma_bounds'][0]:.2f}, {y_linewidth_constraints['sigma_bounds'][1]:.2f})")

            # Apply constraints through enhanced_voigt_fitter directly
            # This will require updating core_integrator to pass constraints through
            result = integrator_instance.fit_peak_voigt_2d(
                x_pos, y_pos, peak_id,
                use_dynamic_optimization=True,
                all_peaks_context=[(x_pos, y_pos)],  # Provide context
                linewidth_constraints={'x': x_linewidth_constraints, 'y': y_linewidth_constraints}
            )

            return result

        except Exception as e:
            logger.error(f"Constrained fitting failed for {peak_id}: {e}")
            return None

    def _extract_linewidth_from_result(self, result: Dict[str, Any], dimension: str) -> Optional[float]:
        """
        Extract linewidth from fitting result for a specific dimension

        Args:
            result: Fitting result dictionary from core_integrator
            dimension: 'x' or 'y' dimension

        Returns:
            Total linewidth (sigma + gamma) or None if not available
        """
        try:
            # Extract fitting parameters for the dimension
            dim_key = f'{dimension}_parameters'
            if dim_key in result:
                params = result[dim_key]
                # Voigt parameters: [amplitude, center, sigma, gamma, baseline]
                if len(params) >= 4:
                    sigma = params[2]  # Gaussian width
                    gamma = params[3]  # Lorentzian width
                    # Total linewidth approximation
                    total_width = sigma + gamma
                    return total_width

            # Fallback: look for direct linewidth keys
            linewidth_key = f'{dimension}_linewidth'
            if linewidth_key in result:
                return result[linewidth_key]

            # Try alternative parameter formats
            if f'{dimension}_fit' in result:
                fit_data = result[f'{dimension}_fit']
                if isinstance(fit_data, dict) and 'parameters' in fit_data:
                    params = fit_data['parameters']
                    if len(params) >= 4:
                        return params[2] + params[3]  # sigma + gamma

            return None

        except Exception as e:
            logger.warning(f"Could not extract {dimension} linewidth from result: {e}")
            return None

    def _print_phase1_summary(self):
        """Print Phase 1 results summary"""
        round1 = self.optimization_rounds[0]
        stats = self.linewidth_stats.get_statistics_summary()

        logger.info(f"\n📊 Phase 1 Survey Complete:")
        logger.info(f"   Total peaks processed: {round1.peaks_processed}")
        logger.info(f"   Excellent (R²>0.95): {round1.quality_distribution[PeakQuality.EXCELLENT]}")
        logger.info(f"   Good (0.85-0.95): {round1.quality_distribution[PeakQuality.GOOD]}")
        logger.info(f"   Poor (0.5-0.85): {round1.quality_distribution[PeakQuality.POOR]}")
        logger.info(f"   Failed (<0.5): {round1.quality_distribution[PeakQuality.FAILED]}")
        logger.info(f"   Success rate: {round1.success_rate:.1f}%")
        logger.info(f"   Processing time: {round1.processing_time:.1f}s")

        logger.info(f"\n📏 Linewidth Statistics:")
        logger.info(f"   Successful peaks for stats: {stats['peak_count']}")
        if stats['x_linewidths']['median']:
            logger.info(f"   X-dimension median: {stats['x_linewidths']['median']:.4f}")
            logger.info(f"   Y-dimension median: {stats['y_linewidths']['median']:.2f}")

    def _print_final_summary(self):
        """Print comprehensive final summary"""
        final_round = self.optimization_rounds[-1]
        initial_round = self.optimization_rounds[0]

        logger.info(f"\n🎯 FINAL OPTIMIZATION SUMMARY")
        logger.info(f"="*50)
        logger.info(f"Total rounds: {len(self.optimization_rounds)}")
        logger.info(f"Total peaks: {self.total_peaks}")
        logger.info(f"")
        logger.info(f"Initial success rate: {initial_round.success_rate:.1f}%")
        logger.info(f"Final success rate: {final_round.success_rate:.1f}%")
        logger.info(f"Improvement: +{final_round.success_rate - initial_round.success_rate:.1f}%")
        logger.info(f"")
        logger.info(f"Final distribution:")
        logger.info(f"  Excellent: {final_round.quality_distribution[PeakQuality.EXCELLENT]}")
        logger.info(f"  Good: {final_round.quality_distribution[PeakQuality.GOOD]}")
        logger.info(f"  Poor: {final_round.quality_distribution[PeakQuality.POOR]}")
        logger.info(f"  Failed: {final_round.quality_distribution[PeakQuality.FAILED]}")

    def _generate_final_report(self) -> Dict[str, Any]:
        """Generate comprehensive final optimization report"""
        return {
            'optimization_summary': {
                'total_peaks': self.total_peaks,
                'total_rounds': len(self.optimization_rounds),
                'convergence_achieved': len(self.optimization_rounds) < self.max_rounds + 1,
                'final_success_rate': self.optimization_rounds[-1].success_rate if self.optimization_rounds else 0.0
            },
            'peak_results': {peak_id: {
                'quality': result.quality.value,
                'x_r_squared': result.x_r_squared,
                'y_r_squared': result.y_r_squared,
                'optimization_round': result.optimization_round,
                'fitting_method': result.fitting_method
            } for peak_id, result in self.peak_results.items()},
            'optimization_rounds': [
                {
                    'round': round.round_number,
                    'peaks_processed': round.peaks_processed,
                    'peaks_improved': round.peaks_improved,
                    'improvement_percentage': round.improvement_percentage,
                    'success_rate': round.success_rate,
                    'processing_time': round.processing_time
                } for round in self.optimization_rounds
            ],
            'linewidth_statistics': self.linewidth_stats.get_statistics_summary()
        }

"""
Overlap Resolver Engine for LunaNMR v0.9

This module provides the main orchestration engine for resolving overlapping NMR peaks
through a comprehensive pipeline that integrates model selection, staged fitting,
uncertainty estimation, and quality assessment.

Architecture:
    OverlapResolverEngine: Main orchestrator for complete overlap resolution
    ModelSelectionEngine: AIC/BIC-based model selection for peak count determination

Workflow:
    1. Determine optimal number of peaks (AIC/BIC model selection)
    2. Execute staged fitting with progressive parameter release
    3. Perform jackknife validation for uncertainties
    4. Analyze parameter correlations for quality assessment
    5. Generate comprehensive report with warnings and recommendations

Author: LunaNMR Development Team
Version: 0.9
"""

import numpy as np
from typing import Dict, List, Optional, Tuple, Callable
import warnings
import time
import logging
from dataclasses import dataclass

# Import the modules we created
from lunaNMR.core.multi_peak_models import (
    MultiPeakVoigtModel, VoigtModelFactory,
    multi_voigt_profile, create_initial_guess
)
from lunaNMR.core.staged_fitting_strategy import StagedFittingStrategy
from lunaNMR.core.jackknife_validator import JackknifeValidator
from lunaNMR.core.parameter_correlation_analyzer import ParameterCorrelationAnalyzer

# Configure logging
logger = logging.getLogger(__name__)


@dataclass
class ResolutionConfig:
    """Configuration for overlap resolution pipeline"""
    enable_staged_fitting: bool = True
    enable_jackknife: bool = True
    enable_correlation_analysis: bool = True
    max_peaks: int = 10
    min_peak_separation: float = 0.01
    model_selection_criterion: str = 'BIC'  # 'AIC' or 'BIC'
    min_r_squared: float = 0.7
    jackknife_resamples: int = 50
    correlation_threshold: float = 0.95
    timeout_seconds: float = 300.0


class ModelSelectionEngine:
    """
    AIC/BIC-based model selection for peak count determination

    Uses information criteria to determine the optimal number of peaks
    that balances model complexity with goodness of fit.
    """

    @staticmethod
    def calculate_aic(n_params: int, n_data: int, rss: float) -> float:
        """
        Calculate Akaike Information Criterion

        AIC = 2k + n*ln(RSS/n)

        Args:
            n_params: Number of model parameters
            n_data: Number of data points
            rss: Residual sum of squares

        Returns:
            AIC value (lower is better)
        """
        if rss <= 0 or n_data <= 0:
            return np.inf

        # Avoid log of zero/negative
        rss_per_point = rss / n_data
        if rss_per_point <= 0:
            return np.inf

        aic = 2 * n_params + n_data * np.log(rss_per_point)

        # Apply small-sample correction (AICc) if n/k < 40
        if n_data / n_params < 40:
            correction = (2 * n_params * (n_params + 1)) / (n_data - n_params - 1)
            aic += correction

        return aic

    @staticmethod
    def calculate_bic(n_params: int, n_data: int, rss: float) -> float:
        """
        Calculate Bayesian Information Criterion

        BIC = k*ln(n) + n*ln(RSS/n)

        Args:
            n_params: Number of model parameters
            n_data: Number of data points
            rss: Residual sum of squares

        Returns:
            BIC value (lower is better)
        """
        if rss <= 0 or n_data <= 0:
            return np.inf

        rss_per_point = rss / n_data
        if rss_per_point <= 0:
            return np.inf

        bic = n_params * np.log(n_data) + n_data * np.log(rss_per_point)
        return bic

    @staticmethod
    def select_best_model(aic_scores: List[float],
                         bic_scores: List[float],
                         criterion: str = 'BIC') -> int:
        """
        Select model index minimizing the chosen criterion

        Args:
            aic_scores: List of AIC values for each model
            bic_scores: List of BIC values for each model
            criterion: 'AIC' or 'BIC'

        Returns:
            Index of best model (0-based)
        """
        if criterion.upper() == 'AIC':
            scores = aic_scores
        elif criterion.upper() == 'BIC':
            scores = bic_scores
        else:
            logger.warning(f"Unknown criterion {criterion}, defaulting to BIC")
            scores = bic_scores

        # Find minimum
        valid_scores = [(i, score) for i, score in enumerate(scores) if np.isfinite(score)]

        if not valid_scores:
            logger.error("No valid model scores available")
            return 0

        best_idx = min(valid_scores, key=lambda x: x[1])[0]
        return best_idx

    @staticmethod
    def calculate_evidence_ratio(scores: List[float], best_idx: int) -> List[float]:
        """
        Calculate evidence ratios (exp(-0.5 * Δ criterion)) for all models

        Values near 1.0 indicate models with similar support from data

        Args:
            scores: Information criterion scores
            best_idx: Index of best (minimum) score

        Returns:
            List of evidence ratios relative to best model
        """
        if not scores or best_idx >= len(scores):
            return []

        best_score = scores[best_idx]
        delta_scores = [score - best_score for score in scores]
        evidence_ratios = [np.exp(-0.5 * delta) for delta in delta_scores]

        return evidence_ratios


class OverlapResolverEngine:
    """
    Main engine for resolving overlapping NMR peaks

    Orchestrates complete pipeline:
    1. Peak count determination (AIC/BIC model selection)
    2. Staged multi-peak fitting
    3. Jackknife uncertainty estimation
    4. Correlation analysis & quality assessment
    5. Comprehensive reporting with warnings and recommendations

    Example:
        >>> engine = OverlapResolverEngine()
        >>> result = engine.resolve_overlapping_peaks(
        ...     x_data=ppm_axis,
        ...     y_data=intensities,
        ...     peak_candidates=[7.5, 7.52, 7.54]
        ... )
        >>> print(f"Resolved {result['n_peaks']} peaks with R² = {result['r_squared']:.3f}")
    """

    def __init__(self, enhanced_fitter=None, config: Optional[Dict] = None):
        """
        Initialize overlap resolver engine

        Args:
            enhanced_fitter: Optional reference to enhanced voigt fitter
            config: Optional configuration dictionary to override defaults
        """
        self.enhanced_fitter = enhanced_fitter

        # Initialize component modules
        self.model_factory = MultiPeakVoigtModel  # Use the class with create_model method
        self.staged_fitter = StagedFittingStrategy()
        self.model_selector = ModelSelectionEngine()

        # Initialize configuration
        default_config = {
            'enable_staged_fitting': True,
            'enable_jackknife': True,
            'enable_correlation_analysis': True,
            'max_peaks': 10,
            'min_peak_separation': 0.01,
            'model_selection_criterion': 'BIC',
            'min_r_squared': 0.7,
            'jackknife_resamples': 50,
            'correlation_threshold': 0.95,
            'timeout_seconds': 300.0
        }

        if config:
            default_config.update(config)
        self.config = default_config

        # Initialize validators (lazy initialization for jackknife due to cost)
        self._jackknife_validator = None
        self.correlation_analyzer = ParameterCorrelationAnalyzer()

        logger.info("OverlapResolverEngine initialized with config: %s", self.config)

    @property
    def jackknife_validator(self):
        """Lazy initialization of jackknife validator"""
        if self._jackknife_validator is None:
            n_resamples = self.config.get('jackknife_resamples', 50)
            self._jackknife_validator = JackknifeValidator(n_resamples=n_resamples)
        return self._jackknife_validator

    def resolve_overlapping_peaks(self,
                                  x_data: np.ndarray,
                                  y_data: np.ndarray,
                                  peak_candidates: List[float],
                                  user_config: Optional[Dict] = None) -> Dict:
        """
        Main entry point for overlap resolution pipeline

        Args:
            x_data: PPM axis (1D array)
            y_data: Intensity data (1D array)
            peak_candidates: Initial peak position estimates (PPM values)
            user_config: Optional configuration overrides for this analysis

        Returns:
            Comprehensive result dictionary containing:
                - success: bool indicating overall success
                - method: string describing method used
                - n_peaks: optimal number of peaks determined
                - fitted_params: fitted parameters array
                - peaks: list of individual peak dictionaries
                - uncertainties: jackknife uncertainty estimates
                - correlation_analysis: correlation diagnostics
                - quality: quality metrics (R², RMSE, category)
                - warnings: list of warning messages
                - recommendations: actionable recommendations
                - model_selection: AIC/BIC scores and selection details
                - execution_time: total processing time
        """
        start_time = time.time()

        # Merge user config with defaults
        active_config = self.config.copy()
        if user_config:
            active_config.update(user_config)

        logger.info(f"Starting overlap resolution with {len(peak_candidates)} peak candidates")

        # Initialize result structure
        result = {
            'success': False,
            'method': 'overlap_resolver_engine',
            'n_peaks': 0,
            'fitted_params': None,
            'peaks': [],
            'baseline': 0.0,  # Will be updated with actual value
            'uncertainties': {},
            'correlation_analysis': {},
            'quality': {},
            'warnings': [],
            'recommendations': [],
            'model_selection': {},
            'execution_time': 0.0
        }

        # Input validation
        try:
            x_data = np.asarray(x_data)
            y_data = np.asarray(y_data)

            if x_data.shape != y_data.shape:
                raise ValueError(f"Shape mismatch: x_data {x_data.shape} vs y_data {y_data.shape}")

            if len(x_data) < 10:
                raise ValueError(f"Insufficient data points: {len(x_data)} < 10")

            if not peak_candidates:
                raise ValueError("No peak candidates provided")

        except Exception as e:
            logger.error(f"Input validation failed: {e}")
            result['warnings'].append(f"Input validation error: {str(e)}")
            result['execution_time'] = time.time() - start_time
            return result

        try:
            # Step 1: Determine optimal number of peaks
            logger.info("Step 1: Model selection to determine optimal peak count")
            model_selection_result = self._select_optimal_peak_count(
                x_data, y_data, peak_candidates, active_config
            )

            result['model_selection'] = model_selection_result

            if not model_selection_result.get('success', False):
                result['warnings'].append("Model selection failed, using all peak candidates")
                optimal_n_peaks = len(peak_candidates)
            else:
                optimal_n_peaks = model_selection_result['optimal_n_peaks']
                logger.info(f"Optimal peak count determined: {optimal_n_peaks}")

            # Use top N peak candidates
            peak_positions = sorted(peak_candidates)[:optimal_n_peaks]
            result['n_peaks'] = len(peak_positions)

            # Step 2: Execute overlap resolution fit
            logger.info(f"Step 2: Executing overlap fit for {len(peak_positions)} peaks")
            fit_result = self._execute_overlap_fit(
                x_data, y_data, peak_positions, active_config
            )

            if not fit_result.get('success', False):
                raise RuntimeError("Overlap fitting failed")

            # Store fitted parameters and results
            result['fitted_params'] = fit_result['fitted_params']
            result['peaks'] = fit_result['peaks']
            result['baseline'] = fit_result.get('baseline', 0.0)
            result['fitted_curve'] = fit_result.get('fitted_curve', None)

            # Step 3: Calculate quality metrics
            logger.info("Step 3: Calculating quality metrics")
            quality_metrics = self._calculate_quality_metrics(
                x_data, y_data, fit_result['fitted_params'], len(peak_positions)
            )
            result['quality'] = quality_metrics
            result['r_squared'] = quality_metrics['r_squared']
            result['rmse'] = quality_metrics['rmse']

            # Step 4: Jackknife validation (if enabled)
            if active_config.get('enable_jackknife', True):
                logger.info("Step 4: Jackknife uncertainty estimation")
                try:
                    uncertainties = fit_result.get('uncertainties', {})
                    result['uncertainties'] = uncertainties

                    if uncertainties.get('jackknife_failed', False):
                        result['warnings'].append("Jackknife validation failed or incomplete")

                except Exception as e:
                    logger.warning(f"Jackknife validation error: {e}")
                    result['warnings'].append(f"Jackknife validation error: {str(e)}")

            # Step 5: Correlation analysis (if enabled)
            if active_config.get('enable_correlation_analysis', True):
                logger.info("Step 5: Parameter correlation analysis")
                try:
                    correlation_result = fit_result.get('correlation_analysis', {})
                    result['correlation_analysis'] = correlation_result

                    # Extract warnings from correlation analysis
                    if 'warnings' in correlation_result:
                        result['warnings'].extend(correlation_result['warnings'])

                except Exception as e:
                    logger.warning(f"Correlation analysis error: {e}")
                    result['warnings'].append(f"Correlation analysis error: {str(e)}")

            # Step 6: Generate comprehensive report
            logger.info("Step 6: Generating comprehensive report")
            report = self._generate_comprehensive_report(
                fit_result, result['uncertainties'],
                result['correlation_analysis'], quality_metrics, active_config
            )

            result['warnings'].extend(report.get('warnings', []))
            result['recommendations'].extend(report.get('recommendations', []))

            # Final quality check
            if quality_metrics['r_squared'] < active_config.get('min_r_squared', 0.7):
                result['warnings'].append(
                    f"Low R² ({quality_metrics['r_squared']:.3f}) indicates poor fit"
                )
                result['recommendations'].append(
                    "Consider adjusting peak candidates or checking data quality"
                )

            result['success'] = True
            logger.info(f"Overlap resolution completed successfully: {optimal_n_peaks} peaks, "
                       f"R² = {quality_metrics['r_squared']:.4f}")

        except Exception as e:
            logger.error(f"Overlap resolution failed: {e}", exc_info=True)
            result['warnings'].append(f"Resolution failed: {str(e)}")
            result['success'] = False

        finally:
            result['execution_time'] = time.time() - start_time
            logger.info(f"Total execution time: {result['execution_time']:.2f} seconds")

        return result

    def _select_optimal_peak_count(self,
                                   x_data: np.ndarray,
                                   y_data: np.ndarray,
                                   peak_candidates: List[float],
                                   config: Dict) -> Dict:
        """
        Determine optimal number of peaks using AIC/BIC model selection

        Fits models with 1, 2, ..., min(len(candidates), max_peaks) peaks
        and selects the one minimizing the chosen information criterion.

        Args:
            x_data: PPM axis
            y_data: Intensity data
            peak_candidates: Initial peak position guesses
            config: Configuration dictionary

        Returns:
            Dictionary with model selection results including AIC/BIC scores
            and optimal peak count
        """
        max_peaks_to_test = min(len(peak_candidates), config.get('max_peaks', 10))
        criterion = config.get('model_selection_criterion', 'BIC')

        aic_scores = []
        bic_scores = []
        fit_results = []

        logger.info(f"Testing models with 1 to {max_peaks_to_test} peaks")

        for n_peaks in range(1, max_peaks_to_test + 1):
            try:
                # Select top n_peaks candidates (sorted by position)
                test_positions = sorted(peak_candidates)[:n_peaks]

                # Create model and initial guess
                model = self.model_factory.create_model(n_peaks)
                p0 = create_initial_guess(test_positions, y_data, x_data)

                # Quick fit without full staging for model selection
                from scipy.optimize import curve_fit

                popt, pcov = curve_fit(
                    model,
                    x_data, y_data,
                    p0=p0,
                    maxfev=1000,
                    method='trf',
                    ftol=1e-6,
                    xtol=1e-6
                )

                # Calculate residuals
                y_fit = model(x_data, *popt)
                residuals = y_data - y_fit
                rss = np.sum(residuals ** 2)

                # Calculate information criteria
                n_params = len(popt)
                n_data = len(x_data)

                aic = self.model_selector.calculate_aic(n_params, n_data, rss)
                bic = self.model_selector.calculate_bic(n_params, n_data, rss)

                aic_scores.append(aic)
                bic_scores.append(bic)
                fit_results.append({
                    'n_peaks': n_peaks,
                    'params': popt,
                    'rss': rss,
                    'success': True
                })

                logger.debug(f"  {n_peaks} peaks: AIC={aic:.2f}, BIC={bic:.2f}, RSS={rss:.2e}")

            except Exception as e:
                logger.warning(f"Failed to fit {n_peaks} peak model: {e}")
                aic_scores.append(np.inf)
                bic_scores.append(np.inf)
                fit_results.append({
                    'n_peaks': n_peaks,
                    'params': None,
                    'rss': np.inf,
                    'success': False
                })

        # Select best model
        best_idx = self.model_selector.select_best_model(aic_scores, bic_scores, criterion)
        optimal_n_peaks = best_idx + 1  # Convert from 0-based index

        # Calculate evidence ratios
        if criterion.upper() == 'AIC':
            evidence_ratios = self.model_selector.calculate_evidence_ratio(aic_scores, best_idx)
        else:
            evidence_ratios = self.model_selector.calculate_evidence_ratio(bic_scores, best_idx)

        result = {
            'success': True,
            'optimal_n_peaks': optimal_n_peaks,
            'aic_scores': aic_scores,
            'bic_scores': bic_scores,
            'evidence_ratios': evidence_ratios,
            'criterion_used': criterion,
            'best_aic': aic_scores[best_idx] if best_idx < len(aic_scores) else np.inf,
            'best_bic': bic_scores[best_idx] if best_idx < len(bic_scores) else np.inf,
            'fit_results': fit_results
        }

        logger.info(f"Model selection complete: {optimal_n_peaks} peaks optimal by {criterion}")

        return result

    def _execute_overlap_fit(self,
                            x_data: np.ndarray,
                            y_data: np.ndarray,
                            peak_positions: List[float],
                            config: Dict) -> Dict:
        """
        Execute complete overlap resolution pipeline:
        1. Staged fitting with progressive parameter release
        2. Jackknife validation for uncertainty estimation
        3. Parameter correlation analysis

        Args:
            x_data: PPM axis
            y_data: Intensity data
            peak_positions: Refined peak positions to fit
            config: Configuration dictionary

        Returns:
            Dictionary with fitted parameters, peaks, uncertainties, and diagnostics
        """
        n_peaks = len(peak_positions)
        result = {
            'success': False,
            'fitted_params': None,
            'peaks': [],
            'uncertainties': {},
            'correlation_analysis': {}
        }

        try:
            # Create model
            model = self.model_factory.create_model(n_peaks)

            # Create initial guess
            p0 = create_initial_guess(peak_positions, y_data, x_data)

            # Execute staged fitting (if enabled)
            if config.get('enable_staged_fitting', True):
                logger.info("Executing staged fitting strategy")

                staged_result = self.staged_fitter.fit_staged(
                    x_data=x_data,
                    y_data=y_data,
                    peak_centers=peak_positions,
                    model_func=model,
                    initial_params=p0,
                    verbose=False
                )

                if not staged_result.success:
                    logger.warning("Staged fitting failed, attempting direct fit")
                    raise RuntimeError("Staged fitting unsuccessful")

                fitted_params = staged_result.final_parameters
                covariance = staged_result.final_covariance

            else:
                # Direct optimization without staging
                logger.info("Executing direct optimization (staging disabled)")
                from scipy.optimize import curve_fit

                fitted_params, covariance = curve_fit(
                    model,
                    x_data, y_data,
                    p0=p0,
                    maxfev=5000,
                    method='trf'
                )

            result['fitted_params'] = fitted_params

            # Parse individual peaks
            peaks = MultiPeakVoigtModel.parse_params(
                params=fitted_params,
                n_peaks=n_peaks,
                constraint_mode='free'
            )
            result['peaks'] = peaks

            # Extract baseline (last parameter in fitted_params)
            result['baseline'] = float(fitted_params[-1])

            # Jackknife validation (if enabled and possible)
            if config.get('enable_jackknife', True):
                try:
                    logger.info("Running jackknife validation")

                    jackknife_result = self.jackknife_validator.estimate_uncertainties(
                        model_func=model,
                        x_data=x_data,
                        y_data=y_data,
                        fitted_params=fitted_params,
                        n_peaks=n_peaks
                    )

                    result['uncertainties'] = jackknife_result

                except Exception as e:
                    logger.warning(f"Jackknife validation failed: {e}")
                    result['uncertainties'] = {'jackknife_failed': True, 'error': str(e)}

            # Correlation analysis (if enabled and covariance available)
            if config.get('enable_correlation_analysis', True) and covariance is not None:
                try:
                    logger.info("Analyzing parameter correlations")

                    correlation_result = self.correlation_analyzer.analyze_correlations(
                        params=fitted_params,
                        covariance=covariance,
                        n_peaks=n_peaks
                    )

                    result['correlation_analysis'] = correlation_result

                except Exception as e:
                    logger.warning(f"Correlation analysis failed: {e}")
                    result['correlation_analysis'] = {'analysis_failed': True, 'error': str(e)}

            # Calculate fitted curve for downstream compatibility
            y_fit = model(x_data, *fitted_params)
            result['fitted_curve'] = y_fit

            result['success'] = True

        except Exception as e:
            logger.error(f"Overlap fit execution failed: {e}", exc_info=True)
            result['success'] = False
            result['error'] = str(e)

        return result

    def _calculate_quality_metrics(self,
                                   x_data: np.ndarray,
                                   y_data: np.ndarray,
                                   fitted_params: np.ndarray,
                                   n_peaks: int) -> Dict:
        """
        Calculate comprehensive quality metrics for the fit

        Args:
            x_data: PPM axis
            y_data: Intensity data
            fitted_params: Fitted parameters
            n_peaks: Number of peaks

        Returns:
            Dictionary with R², RMSE, reduced chi-squared, and quality category
        """
        try:
            # Create model and calculate fitted curve
            model = self.model_factory.create_model(n_peaks)
            y_fit = model(x_data, *fitted_params)

            # Calculate residuals
            residuals = y_data - y_fit
            rss = np.sum(residuals ** 2)

            # R-squared
            ss_tot = np.sum((y_data - np.mean(y_data)) ** 2)
            r_squared = 1.0 - (rss / ss_tot) if ss_tot > 0 else 0.0

            # RMSE
            rmse = np.sqrt(np.mean(residuals ** 2))

            # Reduced chi-squared
            n_params = len(fitted_params)
            n_data = len(x_data)
            degrees_of_freedom = max(1, n_data - n_params)
            reduced_chi_squared = rss / degrees_of_freedom

            # Quality category
            if r_squared > 0.95:
                category = 'excellent'
            elif r_squared > 0.85:
                category = 'good'
            elif r_squared > 0.70:
                category = 'fair'
            else:
                category = 'poor'

            # Normalized RMSE (relative to signal range)
            signal_range = np.ptp(y_data)
            normalized_rmse = rmse / signal_range if signal_range > 0 else np.inf

            return {
                'r_squared': r_squared,
                'rmse': rmse,
                'normalized_rmse': normalized_rmse,
                'reduced_chi_squared': reduced_chi_squared,
                'rss': rss,
                'quality_category': category,
                'degrees_of_freedom': degrees_of_freedom,
                'n_params': n_params,
                'n_data': n_data
            }

        except Exception as e:
            logger.error(f"Quality metrics calculation failed: {e}")
            return {
                'r_squared': 0.0,
                'rmse': np.inf,
                'normalized_rmse': np.inf,
                'reduced_chi_squared': np.inf,
                'quality_category': 'failed',
                'error': str(e)
            }

    def _generate_comprehensive_report(self,
                                      fit_result: Dict,
                                      uncertainties: Dict,
                                      correlations: Dict,
                                      quality: Dict,
                                      config: Dict) -> Dict:
        """
        Compile all results into final report with warnings and recommendations

        Args:
            fit_result: Results from overlap fitting
            uncertainties: Jackknife uncertainty estimates
            correlations: Parameter correlation analysis
            quality: Quality metrics
            config: Configuration used

        Returns:
            Dictionary with comprehensive warnings and recommendations
        """
        warnings_list = []
        recommendations = []

        # Quality-based warnings
        if quality['quality_category'] == 'poor':
            warnings_list.append(
                f"Poor fit quality (R² = {quality['r_squared']:.3f})"
            )
            recommendations.append(
                "Consider revising peak candidates or inspecting data for artifacts"
            )

        elif quality['quality_category'] == 'fair':
            warnings_list.append(
                f"Marginal fit quality (R² = {quality['r_squared']:.3f})"
            )
            recommendations.append(
                "Fit acceptable but may benefit from additional optimization"
            )

        # Normalized RMSE warning
        if quality.get('normalized_rmse', 0) > 0.1:
            warnings_list.append(
                f"High relative error (normalized RMSE = {quality['normalized_rmse']:.3f})"
            )

        # Correlation-based warnings
        if correlations.get('has_high_correlations', False):
            warnings_list.append("High parameter correlations detected")
            recommendations.append(
                "Parameters may be poorly determined due to overlap. "
                "Consider constraining some parameters or using prior knowledge."
            )

        if 'problematic_pairs' in correlations:
            n_problematic = len(correlations['problematic_pairs'])
            if n_problematic > 0:
                warnings_list.append(
                    f"{n_problematic} parameter pair(s) highly correlated (>0.95)"
                )

        # Uncertainty-based warnings
        if uncertainties.get('jackknife_failed', False):
            warnings_list.append("Uncertainty estimation incomplete")
            recommendations.append(
                "Consider manual inspection of fit stability"
            )

        if uncertainties.get('high_uncertainty_peaks', []):
            high_unc_count = len(uncertainties['high_uncertainty_peaks'])
            warnings_list.append(
                f"{high_unc_count} peak(s) with high parameter uncertainties"
            )
            recommendations.append(
                "Peaks with high uncertainty may be poorly resolved or artifacts"
            )

        # Stability warnings from jackknife
        if uncertainties.get('unstable_peaks', []):
            unstable_count = len(uncertainties['unstable_peaks'])
            warnings_list.append(
                f"{unstable_count} peak(s) show fitting instability"
            )
            recommendations.append(
                "Unstable peaks may indicate overfitting or poor initialization"
            )

        # Model selection warnings
        if fit_result.get('staged_stages_completed', 0) < 4:
            warnings_list.append(
                "Staged fitting did not complete all stages"
            )

        # Generate summary recommendation
        if quality['quality_category'] in ['excellent', 'good'] and \
           not correlations.get('has_high_correlations', False):
            recommendations.append(
                "Fit appears reliable. Results can be used with confidence."
            )
        elif quality['quality_category'] in ['excellent', 'good']:
            recommendations.append(
                "Fit quality good but high correlations suggest caution in interpretation."
            )

        return {
            'warnings': warnings_list,
            'recommendations': recommendations,
            'overall_assessment': quality['quality_category']
        }

    def get_configuration(self) -> Dict:
        """Get current configuration"""
        return self.config.copy()

    def update_configuration(self, new_config: Dict) -> None:
        """Update configuration with new values"""
        self.config.update(new_config)
        logger.info(f"Configuration updated: {new_config}")

    def reset_configuration(self) -> None:
        """Reset to default configuration"""
        self.config = {
            'enable_staged_fitting': True,
            'enable_jackknife': True,
            'enable_correlation_analysis': True,
            'max_peaks': 10,
            'min_peak_separation': 0.01,
            'model_selection_criterion': 'BIC',
            'min_r_squared': 0.7,
            'jackknife_resamples': 50,
            'correlation_threshold': 0.95,
            'timeout_seconds': 300.0
        }
        logger.info("Configuration reset to defaults")


# Convenience function for quick overlap resolution
def resolve_overlap(x_data: np.ndarray,
                   y_data: np.ndarray,
                   peak_candidates: List[float],
                   config: Optional[Dict] = None) -> Dict:
    """
    Convenience function for quick overlap resolution

    Args:
        x_data: PPM axis
        y_data: Intensity data
        peak_candidates: Initial peak positions
        config: Optional configuration

    Returns:
        Complete resolution result dictionary
    """
    engine = OverlapResolverEngine(config=config)
    return engine.resolve_overlapping_peaks(x_data, y_data, peak_candidates)


if __name__ == "__main__":
    # Example usage and testing
    import matplotlib.pyplot as plt

    # Create synthetic overlapping peaks
    x = np.linspace(7.0, 8.0, 1000)

    # True parameters for 3 overlapping peaks
    true_params = [
        100.0, 7.3, 0.02, 0.01,  # Peak 1: amplitude, center, sigma, gamma
        80.0, 7.35, 0.025, 0.015,  # Peak 2
        60.0, 7.4, 0.02, 0.01,   # Peak 3
        5.0  # baseline
    ]

    # Generate synthetic data
    y_true = multi_voigt_profile(x, *true_params)
    noise = np.random.normal(0, 2, len(x))
    y_data = y_true + noise

    # Peak candidates (slightly perturbed from true positions)
    peak_candidates = [7.305, 7.355, 7.395]

    # Initialize and run engine
    print("Initializing Overlap Resolver Engine...")
    engine = OverlapResolverEngine()

    print("\nResolving overlapping peaks...")
    result = engine.resolve_overlapping_peaks(x, y_data, peak_candidates)

    # Print results
    print(f"\n{'='*60}")
    print("OVERLAP RESOLUTION RESULTS")
    print(f"{'='*60}")
    print(f"Success: {result['success']}")
    print(f"Optimal peaks: {result['n_peaks']}")
    print(f"R²: {result.get('r_squared', 0):.4f}")
    print(f"RMSE: {result.get('rmse', 0):.4f}")
    print(f"Quality: {result.get('quality', {}).get('quality_category', 'unknown')}")
    print(f"Execution time: {result['execution_time']:.2f} s")

    if result['model_selection']:
        print(f"\nModel Selection ({result['model_selection']['criterion_used']}):")
        for i, (aic, bic) in enumerate(zip(
            result['model_selection']['aic_scores'],
            result['model_selection']['bic_scores']
        ), 1):
            marker = " <-- SELECTED" if i == result['n_peaks'] else ""
            print(f"  {i} peaks: AIC={aic:.2f}, BIC={bic:.2f}{marker}")

    print(f"\nFitted Peaks:")
    for i, peak in enumerate(result.get('peaks', []), 1):
        print(f"  Peak {i}: center={peak['center']:.4f}, "
              f"amplitude={peak['amplitude']:.2f}, "
              f"sigma={peak['sigma']:.4f}, "
              f"gamma={peak['gamma']:.4f}")

    if result.get('warnings'):
        print(f"\nWarnings:")
        for warning in result['warnings']:
            print(f"  - {warning}")

    if result.get('recommendations'):
        print(f"\nRecommendations:")
        for rec in result['recommendations']:
            print(f"  - {rec}")

    # Plot results
    if result['success'] and result['fitted_params'] is not None:
        plt.figure(figsize=(12, 6))

        # Calculate fitted curve
        model = MultiPeakVoigtModel.create_model(result['n_peaks'])
        y_fit = model(x, *result['fitted_params'])

        plt.plot(x, y_data, 'k.', alpha=0.3, label='Data')
        plt.plot(x, y_fit, 'r-', linewidth=2, label='Fitted')
        plt.plot(x, y_true, 'g--', linewidth=1, label='True')

        plt.xlabel('Chemical Shift (ppm)')
        plt.ylabel('Intensity')
        plt.title(f'Overlap Resolution: {result["n_peaks"]} peaks, R²={result["r_squared"]:.4f}')
        plt.legend()
        plt.grid(True, alpha=0.3)

        plt.tight_layout()
        plt.savefig('/tmp/overlap_resolution_test.png', dpi=150)
        print(f"\nPlot saved to /tmp/overlap_resolution_test.png")

    print(f"\n{'='*60}\n")

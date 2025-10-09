"""
Parameter Correlation Analyzer for LunaNMR

Analyzes covariance matrices from Voigt fitting to detect unresolvable peak overlaps
and assess the reliability of multi-peak deconvolution.

Inspired by  covariance analysis approach.
"""

import numpy as np
from typing import Dict, List, Tuple, Optional
import warnings


class ParameterCorrelationAnalyzer:
    """
    Analyze parameter correlations to assess overlap resolution quality

    High amplitude correlation between peaks (r > 0.7) indicates:
    - Peaks are too close to resolve reliably
    - Consider merging or using stronger constraints

    High σ-γ correlation within a peak (r > 0.85) indicates:
    - Voigt profile may be overparameterized
    - Consider simpler profile (Gaussian or Lorentzian)
    """

    def __init__(self,
                 amp_corr_threshold: float = 0.7,
                 voigt_sg_corr_threshold: float = 0.85):
        """
        Initialize correlation analyzer

        Parameters:
        -----------
        amp_corr_threshold : float
            Threshold for inter-peak amplitude correlations (default: 0.7)
        voigt_sg_corr_threshold : float
            Threshold for intra-peak σ-γ correlations (default: 0.85)
        """
        self.amp_correlation_threshold = amp_corr_threshold
        self.voigt_sg_correlation_threshold = voigt_sg_corr_threshold

    def analyze_correlations(self,
                            params: np.ndarray,
                            covariance: np.ndarray,
                            n_peaks: int) -> Dict:
        """
        Comprehensive correlation analysis

        Parameters:
        -----------
        params : np.ndarray
            Fitted parameters [amp1, center1, sigma1, gamma1, ..., baseline]
        covariance : np.ndarray
            Covariance matrix from fit
        n_peaks : int
            Number of peaks in the fit

        Returns:
        --------
        Dict containing:
        {
            'correlation_matrix': np.ndarray,
            'inter_peak_correlations': List[Dict],
            'voigt_sigma_gamma_correlations': List[float],
            'overlap_warnings': List[Dict],
            'resolvability_score': float,  # 0-1 scale
            'resolvability_class': str     # Classification string
        }
        """
        # Compute correlation matrix
        correlation = self._compute_correlation_matrix(covariance)

        # Analyze inter-peak correlations
        inter_peak_corr = self._analyze_inter_peak_correlations(
            correlation, n_peaks
        )

        # Analyze Voigt σ-γ correlations
        voigt_sg_corr = self._analyze_voigt_internal_correlations(
            correlation, n_peaks
        )

        # Identify problematic correlations
        overlap_warnings = self._identify_overlap_warnings(
            inter_peak_corr, voigt_sg_corr, n_peaks
        )

        # Compute overall resolvability score
        resolvability_score = self._compute_resolvability_score(
            inter_peak_corr, voigt_sg_corr
        )

        # Build results dictionary
        results = {
            'correlation_matrix': correlation,
            'inter_peak_correlations': inter_peak_corr,
            'voigt_sigma_gamma_correlations': voigt_sg_corr,
            'overlap_warnings': overlap_warnings,
            'resolvability_score': resolvability_score,
            'resolvability_class': self.assess_resolvability(inter_peak_corr, voigt_sg_corr),
            'n_peaks': n_peaks
        }

        return results

    def _compute_correlation_matrix(self,
                                   covariance: np.ndarray) -> np.ndarray:
        """
        Convert covariance to correlation matrix

        correlation[i,j] = covariance[i,j] / sqrt(var[i] * var[j])

        Parameters:
        -----------
        covariance : np.ndarray
            Covariance matrix

        Returns:
        --------
        np.ndarray : Correlation matrix
        """
        # Extract variances (diagonal elements)
        variances = np.diag(covariance)

        # Handle zero or negative variances
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            std_devs = np.sqrt(np.abs(variances))

        # Create correlation matrix
        correlation = np.zeros_like(covariance)

        for i in range(len(covariance)):
            for j in range(len(covariance)):
                if std_devs[i] > 0 and std_devs[j] > 0:
                    correlation[i, j] = covariance[i, j] / (std_devs[i] * std_devs[j])
                else:
                    # Set to 0 if either variance is zero
                    correlation[i, j] = 0.0

        # Ensure diagonal is exactly 1.0 where possible
        for i in range(len(correlation)):
            if std_devs[i] > 0:
                correlation[i, i] = 1.0

        return correlation

    def _analyze_inter_peak_correlations(self,
                                        correlation: np.ndarray,
                                        n_peaks: int) -> List[Dict]:
        """
        Analyze correlations between different peaks

        Parameter order: [amp, center, sigma, gamma] per peak
        Focus on amplitude-amplitude and position-position correlations

        Parameters:
        -----------
        correlation : np.ndarray
            Correlation matrix
        n_peaks : int
            Number of peaks

        Returns:
        --------
        List[Dict] : List of correlation information dictionaries
        """
        inter_peak_correlations = []

        # Iterate over all peak pairs
        for i in range(n_peaks):
            for j in range(i + 1, n_peaks):
                # Parameter indices for peak i and peak j
                # [amp, center, sigma, gamma] -> indices [0, 1, 2, 3] offset by 4*peak_index
                amp_i = 4 * i
                center_i = 4 * i + 1
                sigma_i = 4 * i + 2
                gamma_i = 4 * i + 3

                amp_j = 4 * j
                center_j = 4 * j + 1
                sigma_j = 4 * j + 2
                gamma_j = 4 * j + 3

                # Extract key correlations
                amp_amp_corr = correlation[amp_i, amp_j]
                center_center_corr = correlation[center_i, center_j]

                # Also check cross-correlations that might indicate issues
                amp_i_center_j = correlation[amp_i, center_j]
                amp_j_center_i = correlation[amp_j, center_i]

                corr_info = {
                    'peak_pair': (i, j),
                    'amplitude_correlation': amp_amp_corr,
                    'position_correlation': center_center_corr,
                    'cross_correlations': {
                        'amp_i_center_j': amp_i_center_j,
                        'amp_j_center_i': amp_j_center_i
                    },
                    'max_correlation': max(abs(amp_amp_corr),
                                          abs(center_center_corr),
                                          abs(amp_i_center_j),
                                          abs(amp_j_center_i))
                }

                inter_peak_correlations.append(corr_info)

        return inter_peak_correlations

    def _analyze_voigt_internal_correlations(self,
                                            correlation: np.ndarray,
                                            n_peaks: int) -> List[float]:
        """
        Analyze sigma-gamma correlation within each peak

        Parameters:
        -----------
        correlation : np.ndarray
            Correlation matrix
        n_peaks : int
            Number of peaks

        Returns:
        --------
        List[float] : List of σ-γ correlations (one per peak)
        """
        sg_correlations = []

        for i in range(n_peaks):
            # Indices for sigma and gamma of peak i
            sigma_idx = 4 * i + 2
            gamma_idx = 4 * i + 3

            # Extract σ-γ correlation
            sg_corr = correlation[sigma_idx, gamma_idx]
            sg_correlations.append(sg_corr)

        return sg_correlations

    def _identify_overlap_warnings(self,
                                  inter_peak_corr: List[Dict],
                                  voigt_sg_corr: List[float],
                                  n_peaks: int) -> List[Dict]:
        """
        Identify problematic correlations and generate warnings

        Parameters:
        -----------
        inter_peak_corr : List[Dict]
            Inter-peak correlation data
        voigt_sg_corr : List[float]
            Voigt σ-γ correlations
        n_peaks : int
            Number of peaks

        Returns:
        --------
        List[Dict] : Warning messages with details
        """
        warnings_list = []

        # Check inter-peak correlations
        for corr_info in inter_peak_corr:
            peak_i, peak_j = corr_info['peak_pair']
            max_corr = corr_info['max_correlation']

            if abs(max_corr) >= self.amp_correlation_threshold:
                warning = {
                    'type': 'inter_peak',
                    'peaks': (peak_i, peak_j),
                    'correlation': max_corr,
                    'severity': self._get_severity(abs(max_corr),
                                                  self.amp_correlation_threshold),
                    'message': f"High correlation between peaks {peak_i} and {peak_j} (r={max_corr:.3f})"
                }
                warnings_list.append(warning)

        # Check Voigt σ-γ correlations
        for i, sg_corr in enumerate(voigt_sg_corr):
            if abs(sg_corr) >= self.voigt_sg_correlation_threshold:
                warning = {
                    'type': 'voigt_sg',
                    'peak': i,
                    'correlation': sg_corr,
                    'severity': self._get_severity(abs(sg_corr),
                                                  self.voigt_sg_correlation_threshold),
                    'message': f"High σ-γ correlation for peak {i} (r={sg_corr:.3f})"
                }
                warnings_list.append(warning)

        return warnings_list

    def _get_severity(self, correlation: float, threshold: float) -> str:
        """
        Classify correlation severity

        Parameters:
        -----------
        correlation : float
            Absolute correlation value
        threshold : float
            Warning threshold

        Returns:
        --------
        str : 'critical', 'high', or 'moderate'
        """
        if correlation >= 0.9:
            return 'critical'
        elif correlation >= threshold + 0.1:
            return 'high'
        else:
            return 'moderate'

    def _compute_resolvability_score(self,
                                    inter_peak_corr: List[Dict],
                                    voigt_sg_corr: List[float]) -> float:
        """
        Compute overall resolvability score (0-1 scale, 1=best)

        Parameters:
        -----------
        inter_peak_corr : List[Dict]
            Inter-peak correlations
        voigt_sg_corr : List[float]
            Voigt σ-γ correlations

        Returns:
        --------
        float : Resolvability score (0-1)
        """
        if not inter_peak_corr and not voigt_sg_corr:
            return 1.0

        # Find maximum inter-peak correlation
        max_inter_peak = 0.0
        if inter_peak_corr:
            max_inter_peak = max([abs(c['max_correlation']) for c in inter_peak_corr])

        # Find maximum σ-γ correlation
        max_sg = 0.0
        if voigt_sg_corr:
            max_sg = max([abs(c) for c in voigt_sg_corr])

        # Score based on highest correlation (1.0 = no correlation, 0.0 = perfect correlation)
        # Weight inter-peak correlations more heavily
        score = 1.0 - (0.7 * max_inter_peak + 0.3 * max_sg)

        return max(0.0, min(1.0, score))

    def assess_resolvability(self,
                           inter_peak_corr: List[Dict],
                           voigt_sg_corr: List[float]) -> str:
        """
        Overall resolvability classification

        Parameters:
        -----------
        inter_peak_corr : List[Dict]
            Inter-peak correlations
        voigt_sg_corr : List[float]
            Voigt σ-γ correlations

        Returns:
        --------
        str : 'well_resolved', 'marginal', 'poorly_resolved', 'unresolvable'

        Logic:
        - Any inter-peak correlation > 0.9: 'unresolvable'
        - Multiple correlations > 0.8: 'poorly_resolved'
        - Some correlations > 0.7: 'marginal'
        - All correlations < 0.7: 'well_resolved'
        """
        # Check for unresolvable correlations
        for corr_info in inter_peak_corr:
            if abs(corr_info['max_correlation']) > 0.9:
                return 'unresolvable'

        # Count high correlations
        high_corr_count = sum(1 for c in inter_peak_corr
                             if abs(c['max_correlation']) > 0.8)

        if high_corr_count >= 2:
            return 'poorly_resolved'

        # Check for marginal correlations
        marginal_corr_count = sum(1 for c in inter_peak_corr
                                 if abs(c['max_correlation']) > 0.7)

        if marginal_corr_count > 0:
            return 'marginal'

        return 'well_resolved'

    def generate_recommendations(self,
                                correlation_analysis: Dict) -> List[str]:
        """
        Generate actionable recommendations based on correlation analysis

        Parameters:
        -----------
        correlation_analysis : Dict
            Results from analyze_correlations()

        Returns:
        --------
        List[str] : Actionable recommendations

        Examples:
        - "Peaks 0-1 highly correlated (r=0.92): consider merging"
        - "Peak 2 σ-γ correlation high (r=0.88): try Gaussian profile"
        - "All peaks well-resolved: fit is reliable"
        """
        recommendations = []

        warnings_list = correlation_analysis['overlap_warnings']
        resolvability_class = correlation_analysis['resolvability_class']
        resolvability_score = correlation_analysis['resolvability_score']

        # Overall assessment
        if resolvability_class == 'well_resolved':
            recommendations.append(
                f"✓ All peaks well-resolved (score: {resolvability_score:.2f}): fit is reliable"
            )
        elif resolvability_class == 'marginal':
            recommendations.append(
                f"⚠ Marginal resolution (score: {resolvability_score:.2f}): interpret with caution"
            )
        elif resolvability_class == 'poorly_resolved':
            recommendations.append(
                f"⚠ Poor resolution (score: {resolvability_score:.2f}): results may be unreliable"
            )
        else:  # unresolvable
            recommendations.append(
                f"✗ Unresolvable overlap (score: {resolvability_score:.2f}): fit parameters not trustworthy"
            )

        # Specific warnings and recommendations
        for warning in warnings_list:
            if warning['type'] == 'inter_peak':
                peak_i, peak_j = warning['peaks']
                corr = warning['correlation']

                if warning['severity'] == 'critical':
                    recommendations.append(
                        f"• Peaks {peak_i}-{peak_j} critically correlated (r={corr:.3f}): "
                        f"consider merging into single peak"
                    )
                elif warning['severity'] == 'high':
                    recommendations.append(
                        f"• Peaks {peak_i}-{peak_j} highly correlated (r={corr:.3f}): "
                        f"consider constraining peak positions or using stronger regularization"
                    )
                else:
                    recommendations.append(
                        f"• Peaks {peak_i}-{peak_j} moderately correlated (r={corr:.3f}): "
                        f"monitor fit quality"
                    )

            elif warning['type'] == 'voigt_sg':
                peak = warning['peak']
                corr = warning['correlation']

                if warning['severity'] == 'critical':
                    recommendations.append(
                        f"• Peak {peak} σ-γ critically correlated (r={corr:.3f}): "
                        f"Voigt profile overparameterized, try pure Gaussian or Lorentzian"
                    )
                elif warning['severity'] == 'high':
                    recommendations.append(
                        f"• Peak {peak} σ-γ highly correlated (r={corr:.3f}): "
                        f"consider simpler lineshape model"
                    )
                else:
                    recommendations.append(
                        f"• Peak {peak} σ-γ moderately correlated (r={corr:.3f}): "
                        f"Voigt profile may be adequate but monitor"
                    )

        return recommendations


class CorrelationVisualizer:
    """Visualization and reporting utilities for correlation analysis"""

    @staticmethod
    def format_correlation_report(correlation_analysis: Dict) -> str:
        """
        Format correlation analysis as readable text report

        Parameters:
        -----------
        correlation_analysis : Dict
            Results from ParameterCorrelationAnalyzer.analyze_correlations()

        Returns:
        --------
        str : Formatted text report
        """
        report_lines = []

        # Header
        report_lines.append("=" * 70)
        report_lines.append("PARAMETER CORRELATION ANALYSIS REPORT")
        report_lines.append("=" * 70)
        report_lines.append("")

        # Overall assessment
        n_peaks = correlation_analysis['n_peaks']
        resolvability_class = correlation_analysis['resolvability_class']
        resolvability_score = correlation_analysis['resolvability_score']

        report_lines.append(f"Number of peaks: {n_peaks}")
        report_lines.append(f"Resolvability class: {resolvability_class.upper()}")
        report_lines.append(f"Resolvability score: {resolvability_score:.3f} (0=poor, 1=excellent)")
        report_lines.append("")

        # Inter-peak correlations
        report_lines.append("-" * 70)
        report_lines.append("INTER-PEAK CORRELATIONS")
        report_lines.append("-" * 70)

        inter_peak_corr = correlation_analysis['inter_peak_correlations']
        if inter_peak_corr:
            for corr_info in inter_peak_corr:
                peak_i, peak_j = corr_info['peak_pair']
                amp_corr = corr_info['amplitude_correlation']
                pos_corr = corr_info['position_correlation']
                max_corr = corr_info['max_correlation']

                report_lines.append(f"  Peaks {peak_i} - {peak_j}:")
                report_lines.append(f"    Amplitude correlation:  {amp_corr:+.3f}")
                report_lines.append(f"    Position correlation:   {pos_corr:+.3f}")
                report_lines.append(f"    Maximum correlation:    {max_corr:+.3f}")
                report_lines.append("")
        else:
            report_lines.append("  Single peak - no inter-peak correlations")
            report_lines.append("")

        # Voigt σ-γ correlations
        report_lines.append("-" * 70)
        report_lines.append("VOIGT σ-γ CORRELATIONS (within each peak)")
        report_lines.append("-" * 70)

        voigt_sg_corr = correlation_analysis['voigt_sigma_gamma_correlations']
        for i, sg_corr in enumerate(voigt_sg_corr):
            status = "OK" if abs(sg_corr) < 0.85 else "HIGH"
            report_lines.append(f"  Peak {i}: r(σ,γ) = {sg_corr:+.3f} [{status}]")
        report_lines.append("")

        # Warnings
        warnings_list = correlation_analysis['overlap_warnings']
        if warnings_list:
            report_lines.append("-" * 70)
            report_lines.append("WARNINGS")
            report_lines.append("-" * 70)
            for warning in warnings_list:
                severity_marker = {
                    'critical': '✗✗✗',
                    'high': '⚠⚠',
                    'moderate': '⚠'
                }
                marker = severity_marker.get(warning['severity'], '•')
                report_lines.append(f"  {marker} {warning['message']}")
            report_lines.append("")

        # Recommendations
        analyzer = ParameterCorrelationAnalyzer()
        recommendations = analyzer.generate_recommendations(correlation_analysis)

        report_lines.append("-" * 70)
        report_lines.append("RECOMMENDATIONS")
        report_lines.append("-" * 70)
        for rec in recommendations:
            report_lines.append(f"  {rec}")

        report_lines.append("")
        report_lines.append("=" * 70)

        return "\n".join(report_lines)

    @staticmethod
    def format_correlation_matrix_summary(correlation_matrix: np.ndarray,
                                         n_peaks: int,
                                         threshold: float = 0.5) -> str:
        """
        Format correlation matrix with parameter labels

        Parameters:
        -----------
        correlation_matrix : np.ndarray
            Correlation matrix
        n_peaks : int
            Number of peaks
        threshold : float
            Only show correlations above this threshold (default: 0.5)

        Returns:
        --------
        str : Formatted matrix summary
        """
        param_labels = []
        for i in range(n_peaks):
            param_labels.extend([
                f'A{i}',      # amplitude
                f'μ{i}',      # center
                f'σ{i}',      # sigma
                f'γ{i}'       # gamma
            ])
        param_labels.append('B')  # baseline

        lines = []
        lines.append("HIGH CORRELATIONS (|r| > {:.2f}):".format(threshold))
        lines.append("")

        # Find high correlations (excluding diagonal)
        high_corr = []
        n = len(correlation_matrix)
        for i in range(n):
            for j in range(i + 1, n):
                if abs(correlation_matrix[i, j]) > threshold:
                    high_corr.append((i, j, correlation_matrix[i, j]))

        # Sort by absolute correlation value
        high_corr.sort(key=lambda x: abs(x[2]), reverse=True)

        if high_corr:
            for i, j, corr in high_corr:
                lines.append(f"  {param_labels[i]:>3s} - {param_labels[j]:>3s}: {corr:+.3f}")
        else:
            lines.append(f"  No correlations above threshold ({threshold:.2f})")

        return "\n".join(lines)

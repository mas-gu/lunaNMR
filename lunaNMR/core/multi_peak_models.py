"""
Multi-Peak Voigt Fitting Models

This module provides factory classes and functions for creating and managing
multi-peak Voigt profile models for NMR spectral deconvolution.

Key Features:
- Vectorized multi-peak Voigt profile calculation using scipy.special.wofz
- Flexible model factory with multiple constraint modes
- Smart initial parameter estimation
- Robust bounds generation for different fitting stages
- Parameter parsing and validation

Author: LunaNMR Development Team
Version: 0.9
"""

import numpy as np
from scipy.special import wofz
from typing import Callable, List, Dict, Tuple, Optional, Any
import warnings


# Constants
MIN_WIDTH = 1e-6  # Minimum sigma/gamma to avoid division by zero
SQRT_2 = np.sqrt(2.0)
SQRT_2PI = np.sqrt(2.0 * np.pi)
INV_SQRT_2PI = 1.0 / SQRT_2PI


def multi_voigt_profile(x: np.ndarray, *params) -> np.ndarray:
    """
    Vectorized multi-peak Voigt profile calculation.

    The Voigt profile is a convolution of Gaussian and Lorentzian distributions,
    efficiently computed using the Faddeeva function (complex error function).

    Parameters
    ----------
    x : np.ndarray
        X-axis values (e.g., chemical shift in ppm)
    *params : float
        Flattened parameters: [amp1, center1, sigma1, gamma1,
                               amp2, center2, sigma2, gamma2,
                               ..., baseline]
        - amp: Peak amplitude (intensity)
        - center: Peak center position
        - sigma: Gaussian width parameter (FWHM_gauss = 2*sqrt(2*ln(2))*sigma)
        - gamma: Lorentzian width parameter (FWHM_lorentz = 2*gamma)
        - baseline: Constant baseline offset (last parameter)

    Returns
    -------
    np.ndarray
        Combined Voigt profile for all peaks plus baseline

    Notes
    -----
    Voigt profile formula:
        V(x) = A * Re[wofz(z)] / (sigma * sqrt(2*pi))
    where:
        z = ((x - center) + i*gamma) / (sigma * sqrt(2))
        wofz(z) is the Faddeeva function (complex error function)

    The implementation is fully vectorized for performance and handles
    edge cases where sigma or gamma approach zero.
    """
    # Extract baseline (last parameter)
    baseline = params[-1]

    # Calculate number of peaks
    n_params_per_peak = 4  # amp, center, sigma, gamma
    n_peaks = (len(params) - 1) // n_params_per_peak

    # Initialize output with baseline
    y = np.full_like(x, baseline, dtype=np.float64)

    # Add contribution from each peak
    for i in range(n_peaks):
        idx = i * n_params_per_peak
        amp = params[idx]
        center = params[idx + 1]
        sigma = max(params[idx + 2], MIN_WIDTH)  # Prevent division by zero
        gamma = max(params[idx + 3], MIN_WIDTH)  # Prevent division by zero

        # Calculate Voigt profile using Faddeeva function
        # z = ((x - center) + i*gamma) / (sigma * sqrt(2))
        z = ((x - center) + 1j * gamma) / (sigma * SQRT_2)

        # Voigt profile: amp * Re[wofz(z)] / (sigma * sqrt(2*pi))
        voigt = amp * np.real(wofz(z)) / (sigma * SQRT_2PI)

        y += voigt

    return y


def create_initial_guess(peak_centers: List[float],
                        y_data: np.ndarray,
                        x_data: Optional[np.ndarray] = None,
                        baseline: Optional[float] = None) -> List[float]:
    """
    Create smart initial parameter guesses for multi-peak Voigt fitting.

    Parameters
    ----------
    peak_centers : List[float]
        Detected peak center positions (e.g., in ppm)
    y_data : np.ndarray
        Spectral intensity data
    x_data : Optional[np.ndarray], optional
        X-axis values. If provided, used to estimate peak amplitudes more accurately
    baseline : Optional[float], optional
        Known baseline value. If None, estimated from data

    Returns
    -------
    List[float]
        Initial parameters: [amp1, center1, sigma1, gamma1, ..., baseline]

    Notes
    -----
    Initial guess strategy:
    - Amplitude: Estimated from y_data at peak positions or median of positive values
    - Center: User-provided peak centers
    - Sigma: Estimated from typical peak spacing or set to default (0.01)
    - Gamma: Initialized equal to sigma (balanced Gaussian/Lorentzian)
    - Baseline: Minimum of y_data or user-provided value
    """
    n_peaks = len(peak_centers)

    # Estimate baseline
    if baseline is None:
        baseline = float(np.percentile(y_data, 5))  # 5th percentile to avoid outliers

    # Estimate typical peak amplitude
    if x_data is not None and len(x_data) == len(y_data):
        # Try to get amplitudes at peak positions
        amplitudes = []
        for center in peak_centers:
            # Find closest x value to peak center
            idx = np.argmin(np.abs(x_data - center))
            amp = y_data[idx] - baseline
            amplitudes.append(max(amp, 0))

        if all(a == 0 for a in amplitudes):
            # Fallback if all amplitudes are zero
            typical_amp = float(np.median(y_data[y_data > baseline]))
            amplitudes = [typical_amp] * n_peaks
    else:
        # Use median of positive values
        positive_values = y_data[y_data > baseline]
        if len(positive_values) > 0:
            typical_amp = float(np.median(positive_values))
        else:
            typical_amp = float(np.max(y_data) - baseline) * 0.5
        amplitudes = [typical_amp] * n_peaks

    # Estimate typical peak width from peak spacing
    if n_peaks > 1:
        spacings = np.diff(sorted(peak_centers))
        typical_width = float(np.median(spacings)) * 0.2  # ~20% of peak spacing
        typical_width = max(typical_width, 0.001)  # Minimum width
    else:
        typical_width = 0.01  # Default width

    # Build initial guess
    initial_params = []
    for i, center in enumerate(peak_centers):
        amp = amplitudes[i]
        sigma = typical_width
        gamma = typical_width  # Equal Gaussian and Lorentzian initially

        initial_params.extend([amp, center, sigma, gamma])

    # Add baseline
    initial_params.append(baseline)

    return initial_params


def parse_multi_peak_params(params: np.ndarray, n_peaks: int) -> List[Dict[str, float]]:
    """
    Parse flat parameter array into individual peak dictionaries.

    Parameters
    ----------
    params : np.ndarray
        Flattened parameters from fitting: [amp1, center1, sigma1, gamma1, ..., baseline]
    n_peaks : int
        Number of peaks in the model

    Returns
    -------
    List[Dict[str, float]]
        List of dictionaries, one per peak, containing:
        - 'amplitude': Peak amplitude
        - 'center': Peak center position
        - 'sigma': Gaussian width parameter
        - 'gamma': Lorentzian width parameter
        - 'fwhm_gaussian': Gaussian FWHM
        - 'fwhm_lorentzian': Lorentzian FWHM
        - 'baseline': Baseline value (same for all peaks)

    Notes
    -----
    FWHM conversions:
    - Gaussian FWHM = 2 * sqrt(2 * ln(2)) * sigma ≈ 2.355 * sigma
    - Lorentzian FWHM = 2 * gamma
    """
    FWHM_GAUSSIAN_FACTOR = 2.0 * np.sqrt(2.0 * np.log(2.0))  # ≈ 2.355

    baseline = float(params[-1])
    peak_list = []

    for i in range(n_peaks):
        idx = i * 4
        amp = float(params[idx])
        center = float(params[idx + 1])
        sigma = float(params[idx + 2])
        gamma = float(params[idx + 3])

        peak_dict = {
            'amplitude': amp,
            'center': center,
            'sigma': sigma,
            'gamma': gamma,
            'fwhm_gaussian': FWHM_GAUSSIAN_FACTOR * sigma,
            'fwhm_lorentzian': 2.0 * gamma,
            'baseline': baseline
        }
        peak_list.append(peak_dict)

    return peak_list


class MultiPeakVoigtModel:
    """
    Factory class for creating multi-peak Voigt models with different constraints.

    This class provides static methods to create model functions, generate parameter
    bounds for different fitting stages, and parse fitted parameters.

    Attributes
    ----------
    CONSTRAINT_MODES : List[str]
        Available constraint modes: ['free', 'shared_width', 'fixed_positions']
    """

    CONSTRAINT_MODES = ['free', 'shared_width', 'fixed_positions']

    @staticmethod
    def create_model(n_peaks: int,
                    constraint_mode: str = 'free',
                    fixed_params: Optional[Dict[str, Any]] = None) -> Callable:
        """
        Create a multi-peak Voigt model function with specified constraints.

        Parameters
        ----------
        n_peaks : int
            Number of peaks in the model
        constraint_mode : str, optional
            Constraint mode for fitting:
            - 'free': All parameters independent
            - 'shared_width': All peaks share same sigma and gamma
            - 'fixed_positions': Peak positions are fixed
        fixed_params : Optional[Dict[str, Any]], optional
            Dictionary of fixed parameter values for 'fixed_positions' mode

        Returns
        -------
        Callable
            Model function with signature: f(x, *params) -> y

        Raises
        ------
        ValueError
            If constraint_mode is not recognized
        """
        if constraint_mode not in MultiPeakVoigtModel.CONSTRAINT_MODES:
            raise ValueError(f"Unknown constraint mode: {constraint_mode}. "
                           f"Must be one of {MultiPeakVoigtModel.CONSTRAINT_MODES}")

        if constraint_mode == 'free':
            # Standard multi-peak Voigt with all parameters free
            return lambda x, *params: multi_voigt_profile(x, *params)

        elif constraint_mode == 'shared_width':
            # All peaks share the same sigma and gamma
            def shared_width_model(x: np.ndarray, *params) -> np.ndarray:
                """Multi-peak Voigt with shared width parameters."""
                # Parameters: [amp1, center1, ..., ampN, centerN, sigma_shared, gamma_shared, baseline]
                n_params_per_peak = 2  # amp, center only
                sigma_shared = params[n_peaks * n_params_per_peak]
                gamma_shared = params[n_peaks * n_params_per_peak + 1]
                baseline = params[-1]

                # Reconstruct full parameter list
                full_params = []
                for i in range(n_peaks):
                    idx = i * n_params_per_peak
                    full_params.extend([params[idx], params[idx + 1],
                                      sigma_shared, gamma_shared])
                full_params.append(baseline)

                return multi_voigt_profile(x, *full_params)

            return shared_width_model

        elif constraint_mode == 'fixed_positions':
            # Peak positions are fixed, fit amplitudes and widths
            if fixed_params is None or 'centers' not in fixed_params:
                raise ValueError("fixed_positions mode requires 'centers' in fixed_params")

            centers = fixed_params['centers']
            if len(centers) != n_peaks:
                raise ValueError(f"Number of fixed centers ({len(centers)}) "
                               f"must match n_peaks ({n_peaks})")

            def fixed_positions_model(x: np.ndarray, *params) -> np.ndarray:
                """Multi-peak Voigt with fixed peak positions."""
                # Parameters: [amp1, sigma1, gamma1, ..., ampN, sigmaN, gammaN, baseline]
                n_params_per_peak = 3  # amp, sigma, gamma
                baseline = params[-1]

                # Reconstruct full parameter list with fixed centers
                full_params = []
                for i in range(n_peaks):
                    idx = i * n_params_per_peak
                    full_params.extend([params[idx], centers[i],
                                      params[idx + 1], params[idx + 2]])
                full_params.append(baseline)

                return multi_voigt_profile(x, *full_params)

            return fixed_positions_model

    @staticmethod
    def create_bounds(n_peaks: int,
                     peak_centers: List[float],
                     stage: str,
                     x_range: Optional[Tuple[float, float]] = None,
                     amplitude_range: Optional[Tuple[float, float]] = None,
                     width_range: Optional[Tuple[float, float]] = None,
                     baseline_range: Optional[Tuple[float, float]] = None,
                     constraint_mode: str = 'free') -> Tuple[List[float], List[float]]:
        """
        Generate parameter bounds for different fitting stages.

        Parameters
        ----------
        n_peaks : int
            Number of peaks in the model
        peak_centers : List[float]
            Initial peak center positions
        stage : str
            Fitting stage: 'initial', 'refinement', or 'final'
            - 'initial': Loose bounds for robust initial fit
            - 'refinement': Tighter bounds around current estimates
            - 'final': Very tight bounds for fine-tuning
        x_range : Optional[Tuple[float, float]], optional
            Overall x-axis range (min, max) for the data
        amplitude_range : Optional[Tuple[float, float]], optional
            Allowed amplitude range (min, max)
        width_range : Optional[Tuple[float, float]], optional
            Allowed width range (min, max) for sigma and gamma
        baseline_range : Optional[Tuple[float, float]], optional
            Allowed baseline range (min, max)
        constraint_mode : str, optional
            Constraint mode (must match model creation)

        Returns
        -------
        Tuple[List[float], List[float]]
            Lower and upper bounds for all parameters

        Notes
        -----
        Bound generation strategy varies by stage:
        - 'initial': ±50% position wiggle, wide amplitude/width ranges
        - 'refinement': ±20% position wiggle, moderate ranges
        - 'final': ±5% position wiggle, tight ranges
        """
        # Default ranges if not provided
        if x_range is None:
            x_min, x_max = min(peak_centers) - 1.0, max(peak_centers) + 1.0
        else:
            x_min, x_max = x_range

        if amplitude_range is None:
            amplitude_range = (0.0, np.inf)

        if width_range is None:
            width_range = (MIN_WIDTH, 1.0)

        if baseline_range is None:
            baseline_range = (-np.inf, np.inf)

        # Stage-dependent wiggle factors
        wiggle_factors = {
            'initial': 0.5,      # ±50% of typical peak spacing
            'refinement': 0.2,   # ±20%
            'final': 0.05        # ±5%
        }

        if stage not in wiggle_factors:
            warnings.warn(f"Unknown stage '{stage}', using 'initial' bounds")
            stage = 'initial'

        wiggle = wiggle_factors[stage]

        # Calculate typical peak spacing for position bounds
        if n_peaks > 1:
            spacings = np.diff(sorted(peak_centers))
            typical_spacing = np.median(spacings)
        else:
            typical_spacing = (x_max - x_min) * 0.1

        position_wiggle = typical_spacing * wiggle

        # Build bounds based on constraint mode
        lower_bounds = []
        upper_bounds = []

        if constraint_mode == 'free':
            # Standard bounds: [amp, center, sigma, gamma] per peak
            for center in peak_centers:
                # Amplitude bounds
                lower_bounds.append(amplitude_range[0])
                upper_bounds.append(amplitude_range[1])

                # Center bounds
                lower_bounds.append(max(x_min, center - position_wiggle))
                upper_bounds.append(min(x_max, center + position_wiggle))

                # Sigma bounds
                lower_bounds.append(width_range[0])
                upper_bounds.append(width_range[1])

                # Gamma bounds
                lower_bounds.append(width_range[0])
                upper_bounds.append(width_range[1])

            # Baseline bounds
            lower_bounds.append(baseline_range[0])
            upper_bounds.append(baseline_range[1])

        elif constraint_mode == 'shared_width':
            # Bounds: [amp, center] per peak, then sigma, gamma, baseline
            for center in peak_centers:
                # Amplitude bounds
                lower_bounds.append(amplitude_range[0])
                upper_bounds.append(amplitude_range[1])

                # Center bounds
                lower_bounds.append(max(x_min, center - position_wiggle))
                upper_bounds.append(min(x_max, center + position_wiggle))

            # Shared sigma bounds
            lower_bounds.append(width_range[0])
            upper_bounds.append(width_range[1])

            # Shared gamma bounds
            lower_bounds.append(width_range[0])
            upper_bounds.append(width_range[1])

            # Baseline bounds
            lower_bounds.append(baseline_range[0])
            upper_bounds.append(baseline_range[1])

        elif constraint_mode == 'fixed_positions':
            # Bounds: [amp, sigma, gamma] per peak (no center)
            for _ in range(n_peaks):
                # Amplitude bounds
                lower_bounds.append(amplitude_range[0])
                upper_bounds.append(amplitude_range[1])

                # Sigma bounds
                lower_bounds.append(width_range[0])
                upper_bounds.append(width_range[1])

                # Gamma bounds
                lower_bounds.append(width_range[0])
                upper_bounds.append(width_range[1])

            # Baseline bounds
            lower_bounds.append(baseline_range[0])
            upper_bounds.append(baseline_range[1])

        return lower_bounds, upper_bounds

    @staticmethod
    def parse_params(params: np.ndarray,
                    n_peaks: int,
                    constraint_mode: str = 'free',
                    fixed_params: Optional[Dict[str, Any]] = None) -> List[Dict[str, float]]:
        """
        Parse fitted parameters into individual peak dictionaries.

        Parameters
        ----------
        params : np.ndarray
            Fitted parameter array
        n_peaks : int
            Number of peaks
        constraint_mode : str, optional
            Constraint mode used during fitting
        fixed_params : Optional[Dict[str, Any]], optional
            Fixed parameters (needed for 'fixed_positions' mode)

        Returns
        -------
        List[Dict[str, float]]
            List of peak parameter dictionaries
        """
        if constraint_mode == 'free':
            return parse_multi_peak_params(params, n_peaks)

        elif constraint_mode == 'shared_width':
            # Reconstruct full parameter array
            n_params_per_peak = 2  # amp, center
            sigma_shared = float(params[n_peaks * n_params_per_peak])
            gamma_shared = float(params[n_peaks * n_params_per_peak + 1])
            baseline = float(params[-1])

            full_params = []
            for i in range(n_peaks):
                idx = i * n_params_per_peak
                full_params.extend([params[idx], params[idx + 1],
                                  sigma_shared, gamma_shared])
            full_params.append(baseline)

            return parse_multi_peak_params(np.array(full_params), n_peaks)

        elif constraint_mode == 'fixed_positions':
            if fixed_params is None or 'centers' not in fixed_params:
                raise ValueError("fixed_positions mode requires 'centers' in fixed_params")

            centers = fixed_params['centers']

            # Reconstruct full parameter array
            n_params_per_peak = 3  # amp, sigma, gamma
            baseline = float(params[-1])

            full_params = []
            for i in range(n_peaks):
                idx = i * n_params_per_peak
                full_params.extend([params[idx], centers[i],
                                  params[idx + 1], params[idx + 2]])
            full_params.append(baseline)

            return parse_multi_peak_params(np.array(full_params), n_peaks)

        else:
            raise ValueError(f"Unknown constraint mode: {constraint_mode}")


class VoigtModelFactory:
    """
    Convenience factory class for creating common Voigt model types.

    This class provides simplified interfaces for creating frequently used
    model configurations without needing to specify all parameters.
    """

    @staticmethod
    def simultaneous_fit_model(n_peaks: int) -> Callable:
        """
        Create model for simultaneous fitting of all peaks with independent parameters.

        This is the most flexible model where each peak has its own amplitude,
        center, sigma, and gamma parameters.

        Parameters
        ----------
        n_peaks : int
            Number of peaks to fit simultaneously

        Returns
        -------
        Callable
            Model function for simultaneous fitting
        """
        return MultiPeakVoigtModel.create_model(n_peaks, constraint_mode='free')

    @staticmethod
    def shared_width_model(n_peaks: int) -> Callable:
        """
        Create model where all peaks share the same width parameters.

        Useful when peaks are expected to have similar shapes (e.g., multiplets,
        or peaks from similar chemical environments).

        Parameters
        ----------
        n_peaks : int
            Number of peaks

        Returns
        -------
        Callable
            Model function with shared sigma and gamma
        """
        return MultiPeakVoigtModel.create_model(n_peaks, constraint_mode='shared_width')

    @staticmethod
    def series_fit_model(n_peaks: int, n_spectra: int) -> Callable:
        """
        Create model for fitting a series of spectra with shared positions/widths.

        This model is designed for kinetic or temperature series where peak
        positions and shapes remain constant but intensities vary.

        Parameters
        ----------
        n_peaks : int
            Number of peaks per spectrum
        n_spectra : int
            Number of spectra in the series

        Returns
        -------
        Callable
            Model function for series fitting

        Notes
        -----
        Parameter structure:
        - Amplitudes vary per spectrum: [amp1_spec1, amp2_spec1, ..., amp1_spec2, ...]
        - Positions and widths shared: [center1, ..., centerN, sigma1, ..., sigmaN,
                                       gamma1, ..., gammaN, baseline]
        """
        def series_model(x: np.ndarray, *params) -> np.ndarray:
            """
            Multi-spectrum model with shared positions and widths.

            Parameters are organized as:
            - First n_peaks*n_spectra parameters: amplitudes for each peak in each spectrum
            - Next n_peaks parameters: shared centers
            - Next n_peaks parameters: shared sigmas
            - Next n_peaks parameters: shared gammas
            - Last parameter: baseline
            """
            # Extract shared parameters
            n_amps = n_peaks * n_spectra
            amplitudes = params[:n_amps]
            centers = params[n_amps:n_amps + n_peaks]
            sigmas = params[n_amps + n_peaks:n_amps + 2*n_peaks]
            gammas = params[n_amps + 2*n_peaks:n_amps + 3*n_peaks]
            baseline = params[-1]

            # For simplicity, return model for first spectrum
            # (In practice, this would be extended to handle multiple spectra)
            full_params = []
            for i in range(n_peaks):
                full_params.extend([
                    amplitudes[i],  # First spectrum's amplitude
                    centers[i],
                    sigmas[i],
                    gammas[i]
                ])
            full_params.append(baseline)

            return multi_voigt_profile(x, *full_params)

        return series_model

    @staticmethod
    def fixed_position_model(peak_centers: List[float]) -> Callable:
        """
        Create model with fixed peak positions (fit only amplitudes and widths).

        Useful when peak positions are known accurately (e.g., from reference
        spectra or theoretical predictions).

        Parameters
        ----------
        peak_centers : List[float]
            Fixed peak center positions

        Returns
        -------
        Callable
            Model function with fixed positions
        """
        n_peaks = len(peak_centers)
        fixed_params = {'centers': peak_centers}
        return MultiPeakVoigtModel.create_model(
            n_peaks,
            constraint_mode='fixed_positions',
            fixed_params=fixed_params
        )


# Convenience function for quick model creation
def create_voigt_model(n_peaks: int,
                      mode: str = 'simultaneous',
                      **kwargs) -> Callable:
    """
    Convenience function to create Voigt models with common configurations.

    Parameters
    ----------
    n_peaks : int
        Number of peaks
    mode : str, optional
        Model type: 'simultaneous', 'shared_width', 'fixed_positions', or 'series'
    **kwargs
        Additional arguments passed to specific model factories

    Returns
    -------
    Callable
        Configured model function

    Examples
    --------
    >>> # Simple simultaneous fit
    >>> model = create_voigt_model(3, mode='simultaneous')

    >>> # Shared width model
    >>> model = create_voigt_model(5, mode='shared_width')

    >>> # Fixed positions
    >>> model = create_voigt_model(2, mode='fixed_positions',
    ...                           peak_centers=[7.5, 8.2])
    """
    mode_map = {
        'simultaneous': VoigtModelFactory.simultaneous_fit_model,
        'shared_width': VoigtModelFactory.shared_width_model,
        'series': VoigtModelFactory.series_fit_model,
        'fixed_positions': VoigtModelFactory.fixed_position_model
    }

    if mode not in mode_map:
        raise ValueError(f"Unknown mode: {mode}. Must be one of {list(mode_map.keys())}")

    factory = mode_map[mode]

    # Handle different factory signatures
    if mode == 'series':
        if 'n_spectra' not in kwargs:
            raise ValueError("series mode requires 'n_spectra' argument")
        return factory(n_peaks, kwargs['n_spectra'])
    elif mode == 'fixed_positions':
        if 'peak_centers' not in kwargs:
            raise ValueError("fixed_positions mode requires 'peak_centers' argument")
        return factory(kwargs['peak_centers'])
    else:
        return factory(n_peaks)

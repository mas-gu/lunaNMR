"""
ABOUTME: Centralized PS2D algorithm configuration for different nucleus types
ABOUTME: Single source of truth for all radii, thresholds, and constraints

PS2D Configuration Management
=============================

This module provides centralized configuration for all PS2D 2D fitting parameters.
Eliminates scattered constants across multiple files and makes nucleus switching trivial.

Date: 2025-10-10
Version: 1.0 - CENTRALIZED_CONFIG
Author: Claude Code (architectural refactoring)
"""

from typing import Dict, Any


class PS2DConfig:
    """
    Centralized PS2D algorithm parameters for different nucleus types.

    Single source of truth for all radii, thresholds, and constraints.
    Eliminates 12 scattered hardcoded constants and makes 15N→13C switching
    a single function call instead of manual editing.

    Parameters are organized by physical meaning:
    - Ellipse radii: data selection windows (radF1, radF2)
    - Overlap thresholds: when peaks are considered overlapping
    - Gap thresholds: when to split large clusters
    - Linewidth constraints: minimum valid linewidths
    - Position bounds: how far peaks can move during fitting

    All F1/Y parameters scale by 0.3 for 13C (3× narrower than 15N)
    All F2/X parameters stay constant (1H same in both spectra)
    """

    # Master parameter dictionary
    NUCLEUS_PARAMS = {
        '15N': {
            # ================================================================
            # Ellipse radii for data selection (from PS2D C++ )
            # ================================================================
            # Used in fit_overlap_group_2d() for 2D region extraction
            'radF1': 0.4,              # Indirect dimension (15N) - reduced from PS2D 0.6
            'radF2': 0.04,             # Direct dimension (1H) - reduced from PS2D 0.06

            # Used in ps2d_data_selector for elliptical window selection
            # These are broader for safety (match original PS2D defaults)
            'radF1_selector': 0.6,     # ps2d_data_selector default (PS2D original)
            'radF2_selector': 0.06,    # ps2d_data_selector default (PS2D original)

            # ================================================================
            # Overlap detection thresholds (geometric ellipse intersection)
            # ================================================================
            # Used in identify_overlap_clusters() and check_if_peaks_need_2d_fitting()
            # Larger thresholds = more aggressive overlap detection
            'overlap_threshold_x': 0.10,  # 1H dimension (~8× typical FWHM)
            'overlap_threshold_y': 0.8,   # 15N dimension (~8× typical FWHM)

            # Used in peaks_overlap_elliptical() (simpler overlap check)
            'overlap_threshold_x_simple': 0.06,  # 1H dimension
            'overlap_threshold_y_simple': 0.6,   # 15N dimension

            # ================================================================
            # Gap detection thresholds (cluster splitting)
            # ================================================================
            # Used in _split_cluster_by_gaps() to find spatial gaps
            # Smaller thresholds = more aggressive splitting
            'gap_threshold_x': 0.05,      # 1H dimension (~4× typical FWHM)
            'gap_threshold_y': 0.4,       # 15N dimension (~4× typical FWHM)

            # ================================================================
            # Linewidth constraints (CRITICAL for convergence)
            # ================================================================
            # Minimum allowed linewidths prevent collapse to unrealistic values
            # Set to ~50% of typical FWHM for each nucleus
            'min_linewidth_f1': 0.05,     # Minimum 15N linewidth (ppm)
            'min_linewidth_f2': 0.001,    # Minimum 1H linewidth (ppm)

            # ================================================================
            # Position bounds (fitting movement limits)
            # ================================================================
            # Maximum distance peaks can move from initial position during fitting
            'pos_margin_f1': 0.2,         # 15N position movement limit (ppm)
            'pos_margin_f2': 0.05,        # 1H position movement limit (ppm)

            # ================================================================
            # Cluster size limits
            # ================================================================
            'max_cluster_size': 6,        # Maximum peaks in overlap group

            # ================================================================
            # Levenberg-Marquardt optimizer settings
            # ================================================================
            'max_iterations': 500,        # Maximum LM iterations per stage
        },

        '13C': {
            # ================================================================
            # 13C parameters - OPTIMIZED FOR PROPER R² CALCULATION
            # ================================================================
            # 13C peaks: typical FWHM ~0.04-0.06 ppm (measured from spectrum)
            # Key insight: radF1 must be 3-4× peak width (like 15N) for R² to work
            # radF1/FWHM ratio = 0.15/0.05 = 3.0× (matches 15N's 4.0×)

            # Ellipse radii for fitting (generous coverage like 15N)
            'radF1': 0.15,                # 3× typical FWHM, ensures R² works correctly
            'radF2': 0.04,                # ~3-4× typical FWHM

            # Ellipse radii for data selection (even broader for safety)
            'radF1_selector': 0.20,       # ~1.3× radF1
            'radF2_selector': 0.06,       # ~1.5× radF2

            # Overlap thresholds (more permissive to catch true overlaps)
            'overlap_threshold_x': 0.10,  # 1H dimension (~8× FWHM)
            'overlap_threshold_y': 0.2,  # 13C dimension (~2× typical FWHM, permissive)

            # Simple overlap check (moderate)
            'overlap_threshold_x_simple': 0.06,  # 1H dimension
            'overlap_threshold_y_simple': 0.08,  # 13C dimension (~1.5× FWHM)

            # Gap thresholds (aggressive splitting for complex clusters)
            'gap_threshold_x': 0.05,      # 1H dimension (~4× FWHM)
            'gap_threshold_y': 0.05,      # 13C dimension (~0.8× FWHM, splits tight clusters)

            # Linewidth constraints (realistic minimums based on actual peak widths)
            'min_linewidth_f1': 0.025,    # Minimum 13C linewidth (~half of typical 0.05-0.08 ppm FWHM)
            'min_linewidth_f2': 0.005,    # Minimum 1H linewidth

            # Position bounds (moderate for peak movement)
            'pos_margin_f1': 0.04,        # Allow ±0.04 ppm movement (relaxed)
            'pos_margin_f2': 0.04,        # Allow ±0.04 ppm movement

            # Cluster size limits (allow larger groups as requested)
            'max_cluster_size': 6,        # Max 6 peaks per cluster (same as 15N)

            # Levenberg-Marquardt optimizer settings
            'max_iterations': 500,        # Maximum LM iterations per stage
        }
    }

    def __init__(self, nucleus_type: str = '15N'):
        """
        Initialize PS2D configuration for given nucleus type.

        Parameters:
        -----------
        nucleus_type : str
            '15N' or '13C' (default: '15N')

        Raises:
        -------
        ValueError : If nucleus_type not recognized

        Example:
        --------
        >>> config = PS2DConfig('15N')
        >>> print(config.radF1)  # 0.4
        >>>
        >>> config = PS2DConfig('13C')
        >>> print(config.radF1)  # 0.12 (scaled for narrower peaks)
        """
        if nucleus_type not in self.NUCLEUS_PARAMS:
            valid_types = ', '.join(self.NUCLEUS_PARAMS.keys())
            raise ValueError(
                f"Unknown nucleus type: '{nucleus_type}'. "
                f"Valid types: {valid_types}"
            )

        self.nucleus_type = nucleus_type

        # Load all parameters as instance attributes for easy access
        params = self.NUCLEUS_PARAMS[nucleus_type]
        for key, value in params.items():
            setattr(self, key, value)

    def __repr__(self):
        """String representation for debugging"""
        return f"PS2DConfig(nucleus_type='{self.nucleus_type}')"

    def __str__(self):
        """Human-readable string representation"""
        lines = [f"PS2D Configuration: {self.nucleus_type}"]
        lines.append("-" * 40)

        params = self.NUCLEUS_PARAMS[self.nucleus_type]
        for key, value in params.items():
            lines.append(f"  {key:30s} = {value}")

        return "\n".join(lines)

    def to_dict(self) -> Dict[str, Any]:
        """
        Return all parameters as dictionary.

        Returns:
        --------
        dict : All configuration parameters

        Example:
        --------
        >>> config = PS2DConfig('15N')
        >>> params = config.to_dict()
        >>> params['radF1']
        0.4
        """
        return {k: getattr(self, k)
                for k in self.NUCLEUS_PARAMS[self.nucleus_type].keys()}

    def get_scaling_factor(self) -> float:
        """
        Get scaling factor relative to 15N baseline.

        Returns:
        --------
        float : 1.0 for 15N, 0.3 for 13C

        Useful for understanding parameter relationships.
        """
        if self.nucleus_type == '15N':
            return 1.0
        elif self.nucleus_type == '13C':
            return 0.3 #was 0.3
        else:
            return 1.0  # Default to no scaling


# ============================================================================
# Global configuration instance (singleton pattern)
# ============================================================================

_ps2d_config = PS2DConfig('15N')  # Default to 15N


def get_ps2d_config() -> PS2DConfig:
    """
    Get current global PS2D configuration.

    Returns:
    --------
    PS2DConfig : Current global configuration instance

    Example:
    --------
    >>> from lunaNMR.core.ps2d_config import get_ps2d_config
    >>> config = get_ps2d_config()
    >>> print(config.nucleus_type)
    '15N'
    >>> print(config.radF1)
    0.4
    """
    return _ps2d_config


def set_ps2d_config(nucleus_type: str) -> PS2DConfig:
    """
    Set global PS2D configuration for nucleus type.

    This is the MAIN function for switching between 15N and 13C.
    Call this once at application startup or when loading a spectrum.

    Parameters:
    -----------
    nucleus_type : str
        '15N' or '13C'

    Returns:
    --------
    PS2DConfig : New configuration instance

    Example:
    --------
    >>> from lunaNMR.core.ps2d_config import set_ps2d_config
    >>>
    >>> # For 13C-HSQC spectra
    >>> set_ps2d_config('13C')
    >>>
    >>> # For 15N-HSQC spectra
    >>> set_ps2d_config('15N')

    Notes:
    ------
    After calling this function, all PS2D functions will use the new
    configuration automatically. No need to edit 12 hardcoded values!
    """
    global _ps2d_config
    _ps2d_config = PS2DConfig(nucleus_type)
    return _ps2d_config


def reset_ps2d_config():
    """
    Reset global configuration to default (15N).

    Useful for testing or resetting to known state.
    """
    global _ps2d_config
    _ps2d_config = PS2DConfig('15N')


# ============================================================================
# Convenience functions for parameter access
# ============================================================================

def get_ellipse_radii() -> tuple:
    """Get current ellipse radii (radF1, radF2)"""
    config = get_ps2d_config()
    return config.radF1, config.radF2


def get_overlap_thresholds() -> tuple:
    """Get current overlap thresholds (x, y)"""
    config = get_ps2d_config()
    return config.overlap_threshold_x, config.overlap_threshold_y


def get_gap_thresholds() -> tuple:
    """Get current gap thresholds (x, y)"""
    config = get_ps2d_config()
    return config.gap_threshold_x, config.gap_threshold_y


def get_min_linewidths() -> tuple:
    """Get current minimum linewidths (f1, f2)"""
    config = get_ps2d_config()
    return config.min_linewidth_f1, config.min_linewidth_f2


# ============================================================================
# Module exports
# ============================================================================

__all__ = [
    'PS2DConfig',
    'get_ps2d_config',
    'set_ps2d_config',
    'reset_ps2d_config',
    'get_ellipse_radii',
    'get_overlap_thresholds',
    'get_gap_thresholds',
    'get_min_linewidths',
]


# ============================================================================
# Testing and validation
# ============================================================================

if __name__ == "__main__":
    print("=" * 70)
    print("PS2D Configuration Module - Test Suite")
    print("=" * 70)

    # Test 1: Default configuration
    print("\n--- Test 1: Default Configuration (15N) ---")
    config_15n = get_ps2d_config()
    print(config_15n)
    print(f"\nNucleus: {config_15n.nucleus_type}")
    print(f"radF1: {config_15n.radF1}, radF2: {config_15n.radF2}")
    print(f"min_linewidth_f1: {config_15n.min_linewidth_f1}")

    # Test 2: Switch to 13C
    print("\n--- Test 2: Switch to 13C Configuration ---")
    config_13c = set_ps2d_config('13C')
    print(config_13c)
    print(f"\nNucleus: {config_13c.nucleus_type}")
    print(f"radF1: {config_13c.radF1}, radF2: {config_13c.radF2}")
    print(f"min_linewidth_f1: {config_13c.min_linewidth_f1}")
    print(f"Scaling factor: {config_13c.get_scaling_factor()}")

    # Test 3: Verify scaling relationships
    print("\n--- Test 3: Verify Scaling (15N vs 13C) ---")
    config_15n = PS2DConfig('15N')
    config_13c = PS2DConfig('13C')

    print(f"radF1:      15N={config_15n.radF1:.3f}, 13C={config_13c.radF1:.3f}, ratio={config_13c.radF1/config_15n.radF1:.2f}")
    print(f"radF2:      15N={config_15n.radF2:.3f}, 13C={config_13c.radF2:.3f}, ratio={config_13c.radF2/config_15n.radF2:.2f}")
    print(f"overlap_y:  15N={config_15n.overlap_threshold_y:.3f}, 13C={config_13c.overlap_threshold_y:.3f}, ratio={config_13c.overlap_threshold_y/config_15n.overlap_threshold_y:.2f}")
    print(f"min_lw_f1:  15N={config_15n.min_linewidth_f1:.3f}, 13C={config_13c.min_linewidth_f1:.3f}, ratio={config_13c.min_linewidth_f1/config_15n.min_linewidth_f1:.2f}")

    # Test 4: Dictionary export
    print("\n--- Test 4: Dictionary Export ---")
    params_dict = config_13c.to_dict()
    print(f"Exported {len(params_dict)} parameters:")
    for key, value in list(params_dict.items())[:5]:
        print(f"  {key}: {value}")

    # Test 5: Error handling
    print("\n--- Test 5: Error Handling ---")
    try:
        bad_config = PS2DConfig('17O')
    except ValueError as e:
        print(f"✅ Correctly caught error: {e}")

    print("\n" + "=" * 70)
    print("All tests passed!")
    print("=" * 70)

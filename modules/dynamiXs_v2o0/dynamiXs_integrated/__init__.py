"""
DynamiXs Integrated Spectral Density Analysis Module

This module provides an automated workflow for complete spectral density analysis
from raw NMR relaxation data (T1, T2, hetNOE intensities) to model-free parameters.
"""

from .integrated_analysis import (
    IntegratedAnalysisParameters,
    IntegratedAnalysisPipeline
)

from .data_converters import (
    convert_relaxation_times_to_rates,
    calculate_hetnoe_from_intensities,
    parse_intensity_csv,
    create_spectral_density_input_csv
)

from .data_validation import (
    merge_datasets_by_residue,
    validate_relaxation_rates,
    validate_field_consistency,
    check_minimum_residue_count,
    generate_validation_report,
    filter_invalid_residues
)

from .fitting_wrapper import (
    run_t1_fitting,
    run_t2_fitting
)

__all__ = [
    # Main pipeline classes
    'IntegratedAnalysisParameters',
    'IntegratedAnalysisPipeline',

    # Data converters
    'convert_relaxation_times_to_rates',
    'calculate_hetnoe_from_intensities',
    'parse_intensity_csv',
    'create_spectral_density_input_csv',

    # Data validation
    'merge_datasets_by_residue',
    'validate_relaxation_rates',
    'validate_field_consistency',
    'check_minimum_residue_count',
    'generate_validation_report',
    'filter_invalid_residues',

    # Fitting wrappers
    'run_t1_fitting',
    'run_t2_fitting',
]

__version__ = '1.0.0'
__author__ = 'DynamiXs Development Team'

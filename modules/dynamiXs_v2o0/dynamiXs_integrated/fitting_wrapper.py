#!/usr/bin/env python3
"""
Fitting Wrapper for T1/T2 Exponential Decay Analysis

This module provides a high-level wrapper around the existing fit_Tx_NMRRE.py
scripts to enable programmatic fitting without command-line interaction.

Author: DynamiXs Development Team
Date: 2025-01-23
"""

import sys
import os
import pandas as pd
import numpy as np
from pathlib import Path
from typing import Dict, Tuple, Optional
import tempfile


def fit_relaxation_data(input_csv: str,
                        experiment_type: str = 'T1',
                        initial_amplitude: float = 5.0,
                        initial_time_constant: float = 100.0,
                        n_bootstrap: int = 1000,
                        error_method: str = 'analytical',
                        output_prefix: Optional[str] = None,
                        json_folder: Optional[str] = None,
                        field_name: str = 'field1',
                        field_freq: float = 600.0) -> Tuple[Dict[str, Dict], str]:
    """
    Fit T1 or T2 relaxation data from a CSV file.

    This is a wrapper around the existing dynamiXs_T1_T2 fitting scripts,
    providing a programmatic interface for use in the integrated workflow.

    Parameters
    ----------
    input_csv : str
        Path to input CSV file with format:
        , time1, time2, time3, ...
        residue1, intensity1, intensity2, intensity3, ...
        residue2, intensity1, intensity2, intensity3, ...
    experiment_type : str
        Type of experiment: 'T1', 'T2', or 'T1rho'
    initial_amplitude : float
        Initial guess for amplitude parameter
    initial_time_constant : float
        Initial guess for relaxation time (ms)
    n_bootstrap : int
        Number of bootstrap iterations for error estimation
    output_prefix : str, optional
        Prefix for output files (if None, uses temp file)

    Returns
    -------
    tuple
        (fitted_params, output_file)
        - fitted_params: dict mapping residue_id to {'value': T, 'error': T_err}
        - output_file: path to results text file

    Examples
    --------
    >>> t1_params, results_file = fit_relaxation_data('T1_data.csv', 'T1')
    >>> print(f"Fitted {len(t1_params)} residues")
    >>> print(f"Results saved to {results_file}")
    """
    # Import the fitting module (prefer multicore version for performance)
    fitting_module_path = Path(__file__).parent.parent / 'dynamiXs_T1_T2'
    sys.path.insert(0, str(fitting_module_path))

    try:
        # Try multicore version first (5-8x faster)
        from fitmulti__Tx_NMRRE import run_analysis_with_params
        print(f"  ⚡ Using multicore fitting ({experiment_type})")
    except ImportError:
        # Fallback to single-core version if multicore not available
        try:
            from fit_Tx_NMRRE import run_analysis_with_params
            print(f"  ⚠ Using single-core fitting ({experiment_type}) - consider using multicore version for better performance")
        except ImportError as e:
            raise ImportError(
                f"Could not import fitting module from {fitting_module_path}. "
                f"Ensure fit_Tx_NMRRE.py or fitmulti__Tx_NMRRE.py exists. Error: {e}"
            )

    # Create output file paths
    if output_prefix is None:
        temp_dir = Path(__file__).parent.parent / 'temp'
        temp_dir.mkdir(exist_ok=True)
        output_prefix = str(temp_dir / f"temp_{experiment_type}_fit")

    results_file = f"{output_prefix}_fit_results.txt"

    # Prepare parameters dictionary for the fitting function
    params = {
        'input_csv_file': input_csv,
        'output_prefix': output_prefix,
        'results_txt_file': results_file,
        'experiment_type': experiment_type,
        'time_units': 'ms',
        'signal_units': 'Intensity',
        'initial_A': initial_amplitude,
        'initial_t2': initial_time_constant,
        'n_bootstrap': n_bootstrap,
        'error_method': error_method,
        'n_plots_per_figure': 20,
        'json_folder': json_folder,
        'field_name': field_name,
        'field_freq': field_freq
    }

    # Run the fitting analysis
    run_analysis_with_params(params)

    # Parse the results file to extract fitted parameters
    fitted_params = _parse_fit_results_file(results_file, experiment_type)

    return fitted_params, results_file


def _parse_fit_results_file(results_file: str, param_name: str = 'T1') -> Dict[str, Dict]:
    """
    Parse the tab-delimited fit results file.

    Expected format:
        Residue\tA\tT1\tA_err\tT1_err\tSuccess
        1.0\t2.244e+06\t7.875e+02\t6.488e+04\t5.824e+01\tYes

    Parameters
    ----------
    results_file : str
        Path to results file
    param_name : str
        Name of the relaxation parameter ('T1', 'T2', 'T1rho')

    Returns
    -------
    dict
        Dictionary mapping residue_id to {'value': T, 'error': T_err}
    """
    if not Path(results_file).exists():
        raise FileNotFoundError(f"Results file not found: {results_file}")

    # Read tab-delimited file
    df = pd.read_csv(results_file, sep='\t')

    fitted_params = {}

    for _, row in df.iterrows():
        residue = str(row['Residue'])
        success = row['Success']

        if success == 'Yes':
            # Extract the relaxation time and error
            value = float(row[param_name])
            error = float(row[f'{param_name}_err'])

            fitted_params[residue] = {
                'value': value,
                'error': error
            }

    return fitted_params


def run_t1_fitting(input_csv: str,
                   initial_amplitude: float = 5.0,
                   initial_t1: float = 800.0,
                   n_bootstrap: int = 1000,
                   error_method: str = 'analytical',
                   output_prefix: Optional[str] = None,
                   json_folder: Optional[str] = None,
                   field_name: str = 'field1',
                   field_freq: float = 600.0) -> Tuple[Dict[str, Dict], str]:
    """
    Convenience function for T1 fitting.

    Parameters
    ----------
    input_csv : str
        Path to T1 data CSV
    initial_amplitude : float
        Initial amplitude guess
    initial_t1 : float
        Initial T1 guess (ms)
    n_bootstrap : int
        Bootstrap iterations
    error_method : str
        Error estimation method ('analytical' or 'bootstrap')
    output_prefix : str, optional
        Output file prefix
    json_folder : str, optional
        Folder for JSON fit data
    field_name : str
        Field identifier ('field1' or 'field2')
    field_freq : float
        Magnetic field frequency in MHz

    Returns
    -------
    tuple
        (t1_params, results_file)
    """
    return fit_relaxation_data(
        input_csv=input_csv,
        experiment_type='T1',
        initial_amplitude=initial_amplitude,
        initial_time_constant=initial_t1,
        n_bootstrap=n_bootstrap,
        error_method=error_method,
        output_prefix=output_prefix,
        json_folder=json_folder,
        field_name=field_name,
        field_freq=field_freq
    )


def run_t2_fitting(input_csv: str,
                   initial_amplitude: float = 5.0,
                   initial_t2: float = 100.0,
                   n_bootstrap: int = 1000,
                   error_method: str = 'analytical',
                   output_prefix: Optional[str] = None,
                   json_folder: Optional[str] = None,
                   field_name: str = 'field1',
                   field_freq: float = 600.0) -> Tuple[Dict[str, Dict], str]:
    """
    Convenience function for T2 fitting.

    Parameters
    ----------
    input_csv : str
        Path to T2 data CSV
    initial_amplitude : float
        Initial amplitude guess
    initial_t2 : float
        Initial T2 guess (ms)
    n_bootstrap : int
        Bootstrap iterations
    error_method : str
        Error estimation method ('analytical' or 'bootstrap')
    output_prefix : str, optional
        Output file prefix
    json_folder : str, optional
        Folder for JSON fit data
    field_name : str
        Field identifier ('field1' or 'field2')
    field_freq : float
        Magnetic field frequency in MHz

    Returns
    -------
    tuple
        (t2_params, results_file)
    """
    return fit_relaxation_data(
        input_csv=input_csv,
        experiment_type='T2',
        initial_amplitude=initial_amplitude,
        initial_time_constant=initial_t2,
        n_bootstrap=n_bootstrap,
        error_method=error_method,
        output_prefix=output_prefix,
        json_folder=json_folder,
        field_name=field_name,
        field_freq=field_freq
    )


# NOTE: The current implementation assumes the fitting scripts output
# results to a text file. For a fully integrated solution, we would need
# to either:
# 1. Refactor the existing fit_Tx_NMRRE.py to accept parameters and return results
# 2. Or use subprocess to call the script and parse outputs
#
# For now, this provides the interface structure for the integrated workflow.

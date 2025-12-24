#!/usr/bin/env python3
"""
Data Converters for Integrated Spectral Density Analysis

This module provides functions to convert between different data formats used in
the NMR relaxation analysis pipeline:
- T1/T2 times → R1/R2 rates with error propagation
- hetNOE calculation from saturated/unsaturated intensities
- Data format transformations for spectral density analysis

Author: DynamiXs Development Team
Date: 2025-01-23
"""

import numpy as np
import pandas as pd
from typing import Dict, Tuple, Optional
from pathlib import Path


def convert_relaxation_times_to_rates(t_values: Dict[str, Dict[str, float]],
                                      time_units: str = 'ms') -> Dict[str, Dict[str, float]]:
    """
    Convert relaxation times (T1 or T2) to relaxation rates (R1 or R2) with error propagation.

    Formula:
        R = 1000 / T    (if T in ms, R in s⁻¹)
        R_err = 1000 * T_err / T²

    Parameters
    ----------
    t_values : dict
        Dictionary mapping residue IDs to {'value': T, 'error': T_err}
    time_units : str
        Units of input time ('ms', 's', 'us')

    Returns
    -------
    dict
        Dictionary mapping residue IDs to {'value': R, 'error': R_err}

    Examples
    --------
    >>> t1_data = {'3.LYS': {'value': 832.0, 'error': 66.0}}
    >>> r1_data = convert_relaxation_times_to_rates(t1_data, 'ms')
    >>> print(f"R1 = {r1_data['3.LYS']['value']:.3f} ± {r1_data['3.LYS']['error']:.3f} s⁻¹")
    R1 = 1.202 ± 0.095 s⁻¹
    """
    # Conversion factors to seconds
    time_conversion = {
        'ms': 1000.0,  # milliseconds → seconds (R = 1000/T)
        's': 1.0,      # seconds → seconds (R = 1/T)
        'us': 1e6,     # microseconds → seconds (R = 1e6/T)
    }

    if time_units not in time_conversion:
        raise ValueError(f"Unknown time units: {time_units}. Must be 'ms', 's', or 'us'")

    conversion_factor = time_conversion[time_units]

    r_values = {}
    for res_id, data in t_values.items():
        T = data['value']
        T_err = data['error']

        if T <= 0:
            # Skip invalid values
            continue

        # Convert time to rate
        R = conversion_factor / T

        # Error propagation: dR/dT = -conversion_factor / T²
        # R_err = |dR/dT| * T_err = conversion_factor * T_err / T²
        R_err = conversion_factor * T_err / (T ** 2)

        r_values[res_id] = {
            'value': R,
            'error': R_err
        }

    return r_values


def calculate_hetnoe_from_intensities(saturated_data: Dict[str, float],
                                      unsaturated_data: Dict[str, float],
                                      saturated_errors: Optional[Dict[str, float]] = None,
                                      unsaturated_errors: Optional[Dict[str, float]] = None) -> Dict[str, Dict[str, float]]:
    """
    Calculate heteronuclear NOE values from saturated and unsaturated peak intensities.

    Formula:
        hetNOE = I_saturated / I_unsaturated
        hetNOE_err = hetNOE * sqrt((I_sat_err/I_sat)² + (I_unsat_err/I_unsat)²)

    Parameters
    ----------
    saturated_data : dict
        Peak intensities with ¹H saturation {residue_id: intensity}
    unsaturated_data : dict
        Peak intensities without ¹H saturation {residue_id: intensity}
    saturated_errors : dict, optional
        Errors on saturated intensities
    unsaturated_errors : dict, optional
        Errors on unsaturated intensities

    Returns
    -------
    dict
        Dictionary mapping residue IDs to {'value': hetNOE, 'error': hetNOE_err}

    Examples
    --------
    >>> sat = {'8.TYR': 14585238}
    >>> unsat = {'8.TYR': 14922886}
    >>> noe = calculate_hetnoe_from_intensities(sat, unsat)
    >>> print(f"hetNOE = {noe['8.TYR']['value']:.3f}")
    hetNOE = 0.977
    """
    noe_values = {}

    # Find common residues
    common_residues = set(saturated_data.keys()) & set(unsaturated_data.keys())

    for res_id in common_residues:
        I_sat = saturated_data[res_id]
        I_unsat = unsaturated_data[res_id]

        if I_unsat <= 0:
            # Can't calculate hetNOE with zero or negative reference intensity
            continue

        # Calculate hetNOE
        noe = I_sat / I_unsat

        # Calculate error if provided
        if saturated_errors is not None and unsaturated_errors is not None:
            I_sat_err = saturated_errors.get(res_id, 0)
            I_unsat_err = unsaturated_errors.get(res_id, 0)

            # Error propagation for ratio
            rel_err_sat = I_sat_err / I_sat if I_sat != 0 else 0
            rel_err_unsat = I_unsat_err / I_unsat if I_unsat != 0 else 0

            noe_err = noe * np.sqrt(rel_err_sat**2 + rel_err_unsat**2)
        else:
            # Estimate 2% error if not provided (typical for well-measured peaks)
            noe_err = noe * 0.02

        noe_values[res_id] = {
            'value': noe,
            'error': noe_err
        }

    return noe_values


def parse_intensity_csv(csv_file: str) -> Tuple[Dict[str, float], Optional[Dict[str, float]]]:
    """
    Parse a CSV file containing peak intensities.

    Expected formats:
    1. Single column: residue, intensity
    2. With errors: residue, intensity, intensity_error

    Parameters
    ----------
    csv_file : str
        Path to CSV file

    Returns
    -------
    tuple
        (intensities_dict, errors_dict) where errors_dict is None if not provided

    Examples
    --------
    >>> intensities, errors = parse_intensity_csv('noe_saturated.csv')
    >>> print(f"Found {len(intensities)} residues")
    """
    df = pd.read_csv(csv_file, header=None)

    # Identify columns
    if len(df.columns) < 2:
        raise ValueError(f"CSV must have at least 2 columns (residue, intensity). Found: {df.columns.tolist()}")

    # First column is residue ID
    residue_col = df.columns[0]
    intensity_col = df.columns[1]

    intensities = {}
    errors = None

    for _, row in df.iterrows():
        res_id = str(row[residue_col])
        intensity = float(row[intensity_col])
        intensities[res_id] = intensity

    # Check if error column exists
    if len(df.columns) >= 3:
        error_col = df.columns[2]
        errors = {}
        for _, row in df.iterrows():
            res_id = str(row[residue_col])
            error = float(row[error_col])
            errors[res_id] = error

    return intensities, errors


def create_spectral_density_input_csv(r1_data: Dict[str, Dict[str, float]],
                                       r2_data: Dict[str, Dict[str, float]],
                                       noe_data: Dict[str, Dict[str, float]],
                                       output_file: str,
                                       field_label: Optional[str] = None) -> int:
    """
    Create a CSV file formatted for spectral density analysis scripts.

    Output format:
        Residue,R1,R1err,R2,R2err,hetNOE,hetNOEerr
        8,1.32993,0.06953,9.35377,0.28993,0.84489,0.01984
        ...

    Parameters
    ----------
    r1_data : dict
        R1 values and errors {residue_id: {'value': R1, 'error': R1_err}}
    r2_data : dict
        R2 values and errors {residue_id: {'value': R2, 'error': R2_err}}
    noe_data : dict
        hetNOE values and errors {residue_id: {'value': NOE, 'error': NOE_err}}
    output_file : str
        Path to output CSV file
    field_label : str, optional
        Label to add to column headers (e.g., "f1", "f2")

    Returns
    -------
    int
        Number of residues written to file

    Examples
    --------
    >>> n_residues = create_spectral_density_input_csv(r1, r2, noe, 'field1_input.csv')
    >>> print(f"Created input file with {n_residues} residues")
    """
    # Find common residues (intersection)
    common_residues = set(r1_data.keys()) & set(r2_data.keys()) & set(noe_data.keys())

    # Sort numerically (extract numeric part for sorting)
    import re
    def extract_numeric(res_id):
        """Extract numeric part from residue ID for sorting (e.g., '3.LYS' -> 3, '10' -> 10)"""
        match = re.search(r'\d+', str(res_id))
        return int(match.group()) if match else 0

    common_residues = sorted(common_residues, key=extract_numeric)

    if len(common_residues) == 0:
        raise ValueError("No common residues found across R1, R2, and hetNOE datasets")

    # Create DataFrame
    rows = []
    for res_id in common_residues:
        row = {
            'Residue': res_id,
            'R1': r1_data[res_id]['value'],
            'R1err': r1_data[res_id]['error'],
            'R2': r2_data[res_id]['value'],
            'R2err': r2_data[res_id]['error'],
            'hetNOE': noe_data[res_id]['value'],
            'hetNOEerr': noe_data[res_id]['error']
        }
        rows.append(row)

    df = pd.DataFrame(rows)

    # Add field label to column names if provided
    if field_label:
        rename_map = {
            'R1': f'R1_{field_label}',
            'R1err': f'R1err_{field_label}',
            'R2': f'R2_{field_label}',
            'R2err': f'R2err_{field_label}',
            'hetNOE': f'hetNOE_{field_label}',
            'hetNOEerr': f'hetNOEerr_{field_label}'
        }
        df = df.rename(columns=rename_map)

    # Write to CSV
    df.to_csv(output_file, index=False)

    return len(common_residues)


def extract_residue_id(residue_str: str) -> str:
    """
    Normalize residue identifiers to handle different formats.

    Handles formats like:
    - "3.LYS" → "3"
    - "LYS3" → "3"
    - "3" → "3"
    - "ALA42" → "42"

    Parameters
    ----------
    residue_str : str
        Residue identifier in any format

    Returns
    -------
    str
        Normalized residue number

    Examples
    --------
    >>> extract_residue_id("3.LYS")
    '3'
    >>> extract_residue_id("ALA42")
    '42'
    """
    import re

    # Try to extract number from string
    numbers = re.findall(r'\d+', str(residue_str))

    if numbers:
        return numbers[0]  # Return first number found
    else:
        return str(residue_str)  # Return as-is if no number found


# Utility function for file validation
def validate_csv_structure(csv_file: str, expected_columns: list) -> bool:
    """
    Validate that a CSV file has the expected column structure.

    Parameters
    ----------
    csv_file : str
        Path to CSV file
    expected_columns : list
        List of expected column names

    Returns
    -------
    bool
        True if valid, raises ValueError if not

    Raises
    ------
    ValueError
        If file doesn't exist or doesn't have expected columns
    """
    if not Path(csv_file).exists():
        raise ValueError(f"File not found: {csv_file}")

    try:
        df = pd.read_csv(csv_file, nrows=1)  # Read only header
    except Exception as e:
        raise ValueError(f"Error reading CSV file {csv_file}: {e}")

    missing_cols = set(expected_columns) - set(df.columns)
    if missing_cols:
        raise ValueError(f"CSV {csv_file} missing required columns: {missing_cols}")

    return True

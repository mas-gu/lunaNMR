#!/usr/bin/env python3
"""
Data Validation and Quality Control for Integrated Analysis

This module provides functions to validate NMR relaxation data and perform
quality control checks before spectral density analysis.

Author: DynamiXs Development Team
Date: 2025-01-23
"""

import numpy as np
from typing import Dict, List, Tuple, Set, Optional
import warnings


class DataValidationWarning(UserWarning):
    """Custom warning for data validation issues"""
    pass


def merge_datasets_by_residue(r1_data: Dict[str, Dict[str, float]],
                              r2_data: Dict[str, Dict[str, float]],
                              noe_data: Dict[str, Dict[str, float]],
                              strategy: str = 'intersection') -> Tuple[Dict[str, Dict], List[str]]:
    """
    Merge R1, R2, and hetNOE datasets by residue ID.

    Parameters
    ----------
    r1_data : dict
        R1 values {residue_id: {'value': R1, 'error': R1_err}}
    r2_data : dict
        R2 values {residue_id: {'value': R2, 'error': R2_err}}
    noe_data : dict
        hetNOE values {residue_id: {'value': NOE, 'error': NOE_err}}
    strategy : str
        Merging strategy: 'intersection' (only common residues) or
        'union' (all residues, fill missing with zeros)

    Returns
    -------
    tuple
        (merged_data, excluded_residues)
        - merged_data: dict mapping residue_id to full parameter dict
        - excluded_residues: list of residue IDs that were excluded

    Examples
    --------
    >>> r1 = {'8': {'value': 1.33, 'error': 0.07}}
    >>> r2 = {'8': {'value': 9.35, 'error': 0.29}, '10': {'value': 8.17, 'error': 0.57}}
    >>> noe = {'8': {'value': 0.84, 'error': 0.02}}
    >>> merged, excluded = merge_datasets_by_residue(r1, r2, noe, 'intersection')
    >>> print(f"Merged {len(merged)} residues, excluded {len(excluded)}")
    Merged 1 residues, excluded 1
    """
    all_residues_r1 = set(r1_data.keys())
    all_residues_r2 = set(r2_data.keys())
    all_residues_noe = set(noe_data.keys())

    if strategy == 'intersection':
        # Only keep residues present in ALL datasets
        common_residues = all_residues_r1 & all_residues_r2 & all_residues_noe
        all_residues = all_residues_r1 | all_residues_r2 | all_residues_noe
        excluded_residues = list(all_residues - common_residues)

        merged_data = {}
        for res_id in common_residues:
            merged_data[res_id] = {
                'R1': r1_data[res_id]['value'],
                'R1_err': r1_data[res_id]['error'],
                'R2': r2_data[res_id]['value'],
                'R2_err': r2_data[res_id]['error'],
                'hetNOE': noe_data[res_id]['value'],
                'hetNOE_err': noe_data[res_id]['error']
            }

    elif strategy == 'union':
        # Keep all residues, fill missing values with zeros
        all_residues = all_residues_r1 | all_residues_r2 | all_residues_noe
        excluded_residues = []

        merged_data = {}
        for res_id in all_residues:
            merged_data[res_id] = {
                'R1': r1_data.get(res_id, {}).get('value', 0.0),
                'R1_err': r1_data.get(res_id, {}).get('error', 0.0),
                'R2': r2_data.get(res_id, {}).get('value', 0.0),
                'R2_err': r2_data.get(res_id, {}).get('error', 0.0),
                'hetNOE': noe_data.get(res_id, {}).get('value', 0.0),
                'hetNOE_err': noe_data.get(res_id, {}).get('error', 0.0)
            }

    else:
        raise ValueError(f"Unknown merging strategy: {strategy}. Use 'intersection' or 'union'")

    # Sort excluded residues numerically (extract numeric part for sorting)
    def extract_numeric(res_id):
        """Extract numeric part from residue ID for sorting (e.g., '3.LYS' -> 3, '10' -> 10)"""
        import re
        match = re.search(r'\d+', str(res_id))
        return int(match.group()) if match else 0

    return merged_data, sorted(excluded_residues, key=extract_numeric)


def validate_relaxation_rates(merged_data: Dict[str, Dict],
                              warn_threshold: float = 0.5) -> Dict[str, List[str]]:
    """
    Validate physical constraints on relaxation rates.

    Checks:
    1. R1, R2 > 0 (rates must be positive)
    2. R2 >= R1 (physical constraint for isotropic tumbling)
    3. 0 < hetNOE < 1.5 (typical range for proteins)
    4. Relative errors < 50% (by default)

    Parameters
    ----------
    merged_data : dict
        Merged dataset {residue_id: {'R1': ..., 'R2': ..., 'hetNOE': ...}}
    warn_threshold : float
        Relative error threshold for warnings (default: 0.5 = 50%)

    Returns
    -------
    dict
        Dictionary of warnings by category:
        {
            'negative_rates': [list of residue IDs],
            'r2_less_than_r1': [list of residue IDs],
            'noe_out_of_range': [list of residue IDs],
            'large_errors': [list of residue IDs]
        }

    Examples
    --------
    >>> data = {'8': {'R1': 1.33, 'R1_err': 0.07, 'R2': 9.35, 'R2_err': 0.29,
    ...               'hetNOE': 0.84, 'hetNOE_err': 0.02}}
    >>> warnings = validate_relaxation_rates(data)
    >>> if warnings['r2_less_than_r1']:
    ...     print(f"Unphysical R2 < R1 for residues: {warnings['r2_less_than_r1']}")
    """
    warnings_dict = {
        'negative_rates': [],
        'r2_less_than_r1': [],
        'noe_out_of_range': [],
        'large_errors': [],
        'zero_values': []
    }

    for res_id, data in merged_data.items():
        R1 = data['R1']
        R2 = data['R2']
        hetNOE = data['hetNOE']
        R1_err = data['R1_err']
        R2_err = data['R2_err']
        hetNOE_err = data['hetNOE_err']

        # Check for zero values (will be skipped by analysis)
        if R1 == 0 or R2 == 0 or hetNOE == 0:
            warnings_dict['zero_values'].append(res_id)
            continue  # Skip other checks for zero values

        # Check 1: Negative rates
        if R1 < 0 or R2 < 0:
            warnings_dict['negative_rates'].append(res_id)

        # Check 2: R2 < R1 (unphysical for isotropic tumbling)
        if R2 < R1:
            warnings_dict['r2_less_than_r1'].append(res_id)

        # Check 3: hetNOE out of typical range
        if hetNOE < 0 or hetNOE > 1.5:
            warnings_dict['noe_out_of_range'].append(res_id)

        # Check 4: Large relative errors
        rel_err_r1 = R1_err / R1 if R1 != 0 else float('inf')
        rel_err_r2 = R2_err / R2 if R2 != 0 else float('inf')
        rel_err_noe = hetNOE_err / hetNOE if hetNOE != 0 else float('inf')

        if (rel_err_r1 > warn_threshold or
            rel_err_r2 > warn_threshold or
            rel_err_noe > warn_threshold):
            warnings_dict['large_errors'].append(res_id)

    return warnings_dict


def validate_field_consistency(field1_data: Dict[str, Dict],
                               field2_data: Dict[str, Dict]) -> Tuple[Set[str], Set[str], Set[str]]:
    """
    Validate consistency between two field datasets.

    Parameters
    ----------
    field1_data : dict
        First field dataset {residue_id: {...}}
    field2_data : dict
        Second field dataset {residue_id: {...}}

    Returns
    -------
    tuple
        (common_residues, only_in_field1, only_in_field2)

    Examples
    --------
    >>> field1 = {'8': {}, '10': {}, '11': {}}
    >>> field2 = {'8': {}, '10': {}, '12': {}}
    >>> common, f1_only, f2_only = validate_field_consistency(field1, field2)
    >>> print(f"{len(common)} common, {len(f1_only)} only in field 1, {len(f2_only)} only in field 2")
    2 common, 1 only in field 1, 1 only in field 2
    """
    residues_field1 = set(field1_data.keys())
    residues_field2 = set(field2_data.keys())

    common_residues = residues_field1 & residues_field2
    only_in_field1 = residues_field1 - residues_field2
    only_in_field2 = residues_field2 - residues_field1

    return common_residues, only_in_field1, only_in_field2


def check_minimum_residue_count(merged_data: Dict[str, Dict],
                                minimum: int = 10) -> bool:
    """
    Check if dataset has sufficient residues for meaningful analysis.

    Parameters
    ----------
    merged_data : dict
        Merged dataset
    minimum : int
        Minimum number of residues required (default: 10)

    Returns
    -------
    bool
        True if sufficient residues, raises ValueError if not

    Raises
    ------
    ValueError
        If number of residues is below minimum threshold
    """
    n_residues = len(merged_data)

    if n_residues < minimum:
        raise ValueError(
            f"Insufficient residues for analysis: {n_residues} found, "
            f"minimum {minimum} required. Check residue ID matching across datasets."
        )

    return True


def generate_validation_report(warnings_dict: Dict[str, List[str]],
                               excluded_residues: List[str]) -> str:
    """
    Generate a human-readable validation report.

    Parameters
    ----------
    warnings_dict : dict
        Dictionary of warnings from validate_relaxation_rates()
    excluded_residues : list
        List of residues excluded during merging

    Returns
    -------
    str
        Formatted validation report

    Examples
    --------
    >>> warnings = {'negative_rates': ['5'], 'r2_less_than_r1': ['12', '15'],
    ...             'noe_out_of_range': [], 'large_errors': ['8']}
    >>> excluded = ['1', '2', '3']
    >>> report = generate_validation_report(warnings, excluded)
    >>> print(report)
    """
    report_lines = []
    report_lines.append("=" * 60)
    report_lines.append("DATA VALIDATION REPORT")
    report_lines.append("=" * 60)

    # Excluded residues
    if excluded_residues:
        report_lines.append(f"\n⚠ EXCLUDED RESIDUES ({len(excluded_residues)} total):")
        report_lines.append(f"   Residues not present in all datasets (R1, R2, hetNOE):")
        report_lines.append(f"   {', '.join(excluded_residues[:20])}")  # Show first 20
        if len(excluded_residues) > 20:
            report_lines.append(f"   ... and {len(excluded_residues) - 20} more")
    else:
        report_lines.append("\n✓ All residues present in all datasets")

    # Zero values (will be skipped)
    if warnings_dict.get('zero_values'):
        zero_res = warnings_dict['zero_values']
        report_lines.append(f"\nℹ ZERO VALUES ({len(zero_res)} residues):")
        report_lines.append(f"   These residues will be skipped in analysis:")
        report_lines.append(f"   {', '.join(zero_res[:20])}")
        if len(zero_res) > 20:
            report_lines.append(f"   ... and {len(zero_res) - 20} more")

    # Negative rates (critical error)
    if warnings_dict.get('negative_rates'):
        neg_res = warnings_dict['negative_rates']
        report_lines.append(f"\n⚠ NEGATIVE RATES ({len(neg_res)} residues - CRITICAL):")
        report_lines.append(f"   Negative R1 or R2 values (unphysical):")
        report_lines.append(f"   {', '.join(neg_res)}")

    # R2 < R1 (unphysical)
    if warnings_dict.get('r2_less_than_r1'):
        unphys_res = warnings_dict['r2_less_than_r1']
        report_lines.append(f"\n⚠ R2 < R1 ({len(unphys_res)} residues - WARNING):")
        report_lines.append(f"   Unphysical for isotropic tumbling:")
        report_lines.append(f"   {', '.join(unphys_res)}")

    # hetNOE out of range
    if warnings_dict.get('noe_out_of_range'):
        noe_res = warnings_dict['noe_out_of_range']
        report_lines.append(f"\n⚠ hetNOE OUT OF RANGE ({len(noe_res)} residues):")
        report_lines.append(f"   hetNOE < 0 or > 1.5 (unusual):")
        report_lines.append(f"   {', '.join(noe_res)}")

    # Large errors
    if warnings_dict.get('large_errors'):
        err_res = warnings_dict['large_errors']
        report_lines.append(f"\nℹ LARGE ERRORS ({len(err_res)} residues):")
        report_lines.append(f"   Relative error > 50%:")
        report_lines.append(f"   {', '.join(err_res[:20])}")
        if len(err_res) > 20:
            report_lines.append(f"   ... and {len(err_res) - 20} more")

    # Summary
    report_lines.append("\n" + "=" * 60)
    total_issues = sum(len(v) for k, v in warnings_dict.items() if k != 'zero_values')
    if total_issues == 0:
        report_lines.append("✓ NO VALIDATION ISSUES DETECTED")
    else:
        report_lines.append(f"⚠ TOTAL ISSUES: {total_issues} residues flagged")
        report_lines.append("   Review warnings above before proceeding")

    report_lines.append("=" * 60)

    return "\n".join(report_lines)


def filter_invalid_residues(merged_data: Dict[str, Dict],
                            remove_negative: bool = True,
                            remove_zeros: bool = True) -> Tuple[Dict[str, Dict], List[str]]:
    """
    Filter out residues with invalid values.

    Parameters
    ----------
    merged_data : dict
        Merged dataset
    remove_negative : bool
        Remove residues with negative R1 or R2 (default: True)
    remove_zeros : bool
        Remove residues with zero values (default: True)

    Returns
    -------
    tuple
        (filtered_data, removed_residues)

    Examples
    --------
    >>> data = {'8': {'R1': 1.33, 'R2': 9.35, 'hetNOE': 0.84},
    ...         '10': {'R1': 0.0, 'R2': 0.0, 'hetNOE': 0.0}}
    >>> filtered, removed = filter_invalid_residues(data, remove_zeros=True)
    >>> print(f"Kept {len(filtered)} residues, removed {len(removed)}")
    Kept 1 residues, removed 1
    """
    filtered_data = {}
    removed_residues = []

    for res_id, data in merged_data.items():
        R1 = data['R1']
        R2 = data['R2']
        hetNOE = data['hetNOE']

        # Check for removal criteria
        should_remove = False

        if remove_negative and (R1 < 0 or R2 < 0):
            should_remove = True

        if remove_zeros and (R1 == 0 or R2 == 0 or hetNOE == 0):
            should_remove = True

        if should_remove:
            removed_residues.append(res_id)
        else:
            filtered_data[res_id] = data

    return filtered_data, removed_residues

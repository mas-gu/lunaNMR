#!/usr/bin/env python3
"""
Enhanced Multi-core T1/T2 NMR Relaxation Fitting Script with Configurable File Handling

This script performs exponential decay fitting for T1/T2 relaxation data with bootstrap error estimation.
MULTICORE RE VERSION: Enhanced with configurable input/output file management and parallel processing.

Features:
- Exponential decay fitting using lmfit
- Bootstrap error estimation
- Multiprocessing support for improved performance
- Configurable file paths in main function
- Multi-figure PDF output with customizable plots per figure
- Comprehensive results output
- Progress tracking and error handling
"""

import os
# Set thread limits before importing numpy/scipy
os.environ['OMP_NUM_THREADS'] = str(int(os.cpu_count() * 0.8))
os.environ['MKL_NUM_THREADS'] = str(int(os.cpu_count() * 0.8))

import numpy as np
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend for thread safety
import matplotlib.pyplot as plt
from lmfit import Model
import pandas as pd
import multiprocessing
from multiprocessing import Pool, cpu_count
from functools import partial
import pickle
import json
import re
from pathlib import Path
import time


def parse_delay_column(col_name):
    """
    Parse a delay column header to extract the delay value in ms.

    Handles formats:
    - Simple numbers: "300" -> 300.0
    - Duplicate markers: "300_2" -> 300.0 (sequence number ignored)
    - Fractional: "50.5" or "50.5_2" -> 50.5

    Args:
        col_name: Column header string

    Returns:
        Delay value as float, or None if not parseable
    """
    col_str = str(col_name).strip()

    # Pattern: number (with optional decimal) optionally followed by _N
    match = re.match(r'^(\d+(?:\.\d+)?)(?:_\d+)?$', col_str)
    if match:
        return float(match.group(1))

    # Fallback: try direct float conversion
    try:
        return float(col_str)
    except ValueError:
        return None


def exp_decay(x, A, t2, C=0.0):
    """Mono-exponential decay with baseline offset: I(t) = A * exp(-t/t2) + C."""
    return A * np.exp(-x / t2) + C


def bootstrap_errors(x, y, model, params, n_bootstrap=1000):
    """
    Classic residual bootstrap error estimation.

    Uses simple residual resampling with replacement - the standard,
    well-understood bootstrap approach for regression. Assumes homoscedastic
    noise (constant variance), which is appropriate for NMR thermal noise.

    Parameters
    ----------
    x : array
        Time points
    y : array
        Signal intensities
    model : lmfit.Model
        Fitting model
    params : lmfit.Parameters
        Initial parameters
    n_bootstrap : int
        Number of bootstrap iterations (default: 1000)

    Returns
    -------
    tuple : (A_error, t2_error)

    Notes
    -----
    Algorithm:
    1. Fit data to get y_fit and residuals
    2. For each bootstrap iteration:
       - Resample residuals WITH REPLACEMENT
       - Create synthetic data: y_synthetic = y_fit + resampled_residuals
       - Refit to get parameter estimates
    3. Return std of bootstrap parameter distributions
    """
    # Initial fit (unweighted)
    result = model.fit(y, params, x=x)
    y_fit = result.best_fit
    residuals = y - y_fit
    n = len(residuals)

    a_values = []
    t2_values = []
    c_values = []

    for _ in range(n_bootstrap):
        # Resample residuals WITH REPLACEMENT
        resampled_indices = np.random.randint(0, n, size=n)
        resampled_residuals = residuals[resampled_indices]

        # Create synthetic data
        y_synthetic = y_fit + resampled_residuals

        # Ensure positive values
        y_synthetic = np.maximum(y_synthetic, 1e-10)

        # Refit (unweighted)
        params_boot = model.make_params(
            A=params['A'].value, t2=params['t2'].value, C=params['C'].value
        )
        params_boot['A'].min = 0
        params_boot['t2'].min = 0

        try:
            res = model.fit(y_synthetic, params_boot, x=x)
            if res.success:
                a_values.append(res.params['A'].value)
                t2_values.append(res.params['t2'].value)
                c_values.append(res.params['C'].value)
        except Exception:
            continue

    # Return standard deviation of bootstrap distribution
    if len(a_values) < 10:
        # Not enough successful fits, fall back to covariance estimate
        a_err = result.params['A'].stderr if result.params['A'].stderr else np.nan
        t2_err = result.params['t2'].stderr if result.params['t2'].stderr else np.nan
        c_err = result.params['C'].stderr if result.params['C'].stderr else np.nan
        return a_err, t2_err, c_err

    return np.std(a_values), np.std(t2_values), np.std(c_values)


def fit_single_residue_parallel(args):
    """
    Wrapper function for parallel processing of single residue fitting.

    Model: I(t) = A * exp(-t/t2) + C

    Parameters
    ----------
    args : tuple
        (x, y, residue_name, initial_A, initial_t2, initial_C,
         n_bootstrap, error_method, idx, total)

    Returns
    -------
    dict : Fitting results with success flag
    """
    (x, y, residue_name, initial_A, initial_t2, initial_C,
     n_bootstrap, error_method, idx, total) = args

    try:
        # Progress indicator
        if idx % 10 == 0:
            print(f"Processing residue {idx+1}/{total}: {residue_name}")

        if initial_C is None:
            initial_C = float(np.min(y))

        model = Model(exp_decay)
        params = model.make_params(A=initial_A, t2=initial_t2, C=initial_C)
        params['A'].min = 0
        params['t2'].min = 0

        # Unweighted least squares fit
        result = model.fit(y, params, x=x)

        if not result.success:
            print(f"Warning: Fit failed for residue {residue_name}")
            return {
                'residue': residue_name,
                'A': np.nan,
                't2': np.nan,
                'C': np.nan,
                'A_err': np.nan,
                't2_err': np.nan,
                'C_err': np.nan,
                'success': False,
                'idx': idx
            }

        a = result.params['A'].value
        t2 = result.params['t2'].value
        c = result.params['C'].value

        # Error estimation based on selected method
        if error_method == 'bootstrap':
            a_err, t2_err, c_err = bootstrap_errors(x, y, model, params, n_bootstrap)
        else:
            # Analytical: use covariance matrix from lmfit
            a_err = result.params['A'].stderr if result.params['A'].stderr else np.nan
            t2_err = result.params['t2'].stderr if result.params['t2'].stderr else np.nan
            c_err = result.params['C'].stderr if result.params['C'].stderr else np.nan

        return {
            'residue': residue_name,
            'A': a,
            't2': t2,
            'C': c,
            'A_err': a_err,
            't2_err': t2_err,
            'C_err': c_err,
            'x': x,
            'y': y,
            'result': result,
            'success': True,
            'idx': idx
        }

    except Exception as e:
        print(f"Error processing residue {residue_name}: {str(e)}")
        return {
            'residue': residue_name,
            'A': np.nan,
            't2': np.nan,
            'C': np.nan,
            'A_err': np.nan,
            't2_err': np.nan,
            'C_err': np.nan,
            'success': False,
            'error': str(e),
            'idx': idx
        }


def create_plots(results_list, output_prefix, n_plots_per_figure=20,
                 experiment_type="T1", time_units="ms", signal_units="Intensity"):
    """
    Create multi-figure PDF plots of fitting results

    Parameters:
    -----------
    results_list : list
        List of fitting results dictionaries
    output_prefix : str
        Prefix for output PDF files
    n_plots_per_figure : int
        Number of plots per figure
    experiment_type : str
        Type of experiment (T1, T2, etc.)
    time_units : str
        Units for time axis
    signal_units : str
        Units for signal axis
    """
    # Filter successful fits for plotting
    successful_results = [r for r in results_list if r.get('success', False)]

    if not successful_results:
        print("Warning: No successful fits to plot")
        return

    n_datasets = len(successful_results)
    n_figures = int(np.ceil(n_datasets / n_plots_per_figure))

    print(f"Creating {n_figures} figure(s) with {n_datasets} successful fits...")

    for fig_idx in range(n_figures):
        fig, axes = plt.subplots(5, 4, figsize=(20, 25), sharex=True)
        axes = axes.flatten()
        start_idx = fig_idx * n_plots_per_figure
        end_idx = min((fig_idx + 1) * n_plots_per_figure, n_datasets)

        for idx, result_idx in enumerate(range(start_idx, end_idx)):
            if result_idx >= len(successful_results):
                break

            result = successful_results[result_idx]
            x = result['x']
            y = result['y']
            residue = result['residue']
            a = result['A']
            t2 = result['t2']
            t2_err = result['t2_err']
            c = result.get('C', 0.0)

            # Generate fit curve
            x_fit = np.linspace(0, max(x)*1.2, 50)
            y_fit = a * np.exp(-x_fit / t2) + c

            ax = axes[idx]
            ax.plot(x, y, 'ko', lw=2, ms=8, label="Data")
            ax.plot(x_fit, y_fit, 'r-', linewidth=2, label="Fit")

            # Add fitting results text
            textstr = f"{experiment_type} = {t2:.2f} ± {t2_err:.2f} {time_units}"
            ax.text(0.05, 0.95, textstr, transform=ax.transAxes,
                    fontsize=10, verticalalignment='top', horizontalalignment='left',
                    bbox=dict(boxstyle="round,pad=0.3", facecolor='white', alpha=0.8))

            ax.set_title(f"Residue: {residue}")
            ax.set_xlabel(f"Time ({time_units})")
            ax.set_ylabel(f"Signal ({signal_units})")
            ax.legend()
            ax.set_xlim(0, max(x)*1.1)
            ax.set_ylim(min(y)*0.9, max(y)*1.1)

        # Hide unused subplots
        for idx in range(end_idx - start_idx, len(axes)):
            axes[idx].set_visible(False)

        plt.tight_layout()
        fig.savefig(f"{output_prefix}_fit_results_fig{fig_idx + 1}.pdf", format='pdf', dpi=300)
        plt.close(fig)


def save_results(results_list, output_file, experiment_type="T1"):
    """
    Save fitting results to text file

    Parameters:
    -----------
    results_list : list
        List of fitting results dictionaries
    output_file : str
        Output filename
    experiment_type : str
        Type of experiment for headers
    """
    with open(output_file, "w") as f:
        f.write(
            f"Residue\tA\t{experiment_type}\tC\tA_err\t{experiment_type}_err\tC_err\tSuccess\n"
        )
        for result in results_list:
            success_flag = "Yes" if result.get('success', False) else "No"
            if result.get('success', False):
                f.write(
                    f"{result['residue']}\t{result['A']:.6e}\t{result['t2']:.6e}\t"
                    f"{result['C']:.6e}\t{result['A_err']:.6e}\t"
                    f"{result['t2_err']:.6e}\t{result['C_err']:.6e}\t{success_flag}\n"
                )
            else:
                f.write(
                    f"{result['residue']}\tNaN\tNaN\tNaN\tNaN\tNaN\tNaN\t{success_flag}\n"
                )


def save_fit_data_json(results_list, output_file, experiment_type, time_units,
                       signal_units, n_bootstrap, field_freq):
    """
    Save complete fitting data as JSON for interactive visualization

    Parameters:
    -----------
    results_list : list
        List of fitting results dictionaries
    output_file : str
        Output JSON filename
    experiment_type : str
        Type of experiment (T1, T2)
    time_units : str
        Units for time axis
    signal_units : str
        Units for signal axis
    n_bootstrap : int
        Number of bootstrap iterations used
    field_freq : float
        Magnetic field frequency in MHz
    """
    # Get time points from first successful result (same for all residues)
    successful_results = [r for r in results_list if r.get('success', False)]
    if not successful_results:
        print(f"Warning: No successful fits to save to JSON")
        return

    time_points = successful_results[0]['x'].tolist()

    # Build fit curve points (dense sampling for smooth plotting)
    max_time = max(time_points)
    fit_time_dense = np.linspace(0, max_time * 1.2, 100)

    fits_data = []
    for result in successful_results:
        c = result.get('C', 0.0)
        c_err = result.get('C_err', np.nan)
        # Calculate fit curve using fitted parameters
        fit_intensity = result['A'] * np.exp(-fit_time_dense / result['t2']) + c

        fits_data.append({
            'residue': str(result['residue']),
            'A': float(result['A']),
            't2': float(result['t2']),
            'C': float(c),
            'A_err': float(result['A_err']),
            't2_err': float(result['t2_err']),
            'C_err': float(c_err),
            'intensities': [float(val) for val in result['y']],
            'fit_curve': {
                'time': [float(t) for t in fit_time_dense],
                'intensity': [float(i) for i in fit_intensity]
            }
        })

    output_data = {
        'metadata': {
            'experiment_type': experiment_type,
            'field_freq': float(field_freq),
            'time_units': time_units,
            'signal_units': signal_units,
            'n_bootstrap': n_bootstrap,
            'n_residues': len(fits_data),
            'time_points': [float(t) for t in time_points]
        },
        'fits': fits_data
    }

    # Ensure output directory exists
    os.makedirs(os.path.dirname(output_file), exist_ok=True)

    with open(output_file, 'w') as f:
        json.dump(output_data, f, indent=2)

    print(f"Saved JSON fit data to: {output_file}")


def main():
    """
    Main function with configurable file handling and multiprocessing
    """

    # ========================
    # CONFIGURATION SECTION
    # ========================

    # Input file configuration
    input_csv_file = "T2_T6D_data.csv"  # Input CSV file path

    # Output file configuration
    output_prefix = "700_T6D_T2"  # Prefix for output files
    results_txt_file = "700_T6D_T2_fit_results.txt"  # Results text file

    # Experiment configuration
    experiment_type = "T2"  # T1, T2, etc. (for labels and headers)
    time_units = "ms"  # Time axis units
    signal_units = "Intensity"  # Signal axis units

    # Fitting parameters
    initial_A = 5  # Initial amplitude estimate
    initial_t2 = 100  # Initial time constant estimate
    n_bootstrap = 1000  # Number of bootstrap iterations for error estimation
    error_method = 'analytical'  # 'analytical' (fast) or 'bootstrap' (robust)

    # Plot configuration
    n_plots_per_figure = 20  # Number of plots per PDF page

    # Multiprocessing configuration
    n_processes = None  # Use None for automatic detection (80% of cores)
    # n_processes = 4  # Or specify a specific number of processes

    # Optional: Time scaling factor (uncomment if needed)
    # time_scaling_factor = 1000  # Multiply time values by this factor

    # ========================
    # END CONFIGURATION
    # ========================

    # Validate input file exists
    if not os.path.exists(input_csv_file):
        raise FileNotFoundError(f"Input file not found: {input_csv_file}")

    # Determine number of processes
    if n_processes is None:
        n_processes = max(1, int(cpu_count() * 0.8))

    print(f"Starting {experiment_type} relaxation fitting analysis...")
    print(f"Input file: {input_csv_file}")
    print(f"Output prefix: {output_prefix}")
    print(f"Error method: {error_method}")
    if error_method == 'bootstrap':
        print(f"Bootstrap iterations: {n_bootstrap}")
    print(f"Number of processes: {n_processes}")

    # Load CSV using pandas to handle mixed data types
    print("Loading data...")
    raw_df = pd.read_csv(input_csv_file, header=None)

    # Detect format: LunaNMR Fit Series vs simple DynamiXs format
    header_row = raw_df.iloc[0].astype(str).tolist()
    lunaNMR_columns = {'Peak_Number', 'Assignment', 'Reference_X', 'Reference_Y'}
    detected_lunaNMR_cols = [col for col in header_row if col in lunaNMR_columns]

    if detected_lunaNMR_cols:
        print(f"Detected LunaNMR Fit Series format (columns: {detected_lunaNMR_cols})")
        delay_start_idx = 0
        for i, col in enumerate(header_row):
            if col not in lunaNMR_columns:
                try:
                    float(col)
                    delay_start_idx = i
                    break
                except ValueError:
                    continue
        assignment_idx = header_row.index('Assignment') if 'Assignment' in header_row else 1
        print(f"  Using 'Assignment' column (index {assignment_idx}) for residue names")
        print(f"  Delay columns start at index {delay_start_idx}")
        residue_names = raw_df.iloc[1:, assignment_idx].to_numpy()
        # Parse delay columns - handles duplicates like "300_2" -> 300.0
        delay_headers = raw_df.iloc[0, delay_start_idx:].tolist()
        x = np.array([parse_delay_column(col) for col in delay_headers])
        y_data = raw_df.iloc[1:, delay_start_idx:].astype(float).to_numpy()
    else:
        print("Detected simple DynamiXs format")
        residue_names = raw_df.iloc[1:, 0].to_numpy()  # Column 1 (residue names), skip header
        # Parse delay columns - handles duplicates like "300_2" -> 300.0
        delay_headers = raw_df.iloc[0, 1:].tolist()
        x = np.array([parse_delay_column(col) for col in delay_headers])
        y_data = raw_df.iloc[1:, 1:].astype(float).to_numpy()  # Values to fit (rows 2+, columns 2+)

    # Filter out dummy residues (case-insensitive)
    valid_mask = np.array([not str(name).lower().startswith('dummy') for name in residue_names])
    n_dummies = np.sum(~valid_mask)
    if n_dummies > 0:
        print(f"Excluding {n_dummies} dummy residue(s) from fitting")
        residue_names = residue_names[valid_mask]
        y_data = y_data[valid_mask]

    print(f"Loaded {len(residue_names)} residues with {len(x)} time points")
    print(f"Time range: {x.min():.3f} to {x.max():.3f} {time_units}")

    # Prepare arguments for parallel processing
    print("Preparing parallel processing...")
    args_list = []
    for idx, residue in enumerate(residue_names):
        y = y_data[idx, :]
        args_list.append((x, y, residue, initial_A, initial_t2, n_bootstrap, error_method, idx, len(residue_names)))

    # Parallel fitting
    print(f"Starting parallel fitting with {n_processes} processes...")
    start_time = time.time()

    with Pool(processes=n_processes) as pool:
        results_list = pool.map(fit_single_residue_parallel, args_list)

    end_time = time.time()
    elapsed_time = end_time - start_time

    # Sort results by original index to maintain order
    results_list = sorted(results_list, key=lambda x: x['idx'])

    # Remove idx from results for consistency
    for result in results_list:
        result.pop('idx', None)

    print(f"Completed fitting in {elapsed_time:.2f} seconds")

    # Statistics
    successful_fits = [r for r in results_list if r.get('success', False)]
    failed_fits = [r for r in results_list if not r.get('success', False)]

    print(f"Successful fits: {len(successful_fits)}/{len(results_list)}")
    if failed_fits:
        print(f"Failed fits: {len(failed_fits)}")
        print("Failed residues:", [r['residue'] for r in failed_fits])

    # Create plots
    if successful_fits:
        print("Generating plots...")
        create_plots(results_list, output_prefix, n_plots_per_figure,
                     experiment_type, time_units, signal_units)

    # Save results
    print(f"Saving results to {results_txt_file}...")
    save_results(results_list, results_txt_file, experiment_type)

    # Summary statistics
    if successful_fits:
        t2_values = [r['t2'] for r in successful_fits]
        t2_errors = [r['t2_err'] for r in successful_fits]

        print(f"\n{experiment_type} Analysis Summary:")
        print(f"  Number of residues fitted: {len(successful_fits)}/{len(results_list)}")
        print(f"  {experiment_type} range: {min(t2_values):.2f} to {max(t2_values):.2f} {time_units}")
        print(f"  Mean {experiment_type}: {np.mean(t2_values):.2f} ± {np.std(t2_values):.2f} {time_units}")
        print(f"  Mean fitting error: {np.mean(t2_errors):.2f} {time_units}")
        print(f"  Processing time: {elapsed_time:.2f} seconds")
        print(f"  Time per residue: {elapsed_time/len(results_list):.2f} seconds")
        print(f"  Results saved to: {results_txt_file}")
        print(f"  Plots saved with prefix: {output_prefix}")

    print("Analysis completed successfully!")


def run_analysis_with_params(params, progress_callback=None):
    """
    Run multicore analysis with parameters provided by GUI

    Parameters
    ----------
    params : dict
        Dictionary containing all analysis parameters
    progress_callback : callable, optional
        Callback function for progress updates.
        Signature: progress_callback(current, total, residue_name, message)
        - current: Current residue number (1-indexed)
        - total: Total number of residues
        - residue_name: Name of the residue being/just fitted
        - message: Status message

    Returns
    -------
    dict : Results summary
    """

    # Input file configuration
    input_csv_file = params['input_csv_file']

    # Output file configuration
    output_prefix = params['output_prefix']
    results_txt_file = params['results_txt_file']

    # Experiment configuration
    experiment_type = params['experiment_type']
    time_units = params.get('time_units', 'ms')
    signal_units = params.get('signal_units', 'Intensity')

    # Fitting parameters
    initial_A = params.get('initial_A', 5)
    initial_t2 = params.get('initial_t2', 100)
    n_bootstrap = params.get('n_bootstrap', 1000)
    error_method = params.get('error_method', 'analytical')

    # Plot configuration
    n_plots_per_figure = params.get('n_plots_per_figure', 20)

    # Multiprocessing configuration
    n_processes = max(1, int(cpu_count() * 0.8))

    # Validate input file exists
    if not os.path.exists(input_csv_file):
        raise FileNotFoundError(f"Input file not found: {input_csv_file}")

    print(f"Starting {experiment_type} relaxation fitting analysis (multicore)...")
    print(f"Input file: {input_csv_file}")
    print(f"Output prefix: {output_prefix}")
    print(f"Error method: {error_method}")
    if error_method == 'bootstrap':
        print(f"Bootstrap iterations: {n_bootstrap}")
    print(f"Number of processes: {n_processes}")

    # Load CSV using pandas to handle mixed data types
    print("Loading data...")
    raw_df = pd.read_csv(input_csv_file, header=None)

    # Detect format: LunaNMR Fit Series vs simple DynamiXs format
    # LunaNMR format has columns: Peak_Number, Assignment, Reference_X, Reference_Y, delay1, delay2, ...
    # Simple format has columns: Residue, delay1, delay2, ...
    header_row = raw_df.iloc[0].astype(str).tolist()

    # Check for LunaNMR Fit Series format
    lunaNMR_columns = {'Peak_Number', 'Assignment', 'Reference_X', 'Reference_Y'}
    detected_lunaNMR_cols = [col for col in header_row if col in lunaNMR_columns]

    if detected_lunaNMR_cols:
        print(f"Detected LunaNMR Fit Series format (columns: {detected_lunaNMR_cols})")

        # Find where the numeric delay columns start
        delay_start_idx = 0
        for i, col in enumerate(header_row):
            if col not in lunaNMR_columns:
                # Check if this looks like a numeric delay value
                try:
                    float(col)
                    delay_start_idx = i
                    break
                except ValueError:
                    continue

        # Find the Assignment column index for residue names
        assignment_idx = header_row.index('Assignment') if 'Assignment' in header_row else 1

        print(f"  Using 'Assignment' column (index {assignment_idx}) for residue names")
        print(f"  Delay columns start at index {delay_start_idx}")

        residue_names = raw_df.iloc[1:, assignment_idx].to_numpy()  # Assignment column, skip header
        # Parse delay columns - handles duplicates like "300_2" -> 300.0
        delay_headers = raw_df.iloc[0, delay_start_idx:].tolist()
        x = np.array([parse_delay_column(col) for col in delay_headers])
        y_data = raw_df.iloc[1:, delay_start_idx:].astype(float).to_numpy()  # Intensity data
    else:
        # Simple DynamiXs format: first column is residue, rest are delays
        print("Detected simple DynamiXs format")
        residue_names = raw_df.iloc[1:, 0].to_numpy()  # Column 1 (residue names), skip header
        # Parse delay columns - handles duplicates like "300_2" -> 300.0
        delay_headers = raw_df.iloc[0, 1:].tolist()
        x = np.array([parse_delay_column(col) for col in delay_headers])
        y_data = raw_df.iloc[1:, 1:].astype(float).to_numpy()  # Values to fit (rows 2+, columns 2+)

    # Filter out dummy residues (case-insensitive)
    valid_mask = np.array([not str(name).lower().startswith('dummy') for name in residue_names])
    n_dummies = np.sum(~valid_mask)
    if n_dummies > 0:
        print(f"Excluding {n_dummies} dummy residue(s) from fitting")
        residue_names = residue_names[valid_mask]
        y_data = y_data[valid_mask]

    total_residues = len(residue_names)
    print(f"Loaded {total_residues} residues with {len(x)} time points")
    print(f"Time range: {x.min():.3f} to {x.max():.3f} {time_units}")

    # Report data loaded via callback
    if progress_callback:
        progress_callback(0, total_residues, "", f"Loaded {total_residues} residues")

    # Prepare arguments for parallel processing
    print("Preparing parallel processing...")
    args_list = []
    for idx, residue in enumerate(residue_names):
        y = y_data[idx, :]
        # initial_C=None lets fit_single_residue_parallel use min(y) per-residue.
        args_list.append((x, y, residue, initial_A, initial_t2, None,
                          n_bootstrap, error_method, idx, total_residues))

    # Parallel fitting with progress reporting
    print(f"Starting parallel fitting with {n_processes} processes...")
    start_time = time.time()

    with Pool(processes=n_processes) as pool:
        # Use imap_unordered for per-residue progress tracking
        results_list = []
        completed = 0

        for result in pool.imap_unordered(fit_single_residue_parallel, args_list):
            results_list.append(result)
            completed += 1
            residue_name = result.get('residue', f'residue_{completed}')

            # Report progress via callback
            if progress_callback:
                progress_callback(completed, total_residues, residue_name,
                                  f"Fitted {residue_name} ({completed}/{total_residues})")

            # Also print to console for non-GUI usage (every 10 or at end)
            if completed % 10 == 0 or completed == total_residues:
                print(f"Progress: {completed}/{total_residues} residues fitted")

    end_time = time.time()
    elapsed_time = end_time - start_time

    # Sort results by original index to maintain order
    results_list = sorted(results_list, key=lambda x: x['idx'])

    # Remove idx from results for consistency
    for result in results_list:
        result.pop('idx', None)

    print(f"Completed fitting in {elapsed_time:.2f} seconds")

    # Statistics
    successful_fits = [r for r in results_list if r.get('success', False)]
    failed_fits = [r for r in results_list if not r.get('success', False)]

    print(f"Successful fits: {len(successful_fits)}/{len(results_list)}")
    if failed_fits:
        print(f"Failed fits: {len(failed_fits)}")
        print("Failed residues:", [r['residue'] for r in failed_fits])

    # Create plots
    if successful_fits:
        print("Generating plots...")
        if progress_callback:
            progress_callback(total_residues, total_residues, "", "Generating plots...")
        create_plots(results_list, output_prefix, n_plots_per_figure,
                     experiment_type, time_units, signal_units)

    # Save results
    print(f"Saving results to {results_txt_file}...")
    if progress_callback:
        progress_callback(total_residues, total_residues, "", "Saving results...")
    save_results(results_list, results_txt_file, experiment_type)

    # Save JSON data for interactive visualization
    json_folder = params.get('json_folder')
    json_file = None
    if json_folder:
        # Create json folder if it doesn't exist
        os.makedirs(json_folder, exist_ok=True)

        # Construct JSON filename: field{N}_{T1|T2}_fit_data.json
        field_name = params.get('field_name', 'field1')
        json_filename = f"{field_name}_{experiment_type}_fit_data.json"
        json_file = os.path.join(json_folder, json_filename)

        print(f"Saving JSON fit data to {json_file}...")
        save_fit_data_json(
            results_list=results_list,
            output_file=json_file,
            experiment_type=experiment_type,
            time_units=time_units,
            signal_units=signal_units,
            n_bootstrap=n_bootstrap,
            field_freq=params.get('field_freq', 600.0)
        )

    # Summary statistics
    if successful_fits:
        t2_values = [r['t2'] for r in successful_fits]
        t2_errors = [r['t2_err'] for r in successful_fits]

        print(f"\n{experiment_type} Analysis Summary:")
        print(f"  Number of residues fitted: {len(successful_fits)}/{len(results_list)}")
        print(f"  {experiment_type} range: {min(t2_values):.2f} to {max(t2_values):.2f} {time_units}")
        print(f"  Mean {experiment_type}: {np.mean(t2_values):.2f} ± {np.std(t2_values):.2f} {time_units}")
        print(f"  Mean fitting error: {np.mean(t2_errors):.2f} {time_units}")
        print(f"  Processing time: {elapsed_time:.2f} seconds")
        print(f"  Time per residue: {elapsed_time/len(results_list):.2f} seconds")
        print(f"  Results saved to: {results_txt_file}")
        print(f"  Plots saved with prefix: {output_prefix}")
        if json_file:
            print(f"  JSON data saved to: {json_file}")

    print("Multicore analysis completed successfully!")

    # Return results summary
    return {
        'n_fitted': len(successful_fits),
        'n_total': len(results_list),
        'results_file': results_txt_file,
        'json_file': json_file,
        'plots_prefix': output_prefix,
        't2_range': (min(t2_values), max(t2_values)) if successful_fits else (0, 0),
        'mean_t2': np.mean(t2_values) if successful_fits else 0,
        'std_t2': np.std(t2_values) if successful_fits else 0,
        'mean_error': np.mean(t2_errors) if successful_fits else 0,
        'n_cores_used': n_processes,
        'processing_time': elapsed_time,
        'failed_residues': [r['residue'] for r in failed_fits]
    }


if __name__ == "__main__":
    main()

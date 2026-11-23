#!/usr/bin/env python3
"""
Enhanced Multi-core T1/T2 NMR Relaxation Fitting Script with Configurable File Handling

This script performs exponential decay fitting for T1/T2 relaxation data with bootstrap error estimation.
MULTICORE RE VERSION: Enhanced with configurable input/output file management and parallel processing.

Features:
- Exponential decay fitting using lmfit
- Bootstrap error estimation
- Multiprocessing support for improved performance
- Multi-figure PDF output with customizable plots per figure
- Comprehensive results output
- Progress tracking and error handling
"""

import os

import numpy as np
# One definition of the fit and its conventions, shared with the single-core fitter.
try:
    from fit_Tx_NMRRE import (DEGENERATE_T2_OVER_TMAX, FIT_BASELINE_BY_DEFAULT,
                              fit_single_residue, read_relaxation_matrix, save_results,
                              unfitted_residue)
    from fit_Tx_NMRRE import create_plots as plot_fit_results
    from fit_Tx_NMRRE import save_fit_data_json as save_fits_as_json
except ImportError:  # imported as dynamiXs_T1_T2.fitmulti__Tx_NMRRE (parent dir on path)
    from dynamiXs_T1_T2.fit_Tx_NMRRE import (DEGENERATE_T2_OVER_TMAX, FIT_BASELINE_BY_DEFAULT,
                                             fit_single_residue, read_relaxation_matrix,
                                             save_results, unfitted_residue)
    from dynamiXs_T1_T2.fit_Tx_NMRRE import create_plots as plot_fit_results
    from dynamiXs_T1_T2.fit_Tx_NMRRE import save_fit_data_json as save_fits_as_json
import multiprocessing
from multiprocessing import cpu_count
import time


def fit_single_residue_parallel(args):
    """
    fit_single_residue for the worker pool: its arguments in one tuple, the residue's
    position added to its result, and an exception returned as a failed fit.

    The baseline convention travels IN `args`. A module-level global cannot work here:
    macOS spawns fresh workers that re-import this module, so each would read the value
    as declared in the source and silently ignore whatever the caller chose.

    Parameters
    ----------
    args : tuple
        (x, y, residue_name, initial_A, initial_t2, initial_C,
         n_bootstrap, error_method, idx, total[, fit_baseline])

    Returns
    -------
    dict : Fitting results with success flag
    """
    (x, y, residue_name, initial_A, initial_t2, initial_C,
     n_bootstrap, error_method, idx, total) = args[:10]
    fit_baseline = bool(args[10]) if len(args) > 10 else FIT_BASELINE_BY_DEFAULT

    try:
        # Progress indicator
        if idx % 10 == 0:
            print(f"Processing residue {idx+1}/{total}: {residue_name}")

        result = fit_single_residue(x, y, residue_name, initial_A, initial_t2, initial_C,
                                    n_bootstrap, error_method, fit_baseline)

    except Exception as e:
        print(f"Error processing residue {residue_name}: {str(e)}")
        result = {**unfitted_residue(residue_name, x, y, fit_baseline), 'error': str(e)}

    result['idx'] = idx
    return result


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
    successful_results = [r for r in results_list if r.get('success', False)]

    if not successful_results:
        print("Warning: No successful fits to plot")
        return

    n_datasets = len(successful_results)
    n_figures = int(np.ceil(n_datasets / n_plots_per_figure))

    print(f"Creating {n_figures} figure(s) with {n_datasets} successful fits...")
    plot_fit_results(successful_results, output_prefix, n_plots_per_figure,
                     experiment_type, time_units, signal_units)


def save_fit_data_json(results_list, output_file, experiment_type, time_units,
                       signal_units, n_bootstrap, field_freq, results_file=None):
    """
    Save the successful fits as JSON for interactive visualization

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
    results_file : str, optional
        The run's results file, recorded relative to the JSON for a refit to update
    """
    successful_results = [r for r in results_list if r.get('success', False)]
    if not successful_results:
        print("Warning: No successful fits to save to JSON")
        return

    save_fits_as_json(successful_results, output_file, experiment_type, time_units,
                      signal_units, n_bootstrap, field_freq, results_file=results_file)


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
    initial_A = params.get('initial_A')
    initial_t2 = params.get('initial_t2')
    n_bootstrap = params.get('n_bootstrap', 1000)
    error_method = params.get('error_method', 'analytical')

    n_plots_per_figure = params.get('n_plots_per_figure', 20)

    # Multiprocessing configuration
    n_processes = max(1, int(cpu_count() * 0.8))

    if not os.path.exists(input_csv_file):
        raise FileNotFoundError(f"Input file not found: {input_csv_file}")

    print(f"Starting {experiment_type} relaxation fitting analysis (multicore)...")
    print(f"Input file: {input_csv_file}")
    print(f"Output prefix: {output_prefix}")
    print(f"Error method: {error_method}")
    if error_method == 'bootstrap':
        print(f"Bootstrap iterations: {n_bootstrap}")
    print(f"Number of processes: {n_processes}")

    residue_names, x, y_data, dropped_columns = read_relaxation_matrix(input_csv_file)

    total_residues = len(residue_names)
    print(f"Loaded {total_residues} residues with {len(x)} time points")
    print(f"Time range: {x.min():.3f} to {x.max():.3f} {time_units}")

    # Report data loaded via callback
    if progress_callback:
        progress_callback(0, total_residues, "", f"Loaded {total_residues} residues")

    print("Preparing parallel processing...")
    # Carried per call rather than read from the module global, so a caller
    # can choose and spawned workers cannot silently disagree.
    fit_baseline = bool(params.get('fit_baseline', FIT_BASELINE_BY_DEFAULT)) \
        if isinstance(params, dict) else FIT_BASELINE_BY_DEFAULT
    args_list = []
    for idx, residue in enumerate(residue_names):
        y = y_data[idx, :]
        # initial_C=None lets fit_single_residue_parallel use min(y) per-residue.
        args_list.append((x, y, residue, initial_A, initial_t2, None,
                          n_bootstrap, error_method, idx, total_residues,
                          fit_baseline))

    # Parallel fitting with progress reporting
    print(f"Starting parallel fitting with {n_processes} processes...")
    start_time = time.time()

    with multiprocessing.get_context('spawn').Pool(processes=n_processes) as pool:
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

    if successful_fits:
        print("Generating plots...")
        if progress_callback:
            progress_callback(total_residues, total_residues, "", "Generating plots...")
        create_plots(results_list, output_prefix, n_plots_per_figure,
                     experiment_type, time_units, signal_units)

    print(f"Saving results to {results_txt_file}...")
    if progress_callback:
        progress_callback(total_residues, total_residues, "", "Saving results...")
    save_results(results_list, results_txt_file, experiment_type)

    json_folder = params.get('json_folder')
    json_file = None
    if json_folder:
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
            field_freq=params.get('field_freq', 600.0),
            results_file=results_txt_file
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

    return {
        'n_fitted': len(successful_fits),
        'dropped_columns': dropped_columns,
        'n_excluded': len(failed_fits),
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
        'failed_residues': [r['residue'] for r in failed_fits],
        # Fit provenance, so a caller can see the convention and the sampling quality
        # without opening the per-residue files.
        'baseline_fixed': not fit_baseline,
        'n_window_marginal': sum(1 for r in successful_fits if r.get('window_marginal')),
    }

#!/usr/bin/env python3
"""
Multi-core T1/T2 NMR Relaxation Fitting

The exponential decay fit of fit_Tx_NMRRE - the same residue fit, results file, plots
and fit JSON - with the residues fitted over a pool of worker processes.
"""

from multiprocessing import cpu_count

# One definition of the fit, its conventions and the analysis run, shared with the single-core fitter.
try:
    from fit_Tx_NMRRE import DEGENERATE_T2_OVER_TMAX
    from fit_Tx_NMRRE import run_analysis_with_params as run_analysis
except ImportError:  # imported as dynamiXs_T1_T2.fitmulti__Tx_NMRRE (parent dir on path)
    from dynamiXs_T1_T2.fit_Tx_NMRRE import DEGENERATE_T2_OVER_TMAX
    from dynamiXs_T1_T2.fit_Tx_NMRRE import run_analysis_with_params as run_analysis


def run_analysis_with_params(params, progress_callback=None):
    """
    Run the analysis with parameters provided by GUI, fitting over 80% of the cores

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
    return run_analysis(params, progress_callback,
                        n_processes=max(1, int(cpu_count() * 0.8)))

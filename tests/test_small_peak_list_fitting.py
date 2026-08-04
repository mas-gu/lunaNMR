# ABOUTME: Regression tests for fitting peak lists too small to be worth parallelising.
# ABOUTME: A one- or two-peak list must still come back fitted, not silently empty.

import numpy as np
import pandas as pd

from lunaNMR.core.core_integrator import EnhancedVoigtIntegrator

PEAKS = (('Peak_101', 9.56, 120.3, 1.0),
         ('Peak_106', 9.44, 120.06, 0.8))


def _synthetic_integrator():
    """An integrator holding a two-peak 15N-HSQC-like spectrum, built in memory."""
    ppm_x = np.linspace(10.0, 6.0, 512)
    ppm_y = np.linspace(130.0, 100.0, 256)
    X, Y = np.meshgrid(ppm_x, ppm_y)

    data = np.random.default_rng(0).normal(0, 0.01, X.shape)
    for _, x0, y0, amplitude in PEAKS:
        data += amplitude * np.exp(-((X - x0) ** 2 / (2 * 0.02 ** 2) +
                                     (Y - y0) ** 2 / (2 * 0.25 ** 2)))

    integrator = EnhancedVoigtIntegrator()
    integrator.nmr_data = data
    integrator.ppm_x_axis = ppm_x
    integrator.ppm_y_axis = ppm_y
    integrator._estimate_noise_level()
    return integrator


def _peak_list(count):
    rows = [{'Assignment': name, 'Position_X': x, 'Position_Y': y}
            for name, x, y, _ in PEAKS[:count]]
    return pd.DataFrame(rows)


def test_two_peaks_are_fitted_although_parallel_is_skipped():
    # Below the parallel threshold the fitter takes its sequential route; a series
    # run of two peaks was coming back with zero fits for every spectrum.
    integrator = _synthetic_integrator()
    peaks = _peak_list(2)

    results, _ = integrator.enhanced_fitter.enhanced_peak_fitting_parallel(
        peaks, use_parallel=True)

    assert len(results) == 2
    assert all(result['success'] for result in results)
    assert [result['assignment'] for result in results] == ['Peak_101', 'Peak_106']


def test_single_peak_is_fitted():
    integrator = _synthetic_integrator()

    result, _ = integrator.enhanced_fitter.enhanced_peak_fitting_parallel(
        _peak_list(1), use_parallel=True)

    assert result['success']
    assert result['assignment'] == 'Peak_101'

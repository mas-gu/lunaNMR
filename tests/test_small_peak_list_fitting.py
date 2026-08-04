# ABOUTME: Regression tests for fitting peak lists too small to be worth parallelising.
# ABOUTME: A one- or two-peak list must come back fitted and clustered like any other.

import numpy as np
import pandas as pd

from lunaNMR.core.core_integrator import EnhancedVoigtIntegrator
from lunaNMR.processors.single_spectrum_processor import SingleSpectrumProcessor
from lunaNMR.utils.parameter_manager import NMRParameterManager

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


def _fit(count):
    """Fit a peak list the way the GUI does with parallel processing switched on."""
    processor = SingleSpectrumProcessor(_synthetic_integrator(), NMRParameterManager())
    results, _ = processor.process_peak_list(_peak_list(count), {'use_parallel': True})
    return processor, results


def test_two_peaks_are_fitted_although_parallel_would_not_pay():
    # A series run of two peaks was reporting zero fits for every spectrum.
    processor, results = _fit(2)

    assert len(results) == 2
    assert all(result['success'] for result in results)
    assert [result['assignment'] for result in results] == ['Peak_101', 'Peak_106']


def test_a_single_peak_is_fitted():
    _, results = _fit(1)

    assert len(results) == 1
    assert results[0]['success']
    assert results[0]['assignment'] == 'Peak_101'


def test_a_small_list_publishes_its_clusters():
    """Series propagation locks spectrum 1's clusters onto every later spectrum. A run
    that fits without publishing them leaves the series to recompute its own."""
    processor, _ = _fit(2)

    assert processor._computed_clusters_by_assignment
    assert sorted(sum(processor._computed_clusters_by_assignment, [])) == ['Peak_101', 'Peak_106']

# ABOUTME: A spectrum has one noise level: the load, both detection modes, the fits' S/N, the
# ABOUTME: parallel workers, the adaptive search and diagnose all read the same number.

"""
The noise level is the standard deviation of the signal-free part of the spectrum where
the peaks are. It is estimated as MAD x 1.4826 over the automatic detection region (1H
at or above 5.8 ppm for 15N, at or below 2.75 ppm for 13C, 0.5 ppm off the indirect
edges), which a few percent of signal barely moves and which holds neither the water
line nor the spectrum's edges, where the noise is not the noise under the peaks.

The reference these tests measure it against is computed here, from the spectrum and its
peak list alone: the median standard deviation of 16 x 32-point blocks inside the listed
peaks' bounding box (widened by 0.5 ppm in 1H and 4 ppm in 15N), each block at least
0.15 ppm (1H) or 1.5 ppm (15N) from every listed peak, outside 4.35-5.05 ppm (water) and
outside the outer 5% of each axis.
"""

import contextlib
import io
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

DATA = Path(__file__).resolve().parents[2] / "data_test"
BINDING = DATA / "data_binding" / "DNAJA1_HSPA8"
HETNOE = next((d for d in sorted(DATA.glob("*/*/hetNOE_800ms"))), None)
T1 = HETNOE.parent / "T1" if HETNOE is not None else None

pytestmark = pytest.mark.skipif(not BINDING.exists() or HETNOE is None,
                                reason="data_test not present")


def _spectra():
    return {
        'titration reference': (BINDING / "2D_DNAJA1_0o0.ft", BINDING / "DNAJA1.txt"),
        'T1 first delay': (next(T1.glob("*_0ms.ft")), next(T1.glob("*assignment.txt"))),
        'hetNOE saturated': (next(HETNOE.glob("*_saturated.ft")),
                             next(HETNOE.glob("*assignment.txt"))),
    }


SPECTRA = ('titration reference', 'T1 first delay', 'hetNOE saturated')


def _peaks(path):
    from lunaNMR.utils.file_manager import NMRFileManager
    with contextlib.redirect_stdout(io.StringIO()):
        return NMRFileManager().load_peak_list(str(path))


def _loaded(path):
    from lunaNMR.core.core_integrator import EnhancedVoigtIntegrator
    it = EnhancedVoigtIntegrator()
    with contextlib.redirect_stdout(io.StringIO()):
        assert it.load_nmr_file(str(path))
    return it


def _signal_free_noise(data, x_ppm, y_ppm, peaks):
    """The reference of the module docstring."""
    px = peaks['Position_X'].to_numpy(float)
    py = peaks['Position_Y'].to_numpy(float)
    ny, nx = data.shape
    edge_y, edge_x = int(np.ceil(0.05 * ny)), int(np.ceil(0.05 * nx))
    x_in = ((x_ppm >= px.min() - 0.5) & (x_ppm <= px.max() + 0.5)
            & ~((x_ppm >= 4.35) & (x_ppm <= 5.05)))
    y_in = (y_ppm >= py.min() - 4.0) & (y_ppm <= py.max() + 4.0)
    x_in[:edge_x] = x_in[nx - edge_x:] = False
    y_in[:edge_y] = y_in[ny - edge_y:] = False
    rows, columns = np.where(y_in)[0], np.where(x_in)[0]
    spreads = []
    for y0 in range(rows.min(), rows.max() - 14, 16):
        for x0 in range(columns.min(), columns.max() - 30, 32):
            if not (y_in[y0:y0 + 16].all() and x_in[x0:x0 + 32].all()):
                continue
            xs, ys = x_ppm[x0:x0 + 32], y_ppm[y0:y0 + 16]
            dx = np.maximum(0, np.maximum(xs.min() - px, px - xs.max()))
            dy = np.maximum(0, np.maximum(ys.min() - py, py - ys.max()))
            if np.any((dx < 0.15) & (dy < 1.5)):
                continue
            spreads.append(data[y0:y0 + 16, x0:x0 + 32].astype(float).std())
    assert len(spreads) > 200
    return float(np.median(spreads))


@pytest.mark.parametrize("name", SPECTRA)
def test_the_noise_at_load_is_that_of_the_signal_free_regions(name):
    spectrum, peak_list = _spectra()[name]
    it = _loaded(spectrum)
    reference = _signal_free_noise(it.nmr_data, it.ppm_x_axis, it.ppm_y_axis, _peaks(peak_list))
    assert it.noise_level == pytest.approx(reference, rel=0.10)


@pytest.mark.parametrize("name", SPECTRA)
def test_a_detected_maximum_reports_its_height_over_that_noise(name):
    it = _loaded(_spectra()[name][0])
    found = it._detect_peaks_by_threshold(10 * it.noise_level)
    assert found
    for peak in found:
        assert peak['snr'] == pytest.approx(abs(peak['intensity']) / it.noise_level, rel=1e-6)
        assert peak['snr'] >= 10


def test_detection_without_a_peak_list_thresholds_on_that_noise():
    it = _loaded(_spectra()['titration reference'][0])
    it.sn_threshold = 10.0
    it.expected_peak_count = 10 ** 6
    with contextlib.redirect_stdout(io.StringIO()):
        found = it.detect_peaks_sn_native()
    weakest = min(abs(p['intensity']) for p in found)
    assert weakest > 10 * it.noise_level
    assert len(found) == len(it._detect_peaks_by_threshold(10 * it.noise_level))


def test_a_manual_detection_region_does_not_change_the_noise():
    it = _loaded(_spectra()['titration reference'][0])
    at_load = it.noise_level
    it.manual_detection_region = [7.5, 9.0, 115.0, 125.0]
    assert it._estimate_noise_level_advanced() == at_load
    it._estimate_noise_level()
    assert it.noise_level == at_load


@pytest.mark.parametrize("name", SPECTRA)
def test_diagnose_reads_the_noise_the_integrator_reads(name):
    from lunaNMR.validation.spectra_check import _capture, peak_heights
    spectrum, peak_list = _spectra()[name]
    listed = _peaks(peak_list)
    peaks = list(zip(listed['Assignment'].astype(str), listed['Position_X'].astype(float),
                     listed['Position_Y'].astype(float)))
    noise = _loaded(spectrum).noise_level
    assert _capture(spectrum, peaks)['noise'] == pytest.approx(noise, rel=1e-6)
    assert peak_heights(spectrum, peaks)[1] == pytest.approx(noise, rel=1e-6)


def _fitted(it, peaks, **gui_params):
    from numba.core.errors import NumbaWarning
    it.gui_params = {'use_ps2d_multi_peak': True, **gui_params}
    with contextlib.redirect_stdout(io.StringIO()), warnings.catch_warnings():
        warnings.simplefilter("ignore", NumbaWarning)
        it.enhanced_fitter.enhanced_peak_fitting_parallel(peaks, parent_integrator=it,
                                                          in_process=True)
    return it.enhanced_fitter.parallel_processor


def test_a_parallel_worker_receives_the_integrators_noise():
    from lunaNMR.core.parallel_voigt_processor import _initialize_worker_fitter
    spectrum, peak_list = _spectra()['titration reference']
    it = _loaded(spectrum)
    processor = _fitted(it, _peaks(peak_list).head(3))
    try:
        worker = _initialize_worker_fitter(processor._prepare_shared_context([]))
    finally:
        processor._cleanup_shared_memory()
    assert worker.noise_level == it.noise_level
    assert worker.enhanced_fitter.noise_level == it.noise_level


def test_the_adaptive_search_scores_against_the_integrators_noise(monkeypatch):
    from lunaNMR.core import parallel_voigt_processor
    from lunaNMR.core.ps2d_config import restore_base_ps2d_config
    given = []
    optimizer = parallel_voigt_processor.AdaptiveOptimizer

    def recording(noise_level, config=None):
        given.append(noise_level)
        return optimizer(noise_level=noise_level, config=config)
    monkeypatch.setattr(parallel_voigt_processor, 'AdaptiveOptimizer', recording)

    spectrum, peak_list = _spectra()['titration reference']
    it = _loaded(spectrum)
    try:
        _fitted(it, _peaks(peak_list), use_adaptive_optimization=True)
    finally:
        restore_base_ps2d_config()
    assert given == [it.noise_level]


def _series_sigmas(out, spectra, peaks, fail_the_loader):
    """detection_sigma of every fitted point of a two-spectrum series, optionally with the
    integrator's own loader failing so each spectrum is read by the series' direct loader."""
    from numba.core.errors import NumbaWarning
    from lunaNMR.cli import _default_series_params
    from lunaNMR.processors.multi_spectrum_processor import MultiSpectrumProcessor
    processor = MultiSpectrumProcessor(_default_series_params())
    if fail_the_loader:
        processor.integrator.load_nmr_file = lambda *args, **kwargs: False
    with contextlib.redirect_stdout(io.StringIO()), warnings.catch_warnings():
        warnings.simplefilter("ignore", NumbaWarning)
        processor.process_nmr_series([str(s) for s in spectra], peaks, str(out),
                                     peak_source_mode="reference", extract_delays=True)
    tidy = pd.read_csv(out / "series_analysis_tidy.csv", dtype={'assignment': str})
    return tidy.set_index(['spectrum_name', 'assignment'])['detection_sigma']


def test_a_spectrum_read_by_the_series_direct_loader_gets_its_own_noise(tmp_path):
    """The two titration points differ in noise by 1.4x: a point that kept the noise of
    the one before it reports every S/N off by that factor."""
    spectra = [BINDING / "2D_DNAJA1_HSPA8_0o5.ft", BINDING / "2D_DNAJA1_HSPA8_1o0.ft"]
    peaks = _peaks(BINDING / "DNAJA1.txt").head(4)
    loaded = _series_sigmas(tmp_path / "loaded", spectra, peaks, fail_the_loader=False)
    direct = _series_sigmas(tmp_path / "direct", spectra, peaks, fail_the_loader=True)
    assert len(direct) == 8
    assert direct.to_dict() == pytest.approx(loaded.to_dict(), rel=1e-6)

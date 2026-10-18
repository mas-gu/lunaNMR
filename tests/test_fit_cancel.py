# ABOUTME: A stop event ends a spectrum fit between clusters, in process and over the pool.
# ABOUTME: Real DNAJA1 0-eq spectrum and its peak list; an unset event leaves every result unchanged.

"""
Cancel used to close the progress dialog and nothing else: the fit ran to its end. The
fit now takes a stop event, checked before each cluster in process and while waiting on
each cluster over the pool, where setting it also ends the pool. A fit that is stopped
raises Cancelled and returns nothing, rather than falling back to a per-peak fit or
handing back the clusters it reached.
"""

import contextlib
import io
import threading
import time
from pathlib import Path

import numpy as np
import pytest

DNAJA1 = Path(__file__).resolve().parents[2] / "data_test" / "data_binding" / "DNAJA1_HSPA8"
SPECTRUM = "2D_DNAJA1_0o0.ft"

pytestmark = pytest.mark.skipif(not DNAJA1.exists(), reason="data_test not present")

MODES = pytest.mark.parametrize("in_process", [True, False], ids=["in-process", "pool"])


@pytest.fixture(scope="module")
def integrator():
    from lunaNMR.core.core_integrator import EnhancedVoigtIntegrator
    from lunaNMR.core.ps2d_config import set_ps2d_config
    from lunaNMR.utils.file_manager import NMRFileManager
    set_ps2d_config('15N')
    it = EnhancedVoigtIntegrator()
    with contextlib.redirect_stdout(io.StringIO()):
        it.load_nmr_file(str(DNAJA1 / SPECTRUM))
    it.noise_level = it._estimate_noise_level_advanced()
    it.peak_list = NMRFileManager().load_peak_list(str(DNAJA1 / "DNAJA1.txt"))
    it.gui_params = {'use_ps2d_multi_peak': True, 'use_adaptive_optimization': False}
    return it


def _fit(integrator, in_process, **kwargs):
    integrator.enhanced_fitter.series_params = None
    results, _ = integrator.enhanced_fitter.enhanced_peak_fitting_parallel(
        integrator.peak_list, parent_integrator=integrator, in_process=in_process, **kwargs)
    return results


@MODES
def test_an_unset_stop_changes_no_result(integrator, in_process):
    plain = _fit(integrator, in_process)
    watched = _fit(integrator, in_process, stop=threading.Event())

    assert len(plain) == len(watched) == len(integrator.peak_list)
    for key in ("center_x", "center_y", "height", "volume", "r_squared"):
        np.testing.assert_array_equal([r.get(key) for r in plain], [r.get(key) for r in watched])


@MODES
def test_a_stop_set_before_the_fit_fits_nothing(integrator, in_process):
    from lunaNMR.utils.cancellation import Cancelled
    stop = threading.Event()
    stop.set()
    fitted = []

    with pytest.raises(Cancelled):
        _fit(integrator, in_process, stop=stop,
             progress_callback=lambda *args: fitted.append(args))

    assert fitted == []


@MODES
def test_a_stop_set_during_the_fit_ends_it_at_once(integrator, in_process):
    """The stop is set when the first cluster reports; nothing else is fitted after it."""
    from lunaNMR.utils.cancellation import Cancelled
    stop = threading.Event()
    reported = []
    set_at = []

    def progress(*args):
        reported.append(args)
        if not stop.is_set():
            stop.set()
            set_at.append(time.monotonic())

    with pytest.raises(Cancelled):
        _fit(integrator, in_process, stop=stop, progress_callback=progress)
    latency = time.monotonic() - set_at[0]

    # One cluster's peaks are reported together, and none after it
    first_cluster = {args[1] for args in reported}
    assert len(first_cluster) == 1
    assert latency < 1.0


def test_a_stopped_spectrum_fit_raises_instead_of_reporting_a_failure(integrator):
    from lunaNMR.processors.single_spectrum_processor import SingleSpectrumProcessor
    from lunaNMR.utils.cancellation import Cancelled
    from lunaNMR.utils.parameter_manager import NMRParameterManager
    stop = threading.Event()
    stop.set()
    processor = SingleSpectrumProcessor(integrator, NMRParameterManager())

    with pytest.raises(Cancelled):
        processor.process_peak_list(integrator.peak_list,
                                    {'use_parallel': False, 'use_voigt_fitting': True,
                                     'skip_series_params': True}, stop=stop)


def test_a_stop_while_the_pool_starts_lands_at_once(integrator):
    """Starting the pool waits on a spawned process importing lunaNMR, seconds on a busy
    machine; a stop set as the first pass is handed out does not wait for it."""
    from lunaNMR.utils.cancellation import Cancelled
    from lunaNMR.utils.output_manager import OutputManager
    stop = threading.Event()
    set_at = []

    def log(progress, task, log_msg=None, failed=False):
        if not stop.is_set() and "Overlap groups" in (log_msg or task or ""):
            stop.set()
            set_at.append(time.monotonic())

    OutputManager.set_callback(log)
    try:
        with pytest.raises(Cancelled):
            _fit(integrator, in_process=False, stop=stop)
    finally:
        OutputManager.set_callback(None)

    assert time.monotonic() - set_at[0] < 0.5

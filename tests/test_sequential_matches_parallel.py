# ABOUTME: Sequential and parallel fitting are one pipeline, run in-process or over a pool.
# ABOUTME: Real DNAJA1 titration: both modes agree exactly, and the workers keep every series bound.

"""
Sequential and parallel were two algorithms. Parallel learned linewidths from the
isolated peaks, ran the adaptive radius search and re-clustered; sequential did none of
it. And the parallel workers, which rebuild the integrator from a serialised copy, lost
three things a series sets on it: the titration territories, the per-peak linewidth caps
and the cascade drift limit. On the DNAJA1 titration the two modes then tracked the same
peak to positions tens of ppb apart.

Sequential is now the parallel pipeline executed in this process, so the two cannot
disagree, and the workers receive what the series set.
"""

import contextlib
import io
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

DNAJA1 = Path(__file__).resolve().parents[2] / "data_test" / "data_binding" / "DNAJA1_HSPA8"
POINTS = ("2D_DNAJA1_0o0.ft", "2D_DNAJA1_HSPA8_0o125.ft", "2D_DNAJA1_HSPA8_0o25.ft")

pytestmark = pytest.mark.skipif(not DNAJA1.exists(), reason="data_test not present")


@pytest.fixture(scope="module")
def integrator():
    from lunaNMR.core.core_integrator import EnhancedVoigtIntegrator
    from lunaNMR.utils.file_manager import NMRFileManager
    it = EnhancedVoigtIntegrator()
    with contextlib.redirect_stdout(io.StringIO()):
        it.load_nmr_file(str(DNAJA1 / POINTS[0]))
    it.noise_level = it._estimate_noise_level_advanced()
    it.peak_list = NMRFileManager().load_peak_list(str(DNAJA1 / "DNAJA1.txt"))
    it.gui_params = {'use_ps2d_multi_peak': True, 'use_adaptive_optimization': False}
    return it


def _worker_for(integrator, **context):
    from lunaNMR.core.parallel_voigt_processor import (
        ParallelVoigtProcessor, _initialize_worker_fitter)
    fitter = integrator.enhanced_fitter
    fitter.nmr_data, fitter.ppm_x_axis, fitter.ppm_y_axis = (
        integrator.nmr_data, integrator.ppm_x_axis, integrator.ppm_y_axis)
    processor = ParallelVoigtProcessor(fitter, enable_adaptive=False)
    try:
        shared = processor._prepare_shared_context([], **context)
        return _initialize_worker_fitter(shared), shared
    finally:
        processor._cleanup_shared_memory()


def test_a_worker_keeps_the_titration_territories(integrator):
    integrator.titration_pos_margins = {'3LysH': (0.05, 0.005)}
    integrator.titration_lw_margins = {'3LysH': (0.3, 0.04)}
    try:
        worker, _ = _worker_for(integrator)
    finally:
        integrator.titration_pos_margins = integrator.titration_lw_margins = None
    assert worker.titration_pos_margins == {'3LysH': (0.05, 0.005)}
    assert worker.titration_lw_margins == {'3LysH': (0.3, 0.04)}


def test_a_worker_keeps_the_cascade_drift_limit(integrator):
    _, shared = _worker_for(integrator, reference_positions={'3LysH': (8.39, 127.56)})
    assert shared['reference_positions'] == {'3LysH': (8.39, 127.56)}


@pytest.mark.parametrize("in_process", [False, True], ids=["pool", "in-process"])
def test_the_drift_limit_holds_a_clustered_peak(integrator, in_process):
    """A peak whose reference sits 8 ppb off in 1H ends within max_drift (5 ppb) of it. (At
    20 ppb the drift window and the fit's own position window do not overlap.)"""
    from lunaNMR.core.ps2d_config import get_ps2d_config
    contexts = [{'assignment': str(r.Assignment), 'x_ppm': r.Position_X, 'y_ppm': r.Position_Y,
                 'pos_x': r.Position_X, 'pos_y': r.Position_Y, 'intensity': None}
                for r in integrator.peak_list.itertuples()]
    member = next(c for cluster in integrator.identify_overlap_clusters(contexts)
                  if len(cluster) > 1 for c in contexts if (c['x_ppm'], c['y_ppm']) == cluster[0])
    offset = member["x_ppm"] + 0.008
    results, _ = integrator.enhanced_fitter.enhanced_peak_fitting_parallel(
        integrator.peak_list, parent_integrator=integrator, in_process=in_process,
        reference_positions={member['assignment']: (offset, member['y_ppm'])})
    fitted = next(r for r in results if str(r.get('assignment')) == member['assignment'])
    assert abs(fitted['center_x'] - offset) <= get_ps2d_config().max_drift_f2 + 1e-9


@pytest.fixture(scope="module")
def series(tmp_path_factory):
    """The first three DNAJA1 titration points, fitted once in each mode."""
    from lunaNMR.cli import main
    root = tmp_path_factory.mktemp("modes")
    (root / "spectra").mkdir()
    for name in POINTS:
        (root / "spectra" / name).symlink_to(DNAJA1 / name)
    out = {}
    for mode, flags in (("sequential", []), ("parallel", ["--parallel"])):
        assert main(["series", "--spectra", str(root / "spectra"),
                     "--peaks", str(DNAJA1 / "DNAJA1.txt"), "--mode", "titration",
                     "--peak-source", "cascade", "--out", str(root / mode), *flags]) == 0
        out[mode] = pd.read_csv(root / mode / "series_analysis_tidy.csv")
    return out


@pytest.mark.parametrize("column", ["ppm_x", "ppm_y", "height", "volume", "r_squared"])
def test_sequential_and_parallel_give_the_same_series(series, column):
    key = ["spectrum_name", "assignment"]
    seq = series["sequential"].set_index(key).sort_index()
    par = series["parallel"].set_index(key).sort_index()
    assert seq.index.equals(par.index)
    np.testing.assert_array_equal(seq[column].to_numpy(), par[column].to_numpy())

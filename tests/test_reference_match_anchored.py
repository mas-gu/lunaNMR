# ABOUTME: Later spectra of a series match each reference only near its listed position.
# ABOUTME: Real 15N T2 planes: a peak that decays below threshold must not take an unlisted neighbour.

"""
The indirect search window is wide enough (0.3 ppm in 15N) to reach a peak one 15N point
off its reference. In a relaxation series a peak decays: once its own maximum falls
below the detection threshold, the wide box reaches a stronger unlisted maximum beside
it. Reference 44LysH (7.713, 120.911) has one 0.063 ppm 1H / 0.14 ppm 15N away, five
times its height at 102 ms, and fitting there reported it Excellent with a T2 eight
times its own.

On every spectrum after the first, a reference may only take a maximum within the
0.070 ppm box around its listed position. The first spectrum, and a single spectrum, keep
the wide window. The box is centred on the list, not on the first spectrum's match: a
hetNOE pair is a two-spectrum series, and anchoring the strong unsaturated plane on the
weak saturated one sent a peak's unsaturated intensity off its own maximum.
"""

import contextlib
import io
import warnings
from pathlib import Path

import pandas as pd
import pytest

DATA = Path(__file__).resolve().parents[2] / "data_test"
T2 = next((d for d in sorted(DATA.glob("*/600_*/T2"))
           if (d / "600_T2_assignment.txt").exists()
           and "44LysH" in (d / "600_T2_assignment.txt").read_text()), None)

pytestmark = pytest.mark.skipif(T2 is None, reason="data_test not present")

LISTED = (7.712760, 120.910511)


def _plane(delay):
    return next(p for p in sorted(T2.glob(f"600_T2_*_{delay}.ft")) if "_b_" not in p.name)


def _peak_list():
    from lunaNMR.utils.file_manager import NMRFileManager
    with contextlib.redirect_stdout(io.StringIO()):
        return NMRFileManager().load_peak_list(str(T2 / "600_T2_assignment.txt"))


@pytest.fixture(autouse=True)
def nucleus_15n():
    from lunaNMR.core.ps2d_config import set_ps2d_config
    set_ps2d_config('15N')
    yield
    set_ps2d_config('15N')


def _detect(plane, anchors=None):
    from lunaNMR.core.core_integrator import EnhancedVoigtIntegrator
    it = EnhancedVoigtIntegrator()
    with contextlib.redirect_stdout(io.StringIO()):
        it.load_nmr_file(str(plane))
        it.peak_list = _peak_list()
        it.peak_list_contour_threshold = None
        it.peak_list_sn_threshold = 10.0
        it.auto_add_dummy_peaks = False
        return {str(p['assignment']): p for p in it._detect_peaks_reference_based(anchors=anchors)}


def test_a_faded_peak_does_not_take_a_neighbour_beyond_its_listed_box():
    listed = {str(r['Assignment']): (float(r['Position_X']), float(r['Position_Y']))
              for _, r in _peak_list().iterrows()}

    late = _detect(_plane("102ms"), anchors=listed)

    assert abs(late['44LysH']['ppm_x'] - LISTED[0]) <= 0.070
    assert abs(late['44LysH']['ppm_y'] - LISTED[1]) <= 0.070


def test_a_single_spectrum_keeps_the_wide_window():
    late = _detect(_plane("102ms"))
    assert late['44LysH']['ppm_y'] - LISTED[1] > 0.070


def test_a_series_keeps_a_faded_peak_on_its_own_position(tmp_path):
    from numba.core.errors import NumbaWarning
    from lunaNMR.cli import _default_series_params
    from lunaNMR.processors.multi_spectrum_processor import MultiSpectrumProcessor
    with contextlib.redirect_stdout(io.StringIO()), warnings.catch_warnings():
        warnings.simplefilter("ignore", NumbaWarning)
        MultiSpectrumProcessor(_default_series_params()).process_nmr_series(
            [str(_plane("8ms")), str(_plane("102ms"))], _peak_list(), str(tmp_path / "out"),
            peak_source_mode="reference", extract_delays=True)
    tidy = pd.read_csv(tmp_path / "out" / "series_analysis_tidy.csv", dtype={'assignment': str})
    rows = tidy[tidy.assignment == '44LysH'].set_index('spectrum_name')
    assert len(rows) == 2
    for _, row in rows.iterrows():
        assert abs(row['ppm_x'] - LISTED[0]) < 0.03
        assert abs(row['ppm_y'] - LISTED[1]) < 0.1

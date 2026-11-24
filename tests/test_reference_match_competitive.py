# ABOUTME: References compete for detected maxima, the pair closest in linewidths settled first.
# ABOUTME: Real 15N T1 and hetNOE planes: maxima a reference sits on, taken by a reference earlier in the list.

"""
References were matched first-come in list order, each taking the closest maximum still
free in its search box. On the first T1 plane at 10 sigma, reference 64.0 has no maximum
of its own and, listed before 78.0, takes the one 78.0 sits exactly on, leaving 78.0 read
at its listed position. A wider indirect window makes this more common: on the saturated
hetNOE plane 70.0 reaches the maximum 98.0 sits on, 0.076 ppm away in 15N, only through
the 0.3 ppm window.

All reference-maximum pairs are now settled together, the pair closest in linewidths
first, so a maximum goes to the reference nearest it and the other is read at its
reference position.
"""

import contextlib
import io
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

DATA = Path(__file__).resolve().parents[2] / "data_test"
HETNOE = next((d for d in sorted(DATA.glob("*/*/hetNOE_800ms"))), None)
T1 = HETNOE.parent / "T1" if HETNOE is not None else None

pytestmark = pytest.mark.skipif(HETNOE is None, reason="data_test not present")


def _peak_list(folder):
    peaks = pd.read_csv(next(folder.glob("*assignment.txt")), skipinitialspace=True)
    peaks.columns = peaks.columns.str.strip()
    peaks['Assignment'] = peaks['Assignment'].astype(str)
    return peaks[['Assignment', 'Position_X', 'Position_Y']]


def _match(plane, contour_threshold=None, window_y=None, folder=HETNOE):
    """Reference detection on one plane with the default windows, at the GUI contour
    threshold or at the CLI's 10 sigma."""
    from lunaNMR.core.core_integrator import EnhancedVoigtIntegrator
    it = EnhancedVoigtIntegrator()
    with contextlib.redirect_stdout(io.StringIO()):
        it.load_nmr_file(str(next(folder.glob(plane))))
    it.peak_list_contour_threshold = (None if contour_threshold is None
                                      else float(np.max(it.nmr_data)) * contour_threshold)
    it.peak_list_sn_threshold = 10.0
    it.auto_add_dummy_peaks = False
    if window_y is not None:
        it.search_window_y = window_y
    it.peak_list = _peak_list(folder)
    with contextlib.redirect_stdout(io.StringIO()):
        return {p['assignment']: p for p in it._detect_peaks_reference_based()}


@pytest.fixture(autouse=True)
def nucleus_15n():
    from lunaNMR.core.ps2d_config import set_ps2d_config
    set_ps2d_config('15N')
    yield
    set_ps2d_config('15N')


def _keeps_its_maximum(matched, owner, x, y, loser, loser_ref):
    assert matched[owner]['detected']
    assert matched[owner]['ppm_x'] == pytest.approx(x, abs=1e-4)
    assert matched[owner]['ppm_y'] == pytest.approx(y, abs=1e-3)
    # the other has no maximum of its own above 10 sigma and is read where it is listed
    assert not matched[loser]['detected']
    assert (matched[loser]['ppm_x'], matched[loser]['ppm_y']) == loser_ref


def test_a_maximum_goes_to_the_reference_nearest_it():
    matched = _match("*_0ms.ft", folder=T1)
    _keeps_its_maximum(matched, '78.0', 8.113845, 121.845543, '64.0', (8.156951, 121.788591))


def test_a_wider_window_does_not_let_a_reference_take_anothers_maximum():
    matched = _match("*_saturated.ft", window_y=0.3)
    _keeps_its_maximum(matched, '98.0', 7.968417, 117.530987, '70.0', (7.936056, 117.607417))


def test_a_peak_one_15n_point_off_its_reference_is_reached():
    """At 600 MHz the planes have 0.074 ppm per 15N point, so a 0.070 ppm indirect window
    cannot reach 20.0's maximum one point off its reference (7.600, 117.111)."""
    matched = _match("*_saturated.ft", contour_threshold=0.05 * 1.05)

    assert matched['20.0']['detected']
    assert matched['20.0']['ppm_x'] == pytest.approx(7.597833, abs=1e-4)
    assert matched['20.0']['ppm_y'] == pytest.approx(116.963144, abs=1e-3)


def test_no_maximum_is_taken_by_two_references():
    for plane, contour in (("*_saturated.ft", 0.05 * 1.05), ("*_unsaturated.ft", None)):
        matched = _match(plane, contour_threshold=contour).values()
        positions = [(p['x_point'], p['y_point']) for p in matched if p['detected']]
        assert len(positions) == len(set(positions))

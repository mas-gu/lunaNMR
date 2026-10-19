# ABOUTME: References compete for detected maxima, the pair closest in linewidths settled first.
# ABOUTME: Real 15N hetNOE planes: maxima a reference sits on, taken by a reference earlier in the list.

"""
References were matched first-come in list order, each taking the closest maximum still
free in its search box. On the unsaturated plane at 10 sigma, reference 24.0 has no
maximum of its own and took the one reference 89.0 sits exactly on, leaving 89.0 read at
its listed position. A wider indirect window makes this more common: at 0.3 ppm, 5.0
takes the maximum 53.0 sits on.

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

pytestmark = pytest.mark.skipif(HETNOE is None, reason="data_test not present")


def _peak_list():
    peaks = pd.read_csv(next(HETNOE.glob("*assignment.txt")), skipinitialspace=True)
    peaks.columns = peaks.columns.str.strip()
    peaks['Assignment'] = peaks['Assignment'].astype(str)
    return peaks[['Assignment', 'Position_X', 'Position_Y']]


def _match(plane, contour_threshold=None, window_y=None):
    """Reference detection on one plane with the default windows, at the GUI contour
    threshold or at the CLI's 10 sigma."""
    from lunaNMR.core.core_integrator import EnhancedVoigtIntegrator
    it = EnhancedVoigtIntegrator()
    with contextlib.redirect_stdout(io.StringIO()):
        it.load_nmr_file(str(next(HETNOE.glob(plane))))
    it.peak_list_contour_threshold = (None if contour_threshold is None
                                      else float(np.max(it.nmr_data)) * contour_threshold)
    it.peak_list_sn_threshold = 10.0
    it.auto_add_dummy_peaks = False
    if window_y is not None:
        it.search_window_y = window_y
    it.peak_list = _peak_list()
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
    matched = _match("*_unsaturated.ft")
    _keeps_its_maximum(matched, '89.0', 8.067644, 114.216473, '24.0', (8.006151, 114.179237))


def test_a_wider_window_does_not_let_a_reference_take_anothers_maximum():
    matched = _match("*_unsaturated.ft", window_y=0.3)
    _keeps_its_maximum(matched, '53.0', 9.129978, 124.395505, '5.0', (9.097587, 124.508188))


@pytest.mark.xfail(reason="needs the wider indirect window (option 2)", strict=True)
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

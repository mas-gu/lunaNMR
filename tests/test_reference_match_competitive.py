# ABOUTME: References compete for detected maxima, nearest pair first, inside a window that reaches one 15N point.
# ABOUTME: Real 15N hetNOE planes: a peak one 15N point off its reference, and a maximum two references can reach.

"""
At 600 MHz the hetNOE planes have 0.074 ppm per 15N point, so a 0.070 ppm indirect window
cannot reach a peak sitting one point away from its reference: on the saturated plane,
reference 20.0 (7.600, 117.111) has its maximum at (7.598, 116.963) and was matched to
a weaker one beside it. The indirect window of the nucleus (0.3 ppm for 15N) reaches it.

A wider window reaches more maxima, and references were matched first-come in list
order: on the unsaturated plane reference 5.0 has no maximum of its own above 10 sigma,
and a 0.3 ppm window let it take the one reference 53.0 sits exactly on. All
reference-maximum pairs are now settled together, the pair closest in linewidths first,
so 53.0 keeps its maximum and 5.0 is read at its reference position.
"""

import contextlib
import io
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

DATA = Path(__file__).resolve().parents[2] / "data_test"
HETNOE = next((d for d in sorted(DATA.glob("*/*/hetNOE_800ms"))), None)


def _peak_list():
    peaks = pd.read_csv(next(HETNOE.glob("*assignment.txt")), skipinitialspace=True)
    peaks.columns = peaks.columns.str.strip()
    peaks['Assignment'] = peaks['Assignment'].astype(str)
    return peaks[['Assignment', 'Position_X', 'Position_Y']]


def _match(plane, contour_threshold=None, window_y=None):
    """Reference detection on one plane with the default windows of a 15N spectrum,
    at the GUI contour threshold or at the CLI's 10 sigma."""
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


@pytest.mark.skipif(HETNOE is None, reason="data_test not present")
def test_a_peak_one_15n_point_off_its_reference_is_reached():
    matched = _match("*_saturated.ft", contour_threshold=0.05 * 1.05)

    assert matched['20.0']['detected']
    assert matched['20.0']['ppm_x'] == pytest.approx(7.597833, abs=1e-4)
    assert matched['20.0']['ppm_y'] == pytest.approx(116.963144, abs=1e-3)


@pytest.mark.skipif(HETNOE is None, reason="data_test not present")
@pytest.mark.parametrize("window_y", [None, 0.3])
def test_a_maximum_goes_to_the_reference_nearest_it(window_y):
    matched = _match("*_unsaturated.ft", window_y=window_y)

    assert matched['53.0']['detected']
    assert matched['53.0']['ppm_x'] == pytest.approx(9.129978, abs=1e-4)
    assert matched['53.0']['ppm_y'] == pytest.approx(124.395505, abs=1e-3)
    # 5.0 has no maximum of its own above 10 sigma and is read where it is listed
    assert not matched['5.0']['detected']
    assert (matched['5.0']['ppm_x'], matched['5.0']['ppm_y']) == (9.097587, 124.508188)


@pytest.mark.skipif(HETNOE is None, reason="data_test not present")
def test_no_maximum_is_taken_by_two_references():
    for plane, contour in (("*_saturated.ft", 0.05 * 1.05), ("*_unsaturated.ft", None)):
        matched = _match(plane, contour_threshold=contour).values()
        positions = [(p['x_point'], p['y_point']) for p in matched if p['detected']]
        assert len(positions) == len(set(positions))

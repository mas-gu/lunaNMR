# ABOUTME: A reference takes the detected maximum closest to it in linewidths, not in raw ppm.
# ABOUTME: Real 15N hetNOE plane where raw ppm swaps two peaks, and a real 13C methyl spectrum.

"""
Within the search box, the closest maximum was the one at the smallest Euclidean
distance in ppm, which counts a 15N ppm as a 1H ppm. A 15N linewidth is ten 1H
linewidths, so a maximum 0.036 ppm away in 1H beat one 0.027 ppm away in 1H and 0.041 in
15N: half a 15N point, a quarter of a linewidth. On the unsaturated plane below, reference
101.0 took the maximum that sits exactly on reference 172.0, and 172.0 was left the one
that belongs to 101.0.

The distance is scaled by the PS2D radii of the spectrum's nucleus (radF2 in 1H, radF1
in the indirect dimension), so a 13C spectrum is measured in 13C linewidths.
"""

import contextlib
import io
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

DATA = Path(__file__).resolve().parents[2] / "data_test"
HETNOE = next(iter(sorted(DATA.glob("*/*/hetNOE_800ms"))), None)
METHYL = DATA / "methyl" / "sampleNBD"


def _integrator(spectrum, nucleus):
    from lunaNMR.core.core_integrator import EnhancedVoigtIntegrator
    from lunaNMR.core.ps2d_config import set_ps2d_config
    set_ps2d_config(nucleus)
    it = EnhancedVoigtIntegrator()
    with contextlib.redirect_stdout(io.StringIO()):
        it.load_nmr_file(str(spectrum))
    it.search_window_x = it.search_window_y = 0.070
    it.peak_list_contour_threshold = None
    it.peak_list_sn_threshold = 10.0
    it.auto_add_dummy_peaks = False
    return it


@pytest.fixture
def restore_nucleus():
    from lunaNMR.core.ps2d_config import set_ps2d_config
    yield
    set_ps2d_config('15N')


def _match(it, rows):
    it.peak_list = pd.DataFrame(rows, columns=['Assignment', 'Position_X', 'Position_Y'])
    with contextlib.redirect_stdout(io.StringIO()):
        return {p['assignment']: p for p in it._detect_peaks_reference_based()}


@pytest.mark.skipif(HETNOE is None, reason="data_test not present")
def test_each_reference_takes_the_maximum_closest_in_linewidths(restore_nucleus):
    it = _integrator(next(HETNOE.glob("*_unsaturated.ft")), '15N')
    matched = _match(it, [('101.0', 7.816826, 117.550025), ('172.0', 7.853040, 117.552897)])

    assert matched['172.0']['ppm_x'] == pytest.approx(7.853040, abs=1e-4)
    assert matched['172.0']['ppm_y'] == pytest.approx(117.552897, abs=1e-3)
    assert matched['101.0']['ppm_x'] == pytest.approx(7.789647, abs=1e-4)
    assert matched['101.0']['ppm_y'] == pytest.approx(117.590594, abs=1e-3)


@pytest.mark.skipif(not METHYL.exists(), reason="data_test not present")
def test_a_13c_spectrum_is_measured_in_13c_linewidths(restore_nucleus):
    from lunaNMR.core.ps2d_config import get_ps2d_config
    it = _integrator(METHYL / "001_T2_NBD_3ms.ft", '13C')
    peaks = pd.read_csv(METHYL / "001_NBDADP_NBD.txt", skipinitialspace=True)
    peaks.columns = peaks.columns.str.strip()
    matched = [p for p in _match(it, peaks[['Assignment', 'Position_X', 'Position_Y']].values).values()
               if p['detected'] and p['distance_from_reference_y'] > 0]
    assert matched, "the 3 ms methyl spectrum has references off their maxima in 13C"

    config = get_ps2d_config()
    assert config.radF1 != 0.4
    for p in matched:
        assert p['distance_from_reference_elliptical'] == pytest.approx(np.hypot(
            p['distance_from_reference_x'] / config.radF2,
            p['distance_from_reference_y'] / config.radF1))

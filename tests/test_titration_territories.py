# ABOUTME: A titration sizes each peak's position window against its nearest neighbour.
# ABOUTME: Without it the shared window exceeds the gap between any two clustered peaks.

import pandas as pd

from lunaNMR.processors.multi_spectrum_processor import MultiSpectrumProcessor


def _processor():
    # _titration_territories reads nothing off the instance.
    return MultiSpectrumProcessor.__new__(MultiSpectrumProcessor)


WINDOW_F1, WINDOW_F2 = 0.3, 0.03


def _peaks(rows):
    return pd.DataFrame(rows, columns=['Assignment', 'Position_X', 'Position_Y'])


def test_territories_are_keyed_by_assignment():
    """The cluster being fitted knows its members by assignment and nothing else."""
    peaks = _peaks([['3', 7.70, 121.50], ['4', 8.90, 118.20]])
    territories = _processor()._titration_territories(peaks, WINDOW_F1, WINDOW_F2)

    assert set(territories) == {'3', '4'}
    assert territories['3'] == (WINDOW_F1, WINDOW_F2), 'distant peaks keep the window'


def test_a_clustered_pair_gives_up_room():
    """0.02 ppm apart in 1H is inside the 0.04 ppm overlap threshold, so these two
    are fitted together - and the 0.03 ppm window would let either cross the gap."""
    peaks = _peaks([['27', 7.700, 121.50], ['60', 7.720, 121.52]])
    territories = _processor()._titration_territories(peaks, WINDOW_F1, WINDOW_F2)

    assert territories['27'][1] < 0.010
    assert territories['27'][1] + territories['60'][1] < 0.020, 'they still meet'


def test_axes_are_read_the_way_the_fitter_reads_them():
    """Position_X is 1H (F2) and Position_Y the indirect dimension (F1). Swapped,
    a peak list of well-separated peaks would come back with tiny 1H margins."""
    peaks = _peaks([['3', 7.70, 121.5], ['4', 7.71, 121.5]])
    territories = _processor()._titration_territories(peaks, WINDOW_F1, WINDOW_F2)

    assert territories['3'][0] < WINDOW_F1 and territories['3'][1] < WINDOW_F2
    # The pair is 0.01 apart in 1H and identical in 15N, so 1H is what constrains it.
    assert territories['3'][1] < 0.005


def test_a_peak_list_without_positions_gives_no_territories():
    """A peak list that cannot be read must leave the shared window in charge
    rather than half-constrain the fit."""
    peaks = pd.DataFrame({'Assignment': ['3', '4']})
    assert _processor()._titration_territories(peaks, WINDOW_F1, WINDOW_F2) is None
    assert _processor()._titration_territories(None, WINDOW_F1, WINDOW_F2) is None

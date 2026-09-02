# ABOUTME: Tests the per-peak position bounds that stop a tracked peak reaching a neighbour.
# ABOUTME: The property under test is disjointness: no position is reachable by two peaks.

import pytest

from lunaNMR.core.ps2d_territories import territory_margins


MARGIN_F1, MARGIN_F2 = 0.3, 0.03          # a typical 15N titration per-step window


def _boxes_are_disjoint(positions, margins):
    """True when no two peaks can reach a common position.

    Two axis-aligned boxes miss each other as soon as they are apart in ONE
    dimension, so that is what this checks per pair.
    """
    for i, ((f1_i, f2_i), (m1_i, m2_i)) in enumerate(zip(positions, margins)):
        for (f1_j, f2_j), (m1_j, m2_j) in list(zip(positions, margins))[i + 1:]:
            apart = (abs(f1_i - f1_j) > m1_i + m1_j
                     or abs(f2_i - f2_j) > m2_i + m2_j)
            if not apart:
                return False
    return True


def test_a_lone_peak_keeps_the_whole_window():
    """With nothing to collide with there is nothing to give up."""
    assert territory_margins([(120.0, 8.0)], MARGIN_F1, MARGIN_F2) == \
           [(MARGIN_F1, MARGIN_F2)]


def test_distant_peaks_keep_the_whole_window():
    """The bound only tightens where a neighbour is within reach."""
    positions = [(120.0, 8.0), (124.0, 9.0)]
    margins = territory_margins(positions, MARGIN_F1, MARGIN_F2)
    assert margins == [(MARGIN_F1, MARGIN_F2)] * 2


def test_a_close_neighbour_takes_the_window_below_half_the_gap():
    """A clustered pair is the case that matters: peaks within the overlap
    threshold (0.04 ppm in 1H) sit closer than the 0.03 ppm per-step window, so
    without this each is free to cross the midpoint in a single step."""
    positions = [(121.50, 7.690), (121.52, 7.710)]      # 0.02 ppm apart in 1H
    margins = territory_margins(positions, MARGIN_F1, MARGIN_F2)

    assert all(m2 <= 0.010 for _, m2 in margins), 'a peak may still reach the other'
    assert _boxes_are_disjoint(positions, margins)


def test_the_window_shrinks_by_the_dimension_that_separates():
    """Separation in either dimension alone is enough to keep peaks apart, and the
    two dimensions are not comparable in ppm - 0.02 ppm is a wide gap in 1H and no
    gap at all in 15N. Scaling by the requested window is what makes them so."""
    inF1 = territory_margins([(121.0, 8.0), (121.2, 8.0)], MARGIN_F1, MARGIN_F2)
    inF2 = territory_margins([(121.0, 8.0), (121.0, 8.02)], MARGIN_F1, MARGIN_F2)

    assert inF1[0][0] < MARGIN_F1 and inF2[0][1] < MARGIN_F2
    # 0.2 ppm in 15N is 2/3 of that window; 0.02 ppm in 1H is 2/3 of that one.
    assert inF1[0][0] / MARGIN_F1 == pytest.approx(inF2[0][1] / MARGIN_F2)


def test_no_position_is_reachable_by_two_peaks():
    """The whole point: a peak cannot arrive where another peak could also be, so
    it cannot take that peak's identity."""
    positions = [(121.50, 7.690), (121.52, 7.710), (121.75, 7.660),
                 (118.20, 8.900), (121.51, 7.700), (124.00, 9.100)]
    margins = territory_margins(positions, MARGIN_F1, MARGIN_F2)

    assert len(margins) == len(positions)
    assert _boxes_are_disjoint(positions, margins)


def test_peaks_driven_together_still_cannot_meet():
    """Bounds are hard clips, so two peaks pushed at each other both end up ON
    their bound. Territories that merely touch would let both land on the shared
    face - coincident, which is exactly what this exists to prevent."""
    positions = [(121.50, 7.700), (121.50, 7.720)]
    (_, m2_a), (_, m2_b) = territory_margins(positions, MARGIN_F1, MARGIN_F2)

    assert m2_a + m2_b < 0.020, 'the territories meet at the midpoint'


def test_coincident_peaks_are_given_no_room():
    """Nothing sensible to hand them - they are already the failure this prevents,
    and the caller refuses them by name. It must not divide by their zero gap."""
    margins = territory_margins([(121.5, 7.7), (121.5, 7.7)], MARGIN_F1, MARGIN_F2)
    assert margins == [(0.0, 0.0), (0.0, 0.0)]


def test_an_empty_peak_list_is_not_an_error():
    assert territory_margins([], MARGIN_F1, MARGIN_F2) == []

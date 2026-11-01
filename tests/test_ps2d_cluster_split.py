# ABOUTME: An oversized cluster is split along its links, so every part stays spatially connected.
# ABOUTME: Peaks on the real integrator; positions only, no spectrum is read.

"""
A cluster over `max_cluster_size` used to be cut by sorting its peaks on 1H and splitting at
a gap, so a part could hold peaks that touch none of each other and skip over peaks of
another cluster. It is now split by removing the longest link (elliptical distance, which
weighs both dimensions) until it falls apart.

Two tight groups 1.4 ppm apart in 15N, joined only by one peak between them: no peak of one
group touches a peak of the other, so no part may hold both.
"""

import pytest

from lunaNMR.core.core_integrator import VoigtIntegrator
from lunaNMR.core.ps2d_config import set_ps2d_config

GROUP_A = [(8.30, 121.00), (8.31, 121.00), (8.30, 121.05), (8.31, 121.05)]
BRIDGE = [(8.305, 121.70)]
GROUP_B = [(8.30, 122.40), (8.31, 122.40), (8.30, 122.45), (8.31, 122.45)]
MAX_SIZE = 6


@pytest.fixture(autouse=True)
def nucleus():
    set_ps2d_config('15N')


def _parts():
    integrator = VoigtIntegrator()
    cluster = GROUP_A + BRIDGE + GROUP_B
    return integrator, cluster, integrator._split_cluster_at_longest_link(
        cluster, max_size=MAX_SIZE)


def _touching_parts(integrator, part):
    """The groups of `part` joined by touching peaks."""
    links = [(i, j) for i in range(len(part)) for j in range(i + 1, len(part))
             if integrator._check_peaks_overlap(part[i], part[j])]
    return integrator._connected_parts(len(part), links)


def test_a_cluster_within_the_cap_is_left_whole():
    integrator = VoigtIntegrator()
    cluster = GROUP_A + BRIDGE
    assert integrator._split_cluster_at_longest_link(cluster, max_size=MAX_SIZE) == [cluster]


def test_no_part_is_larger_than_the_cap_and_every_peak_is_in_one():
    _, cluster, parts = _parts()
    assert max(len(part) for part in parts) <= MAX_SIZE
    assert sorted(peak for part in parts for peak in part) == sorted(cluster)


def test_every_part_is_connected_through_touching_peaks():
    integrator, _, parts = _parts()
    for part in parts:
        assert len(_touching_parts(integrator, part)) == 1


def test_no_part_holds_peaks_of_both_groups():
    _, _, parts = _parts()
    for part in parts:
        assert not (set(part) & set(GROUP_A) and set(part) & set(GROUP_B))

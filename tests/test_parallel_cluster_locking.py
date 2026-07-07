# ABOUTME: Regression tests for parallel-mode locked-cluster assignment matching.
# ABOUTME: Parallel series runs must resolve locked clusters like sequential mode does.

import numpy as np

from lunaNMR.core.parallel_voigt_processor import ParallelVoigtProcessor


def _proc():
    # _convert_locked_clusters_to_positions uses no instance state, so skip __init__.
    return ParallelVoigtProcessor.__new__(ParallelVoigtProcessor)


_CONTEXT = [
    {'assignment': '3', 'x_ppm': 8.33, 'y_ppm': 128.4},
    {'assignment': '4', 'x_ppm': 8.67, 'y_ppm': 122.9},
    {'assignment': '5', 'x_ppm': 9.04, 'y_ppm': 125.4},
]


def test_numeric_assignment_formats_match():
    # Locked clusters store raw DataFrame values (float/int/numpy); the current
    # spectrum's context keys are str(...). These must still match.
    proc = _proc()
    locked = [[3.0], [np.int64(4), 5.0]]  # '3.0', numpy int, float — as stored from a peak list
    clusters = proc._convert_locked_clusters_to_positions(locked, _CONTEXT)
    assert len(clusters) == 2
    assert clusters[0] == [(8.33, 128.4)]
    assert sorted(clusters[1]) == [(8.67, 122.9), (9.04, 125.4)]


def test_string_residue_labels_still_match():
    # Non-numeric assignments (real residue labels) must pass through unchanged.
    proc = _proc()
    context = [{'assignment': 'A12N-H', 'x_ppm': 1.0, 'y_ppm': 2.0}]
    clusters = proc._convert_locked_clusters_to_positions([['A12N-H']], context)
    assert clusters == [[(1.0, 2.0)]]


def test_already_matching_keys_unaffected():
    # The GUI supplies already-matching string keys; the fix must be a no-op there.
    proc = _proc()
    clusters = proc._convert_locked_clusters_to_positions([['3'], ['4']], _CONTEXT)
    assert len(clusters) == 2

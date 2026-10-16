# ABOUTME: The bound flags must survive the CSV boundary and gate the Kd input.
# ABOUTME: A peak at zero is data for I/I0 and not data for CSP; a bound is neither.

"""
Two states look identical in a series and need opposite handling.

A residue whose peak has DECAYED reports a near-zero intensity, and that zero is the
measurement - dropping it from an I/I0 series invents a plateau that is not there. A
residue whose FIT FAILED, pinned on a position or intensity bound, also reports a
near-zero-ish number, and that one carries no information at all: it came from a bound,
which for the intensity floor means it is a product of its neighbours' seeds.

Neither `success` nor R-squared separates them - R-squared is per cluster, so a member
with nothing at its position rides the one its healthy neighbours earn.
"""

import pandas as pd
import pytest

from lunaNMR.processors.multi_spectrum_processor import TIDY_COLUMNS


def test_the_tidy_csv_carries_the_flags_the_analysis_needs():
    """The flags exist on the fit result and used to stop at the CSV boundary, because
    every row is built by explicit key selection."""
    for column in ('position_pinned', 'intensity_pinned', 'detection_sigma'):
        assert column in TIDY_COLUMNS, f'{column} never reaches the analysis'


class TestUsablePoints:
    """`usable_points` is the filter both observables share, differing in one place."""

    @staticmethod
    def _frame(rows):
        return pd.DataFrame(rows)

    def test_a_bound_is_dropped_from_both_observables(self):
        from lunaNMR.utils.quality_categories import usable_points
        rows = self._frame([
            {'assignment': 'A1', 'height': 1.0e9, 'detection_sigma': 40.0,
             'position_pinned': False, 'intensity_pinned': False},
            {'assignment': 'A2', 'height': 4.0e6, 'detection_sigma': 0.4,
             'position_pinned': False, 'intensity_pinned': True},
        ])
        for observable in ('csp', 'intensity'):
            kept = usable_points(rows, observable)
            assert list(kept.assignment) == ['A1'], (
                f'{observable}: a value taken from a bound is not a measurement')

    def test_a_pinned_position_is_dropped_for_csp_and_kept_for_intensity(self):
        """A position on its margin is not a shift anyone measured. The height does not
        depend on it: on a 600 MHz hetNOE pair, widening the 15N margin from 0.020 to 0.05
        ppm moved the freed peaks' heights by a median 0.0%."""
        from lunaNMR.utils.quality_categories import usable_points
        rows = self._frame([
            {'assignment': 'A1', 'height': 1.0e9, 'detection_sigma': 40.0,
             'position_pinned': False, 'intensity_pinned': False},
            {'assignment': 'A2', 'height': 7.0e8, 'detection_sigma': 30.0,
             'position_pinned': True, 'intensity_pinned': False},
        ])
        assert list(usable_points(rows, 'intensity').assignment) == ['A1', 'A2']
        assert list(usable_points(rows, 'csp').assignment) == ['A1']

    def test_a_decayed_peak_is_kept_for_intensity_and_dropped_for_csp(self):
        """The distinction the whole filter exists for."""
        from lunaNMR.utils.quality_categories import usable_points
        rows = self._frame([
            {'assignment': 'A1', 'height': 1.0e9, 'detection_sigma': 40.0,
             'position_pinned': False, 'intensity_pinned': False},
            {'assignment': 'A2', 'height': 2.0e6, 'detection_sigma': 0.8,
             'position_pinned': False, 'intensity_pinned': False},
        ])
        assert list(usable_points(rows, 'intensity').assignment) == ['A1', 'A2'], (
            'zero IS the measurement for I/I0 - dropping it invents a plateau')
        assert list(usable_points(rows, 'csp').assignment) == ['A1'], (
            'a peak that is gone has no position to measure a shift from')

    def test_an_unknown_flag_is_not_read_as_permission(self):
        """A run from before the flags existed carries NaN. It must not be treated as
        'not pinned' silently - it is treated as usable, but only because refusing
        every legacy row would be worse; the caller is told with a count."""
        from lunaNMR.utils.quality_categories import usable_points
        rows = self._frame([{'assignment': 'A1', 'height': 1.0e9,
                             'detection_sigma': float('nan'),
                             'position_pinned': None, 'intensity_pinned': None}])
        kept = usable_points(rows, 'intensity')
        assert len(kept) == 1

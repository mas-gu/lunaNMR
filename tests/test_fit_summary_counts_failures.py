# ABOUTME: The summary of a Fit Spectrum run counts a peak as fitted only when its fit succeeded,
# ABOUTME: and distributes the fits over the labels the one quality rule gave them.

import pytest

from lunaNMR.processors.single_spectrum_processor import SingleSpectrumProcessor


def _summary(results):
    processor = SingleSpectrumProcessor.__new__(SingleSpectrumProcessor)
    processor.stats = {}
    return processor.get_processing_summary(results, len(results))


def _fit(quality, r_squared):
    return {'assignment': quality, 'success': True, 'fitted': True, 'quality': quality,
            'fitting_quality': quality, 'r_squared': r_squared, 'avg_r_squared': r_squared}


FAILED = {'assignment': 'F', 'success': False, 'fitted': False, 'quality': 'Failed',
          'fitting_quality': 'Failed', 'r_squared': 0.0, 'height': 0.0, 'volume': 0.0}


def test_a_peak_whose_fit_failed_is_counted_as_failed():
    """A failed peak keeps a row in the results, so the rows are not the successes."""
    summary = _summary([_fit('Excellent', 0.95), FAILED, _fit('Good', 0.85)])
    assert summary['successful_peaks'] == 2
    assert summary['failed_peaks'] == 1
    assert summary["success_rate"] == pytest.approx(100 * 2 / 3)


def test_the_average_r_squared_is_over_the_fits_that_succeeded():
    summary = _summary([_fit('Excellent', 0.95), FAILED, _fit('Good', 0.85)])
    assert summary['average_r_squared'] == (0.95 + 0.85) / 2


def test_the_distribution_counts_the_label_each_fit_carries():
    """The labels are the fits' own, from the one quality rule: the summary applies no
    thresholds of its own (0.92 is Excellent and 0.7 Fair there)."""
    summary = _summary([_fit('Excellent', 0.92), _fit('Good', 0.87), _fit('Fair', 0.7),
                        _fit('Poor', 0.3), FAILED])
    assert summary['quality_distribution'] == {
        'Excellent': 1, 'Good': 1, 'Fair': 1, 'Poor': 1, 'Failed': 1}

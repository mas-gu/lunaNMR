# ABOUTME: Tests that a singular covariance is captured as a per-fit field, not suppressed.
# ABOUTME: A different OptimizeWarning at the same call site must still reach the log.

import sys
import warnings
from pathlib import Path

import pytest

_DIR = Path(__file__).resolve().parent.parent / "dynamiXs_Kd"
sys.path.insert(0, str(_DIR))

from scipy.optimize import OptimizeWarning

from kd_fit import _singular_covariance


def _caught(category, message):
    """One recorded warning, shaped as catch_warnings(record=True) yields them."""
    with warnings.catch_warnings(record=True) as recorded:
        warnings.simplefilter("always")
        warnings.warn(message, category)
    return list(recorded)


class TestASingularCovarianceIsRecognised:
    """A null Kd_err has several possible causes; only this distinguishes which one."""

    def test_the_covariance_message_sets_the_flag(self):
        caught = _caught(OptimizeWarning,
                         "Covariance of the parameters could not be estimated")
        assert _singular_covariance(caught) is True

    def test_nothing_caught_means_not_singular(self):
        assert _singular_covariance([]) is False

    def test_the_match_is_not_case_sensitive_on_the_leading_letter(self):
        """scipy has emitted both 'Covariance' and 'covariance' across versions."""
        assert _singular_covariance(
            _caught(OptimizeWarning, "covariance of the parameters could not be estimated"))


class TestADifferentWarningStillReachesTheLog:
    """Suppression is how you stop hearing about a real problem later. Anything at these
    call sites that is NOT the expected singular covariance must still surface."""

    def test_another_optimize_warning_is_re_emitted(self):
        caught = _caught(OptimizeWarning, "some entirely different optimiser complaint")
        with warnings.catch_warnings(record=True) as out:
            warnings.simplefilter("always")
            singular = _singular_covariance(caught)
        assert singular is False
        assert [str(w.message) for w in out] == ["some entirely different optimiser complaint"]

    def test_a_non_optimize_warning_is_re_emitted(self):
        caught = _caught(RuntimeWarning, "overflow encountered in exp")
        with warnings.catch_warnings(record=True) as out:
            warnings.simplefilter("always")
            _singular_covariance(caught)
        assert [w.category for w in out] == [RuntimeWarning]

    def test_a_mixed_batch_keeps_both_behaviours(self):
        """The expected one is absorbed into the flag; the unexpected one still shows."""
        caught = (_caught(OptimizeWarning, "Covariance of the parameters could not be estimated")
                  + _caught(RuntimeWarning, "invalid value encountered"))
        with warnings.catch_warnings(record=True) as out:
            warnings.simplefilter("always")
            singular = _singular_covariance(caught)
        assert singular is True
        assert [str(w.message) for w in out] == ["invalid value encountered"]

    def test_the_re_emitted_warning_keeps_its_origin(self):
        """warn_explicit preserves the original file and line, so the log still points
        at the real source rather than at this helper."""
        caught = _caught(OptimizeWarning, "a different complaint")
        origin = (caught[0].filename, caught[0].lineno)
        with warnings.catch_warnings(record=True) as out:
            warnings.simplefilter("always")
            _singular_covariance(caught)
        assert (out[0].filename, out[0].lineno) == origin

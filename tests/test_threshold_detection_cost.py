# ABOUTME: S/N threshold detection costs one pass over the spectrum, not one per candidate peak.
# ABOUTME: Real DNAJA1 spectrum, at the low S/N an automated batch run detects with.

"""
The noise level and each candidate's region were recomputed over the whole spectrum for
every candidate: at S/N 2.5 the DNAJA1 spectrum has ~15000 candidates, and one detection
took over two minutes. `lunaNMR batch` detects five or six times per spectrum.
"""

import time
from pathlib import Path

import pytest

DNAJA1 = Path(__file__).resolve().parents[2] / "data_test" / "data_binding" / "DNAJA1_HSPA8"

pytestmark = pytest.mark.skipif(not DNAJA1.exists(), reason="data_test not present")


@pytest.fixture(scope="module")
def integrator():
    from lunaNMR.core.core_integrator import EnhancedVoigtIntegrator
    it = EnhancedVoigtIntegrator()
    assert it.load_spectrum_only(str(DNAJA1 / "2D_DNAJA1_0o0.ft"))
    return it


def test_the_noise_is_estimated_once_per_detection(integrator, monkeypatch):
    calls = []
    estimate = integrator._estimate_noise_level_advanced
    monkeypatch.setattr(integrator, "_estimate_noise_level_advanced",
                        lambda: calls.append(1) or estimate())
    peaks = integrator._detect_peaks_by_threshold(5 * estimate())
    assert len(peaks) > 100
    assert len(calls) == 1


def test_a_low_threshold_detection_takes_seconds(integrator):
    start = time.time()
    peaks = integrator._detect_peaks_by_threshold(2.5 * integrator._estimate_noise_level_advanced())
    assert len(peaks) > 10000
    assert time.time() - start < 30

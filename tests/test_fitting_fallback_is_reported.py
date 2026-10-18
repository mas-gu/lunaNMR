# ABOUTME: A fit that falls back to a different algorithm says so, on every result and in a warning.
# ABOUTME: Real DNAJA1 reference spectrum; the cluster pipeline is given a cluster it cannot read.

"""
When the cluster pipeline raised, every peak was refitted on its own - no clusters,
territories or drift limits - and the only trace was a traceback on stderr: the results'
processing_mode was overwritten downstream and the run reported success. The fallback
stays; it is now a warning that names it, a mark on each result, and a field the series
summary reports.
"""

from pathlib import Path

import pytest

DNAJA1 = Path(__file__).resolve().parents[2] / "data_test" / "data_binding" / "DNAJA1_HSPA8"

pytestmark = pytest.mark.skipif(not DNAJA1.exists(), reason="data_test not present")

# A locked cluster must be a list of assignments; an integer is one the pipeline cannot read
UNREADABLE_CLUSTERS = [5]


@pytest.fixture(scope="module")
def fallback_run():
    """(fitted results, captured stdout) of a fit whose cluster pipeline raised."""
    import contextlib
    import io
    from lunaNMR.core.core_integrator import EnhancedVoigtIntegrator
    from lunaNMR.processors.single_spectrum_processor import SingleSpectrumProcessor
    from lunaNMR.utils.file_manager import NMRFileManager
    from lunaNMR.utils.parameter_manager import NMRParameterManager

    integrator = EnhancedVoigtIntegrator()
    integrator.load_nmr_file(str(DNAJA1 / "2D_DNAJA1_0o0.ft"))
    peaks = NMRFileManager().load_peak_list(str(DNAJA1 / "DNAJA1.txt")).head(6)
    integrator.peak_list = peaks
    processor = SingleSpectrumProcessor(integrator, NMRParameterManager())
    out = io.StringIO()
    with contextlib.redirect_stdout(out):
        results, _ = processor.process_peak_list(
            peaks, {'use_parallel': False}, locked_clusters_by_assignment=UNREADABLE_CLUSTERS)
    return results, out.getvalue()


def test_every_result_is_marked_with_the_fallback(fallback_run):
    results, _ = fallback_run
    assert len(results) == 6
    assert all(r['processing_mode'] == 'sequential_fallback' for r in results)
    assert all('TypeError' in r['fitting_fallback'] for r in results)


def test_a_warning_names_the_fallback_and_its_reason(fallback_run):
    _, stdout = fallback_run
    warnings = [line for line in stdout.splitlines() if line.startswith('Warning:')]
    assert any('per-peak' in line and 'TypeError' in line for line in warnings), warnings


def test_the_spectrum_reports_which_fallbacks_produced_its_results(fallback_run):
    from lunaNMR.processors.multi_spectrum_processor import MultiSpectrumProcessor
    results, _ = fallback_run
    reported = MultiSpectrumProcessor._fitting_fallbacks(results)
    assert len(reported) == 1 and 'TypeError' in reported[0]


def test_a_normal_fit_reports_no_fallback():
    from lunaNMR.processors.multi_spectrum_processor import MultiSpectrumProcessor
    assert MultiSpectrumProcessor._fitting_fallbacks([{'processing_mode': 'sequential'}, None]) == []


def test_the_series_summary_reports_fallbacks_per_spectrum(dnaja1_titration_series):
    import json
    summary = json.loads((dnaja1_titration_series.parent / "summary.json").read_text())
    summary = summary.get("result", summary)
    assert summary["fitting_fallbacks"] == {}

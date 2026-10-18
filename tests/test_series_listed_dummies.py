# ABOUTME: A dummy that is already in a series' peak list is fitted once per spectrum, not twice.
# ABOUTME: The series dialog's Independent mode hands over the GUI's fitted peaks, dummies included.

"""
Independent mode sends the main window's fitted peaks, which hold the dummy Fit All Peaks
added, as the series' peak list. The first spectrum's detection matched that dummy again
under its own name, the series stored it as a dummy of its own, and every later spectrum
fitted it twice: a two-member cluster of one peak and its copy, sharing its intensity.
"""

import contextlib
import io
import warnings
from pathlib import Path

import pandas as pd
import pytest

DNAJA1 = Path(__file__).resolve().parents[2] / "data_test" / "data_binding" / "DNAJA1_HSPA8"

pytestmark = pytest.mark.skipif(not DNAJA1.exists(), reason="data_test not present")

SPECTRA = ["2D_DNAJA1_0o0.ft", "2D_DNAJA1_HSPA8_0o125.ft"]


def _gui_fitted_peaks():
    """What Fit All Peaks leaves in integrator.fitted_peaks on the first spectrum."""
    from lunaNMR.core.core_integrator import EnhancedVoigtIntegrator
    from lunaNMR.utils.file_manager import NMRFileManager
    integrator = EnhancedVoigtIntegrator()
    integrator.load_nmr_file(str(DNAJA1 / SPECTRA[0]))
    integrator.peak_list = NMRFileManager().load_peak_list(str(DNAJA1 / "DNAJA1.txt"))
    integrator.set_search_window(0.07, 0.07)
    return integrator._detect_peaks_reference_based()


@pytest.fixture(scope="module")
def independent_series(tmp_path_factory):
    from numba.core.errors import NumbaWarning
    from lunaNMR.gui.dialogs.series_integration_dialog import SeriesIntegrationDialog
    from lunaNMR.processors.multi_spectrum_processor import MultiSpectrumProcessor
    from lunaNMR.cli import _default_series_params
    with contextlib.redirect_stdout(io.StringIO()), warnings.catch_warnings():
        warnings.simplefilter("ignore", NumbaWarning)
        fitted = _gui_fitted_peaks()
        reference = SeriesIntegrationDialog._convert_fitted_peaks_to_dataframe(None, fitted)
        assert reference["Assignment"].astype(str).str.startswith("dummy_").any()
        processor = MultiSpectrumProcessor(_default_series_params())
        result = processor.process_nmr_series(
            [str(DNAJA1 / name) for name in SPECTRA], reference,
            str(tmp_path_factory.mktemp("listed_dummies") / "out"),
            peak_source_mode="independent")
    return reference, result


def test_every_spectrum_fits_each_assignment_once(independent_series):
    reference, result = independent_series
    for name, data in result.results.items():
        assignments = pd.Series([str(r["assignment"]) for r in data["fitted_results"] if r])
        assert not assignments.duplicated().any(), (name, assignments[assignments.duplicated()].tolist())
        assert len(assignments) == len(reference), name


def test_a_detection_dummy_does_not_take_a_listed_peaks_name():
    """A peak list holding dummy_001 does not get a second dummy_001 from detection."""
    from lunaNMR.core.core_integrator import EnhancedVoigtIntegrator
    from lunaNMR.utils.file_manager import NMRFileManager
    with contextlib.redirect_stdout(io.StringIO()):
        peaks = NMRFileManager().load_peak_list(str(DNAJA1 / "DNAJA1.txt"))
        # A listed peak under a dummy's name, far from the dummy detection adds
        far = peaks.iloc[[0]].assign(Assignment="dummy_001")
        integrator = EnhancedVoigtIntegrator()
        integrator.load_nmr_file(str(DNAJA1 / SPECTRA[0]))
        integrator.peak_list = pd.concat([peaks, far], ignore_index=True)
        integrator.peak_list.loc[len(peaks), "Position_X"] = 6.0
        integrator.peak_list.loc[len(peaks), "Position_Y"] = 100.0
        integrator.set_search_window(0.07, 0.07)
        detected = integrator._detect_peaks_reference_based()
    names = pd.Series([str(p["assignment"]) for p in detected])
    assert names.str.startswith("dummy_").sum() == 2
    assert not names.duplicated().any(), names[names.duplicated()].tolist()

# ABOUTME: A series is stored slim - fit parameters and region axes - and its fit-detail arrays rebuilt.
# ABOUTME: Real DNAJA1 series; both viewers that draw a fit must get the arrays back.

"""
Stored whole, every fit's surfaces and region grids make batch_results.json ~80 MB for
a 7-point titration, in the run folder and again in every bundle. The parameters
(`all_peaks`) and the region's own ppm axes are enough to rebuild them from the
spectrum, as single-spectrum fit results already are. What matters is that the views
which draw a fit still get one.
"""

import copy
import json
import os
import sys
from pathlib import Path

import numpy as np
import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

DNAJA1 = Path(__file__).resolve().parents[2] / "data_test" / "data_binding" / "DNAJA1_HSPA8"
REFERENCE = DNAJA1 / "2D_DNAJA1_0o0.ft"
SURFACE_KEYS = ("fitted_2d_surface", "individual_surfaces")

pytestmark = pytest.mark.skipif(not DNAJA1.exists(), reason="data_test not present")


@pytest.fixture(scope="module")
def series_dir(tmp_path_factory):
    from lunaNMR.cli import main
    out = tmp_path_factory.mktemp("series") / "DNAJA1"
    assert main(["series", "--spectra", str(DNAJA1), "--peaks", str(DNAJA1 / "DNAJA1.txt"),
                 "--out", str(out), "--mode", "titration", "--peak-source", "cascade"]) == 0
    return out


@pytest.fixture(scope="module")
def stored(series_dir):
    return json.loads((series_dir / "batch_results.json").read_text())


@pytest.fixture(scope="module")
def reference_rows(stored):
    """The reference spectrum's per-peak rows, as stored."""
    result = copy.deepcopy(stored["results"][REFERENCE.name])
    return result["integration_results"]


@pytest.fixture(scope="module")
def integrator():
    from lunaNMR.core.core_integrator import EnhancedVoigtIntegrator
    it = EnhancedVoigtIntegrator()
    it.load_nmr_file(str(REFERENCE))
    return it


@pytest.fixture(scope="module")
def app():
    from PySide6.QtWidgets import QApplication
    return QApplication.instance() or QApplication(sys.argv)


def _fitted(rows):
    return [r for r in rows if r.get("all_peaks")]


def test_the_stored_series_is_small(series_dir):
    size = (series_dir / "batch_results.json").stat().st_size
    assert size < 5e6, f"{size / 1e6:.1f} MB"


def test_no_surface_or_region_grid_is_stored(stored):
    for result in stored["results"].values():
        for key in ("integration_results", "fitted_results"):
            for row in result.get(key, []):
                assert not any(k in row for k in SURFACE_KEYS)
                region = row.get("region_2d")
                if region is not None:
                    assert set(region) == {"f1_ppm", "f2_ppm"}


def test_a_fit_keeps_what_rebuilds_it(reference_rows):
    rows = _fitted(reference_rows)
    assert rows
    for row in rows:
        assert row["region_2d"]["f1_ppm"] and row["region_2d"]["f2_ppm"]


def test_the_region_rebuilds_from_the_spectrum(reference_rows, integrator):
    """Exactly the stored axes, and the spectrum's own values on them."""
    from lunaNMR.utils.project_manager import rebuild_fit_arrays
    row = dict(_fitted(reference_rows)[0])
    assert rebuild_fit_arrays(row, integrator)
    region = row["region_2d"]
    f1, f2 = np.asarray(region["f1_ppm"]), np.asarray(region["f2_ppm"])
    assert region["intensity"].shape == (f1.size, f2.size)
    y0 = int(np.argmin(np.abs(integrator.ppm_y_axis - f1[0])))
    x0 = int(np.argmin(np.abs(integrator.ppm_x_axis - f2[0])))
    np.testing.assert_array_equal(
        region["intensity"], integrator.nmr_data[y0:y0 + f1.size, x0:x0 + f2.size])
    np.testing.assert_allclose(integrator.ppm_y_axis[y0:y0 + f1.size], f1)


def test_the_surface_rebuilds_at_the_fitted_height(reference_rows, integrator):
    from lunaNMR.utils.project_manager import rebuild_fit_arrays
    row = dict(max(_fitted(reference_rows), key=lambda r: abs(r.get("height") or 0)))
    assert rebuild_fit_arrays(row, integrator)
    assert row["fitted_2d_surface"].shape == row["region_2d"]["intensity"].shape
    target = next((p for p in row["all_peaks"]
                   if p.get("assignment") == row.get("assignment")), row["all_peaks"][0])
    own = np.asarray(row["individual_surfaces"][row["all_peaks"].index(target)])
    assert np.max(np.abs(own)) == pytest.approx(abs(target["height"]), rel=0.05)


def test_a_rebuild_without_a_spectrum_leaves_the_row_alone(reference_rows):
    from types import SimpleNamespace
    from lunaNMR.utils.project_manager import rebuild_fit_arrays
    row = dict(_fitted(reference_rows)[0])
    assert not rebuild_fit_arrays(row, SimpleNamespace(nmr_data=None))
    assert "fitted_2d_surface" not in row


def test_the_overlay_viewer_rebuilds_a_fit(app, stored):
    """The overlay's fit-detail lookup returns a result it can draw."""
    from lunaNMR.gui.dialogs.multi_spectrum_viewer_dialog import MultiSpectrumViewerDialog
    result = copy.deepcopy(stored["results"][REFERENCE.name])
    result["spectrum_name"] = REFERENCE.name
    result["spectrum_file"] = str(REFERENCE)
    dialog = MultiSpectrumViewerDialog(all_results=[result])
    try:
        spec = dialog.spectra[0]
        assignment = _fitted(result["integration_results"])[0]["assignment"]
        found = dialog._find_stored_voigt_result(spec, assignment)
        assert found is not None
        assert found["region_2d"]["intensity"].ndim == 2
        assert found["fitted_2d_surface"].shape == found["region_2d"]["intensity"].shape
    finally:
        dialog.deleteLater()


def test_the_spectrum_viewer_rebuilds_a_fit(app, stored):
    from lunaNMR.gui.dialogs.spectrum_viewer_dialog import SpectrumViewerDialog
    result = copy.deepcopy(stored["results"][REFERENCE.name])
    dialog = SpectrumViewerDialog(spectrum_name=REFERENCE.name,
                                  spectrum_file_path=str(REFERENCE), result_data=result)
    try:
        assignment = _fitted(result["integration_results"])[0]["assignment"]
        found = dialog._find_stored_voigt_result(assignment)
        assert found is not None
        assert found["region_2d"]["intensity"].ndim == 2
        assert found["fitted_2d_surface"].shape == found["region_2d"]["intensity"].shape
    finally:
        dialog.deleteLater()


def test_the_overlay_reports_no_fit_when_the_spectrum_is_gone(app, stored, tmp_path):
    """A fit that cannot be rebuilt is 'no stored result', not a KeyError mid-plot."""
    from lunaNMR.gui.dialogs.multi_spectrum_viewer_dialog import MultiSpectrumViewerDialog
    result = copy.deepcopy(stored["results"][REFERENCE.name])
    result["spectrum_name"] = REFERENCE.name
    result["spectrum_file"] = str(tmp_path / "moved_away.ft")
    dialog = MultiSpectrumViewerDialog(all_results=[result])
    try:
        assignment = _fitted(result["integration_results"])[0]["assignment"]
        assert dialog._find_stored_voigt_result(dialog.spectra[0], assignment) is None
    finally:
        dialog.deleteLater()


def test_the_spectrum_viewer_reports_no_fit_when_the_spectrum_is_gone(app, stored, tmp_path):
    from lunaNMR.gui.dialogs.spectrum_viewer_dialog import SpectrumViewerDialog
    result = copy.deepcopy(stored["results"][REFERENCE.name])
    dialog = SpectrumViewerDialog(spectrum_name=REFERENCE.name,
                                  spectrum_file_path=str(tmp_path / "moved_away.ft"),
                                  result_data=result)
    try:
        assignment = _fitted(result["integration_results"])[0]["assignment"]
        assert dialog._find_stored_voigt_result(assignment) is None
    finally:
        dialog.deleteLater()


def _rebuilt_series(stored, integrator):
    """The stored series with the reference spectrum's fits carrying their arrays, as a
    series just fitted in the GUI does."""
    from lunaNMR.utils.project_manager import batch_results_from_dict, rebuild_fit_arrays
    batch = batch_results_from_dict(copy.deepcopy(stored))
    for row in batch.results[REFERENCE.name]["integration_results"]:
        rebuild_fit_arrays(row, integrator)
    return batch


@pytest.mark.parametrize("embed, kept", [(True, True), (False, False)])
def test_embedding_region_data_keeps_a_series_whole(stored, integrator, tmp_path, embed, kept):
    """'Embed raw fit-region data' exists so a bundle opens without its spectra."""
    import types
    from lunaNMR.utils.project_manager import ProjectManager
    window = types.SimpleNamespace(
        integrator=types.SimpleNamespace(peak_list=None),
        saved_series={"DNAJA1": _rebuilt_series(stored, integrator)})
    bundle = tmp_path / "b.lunaNMR"
    assert ProjectManager(window).save_project(bundle, embed_region_data=embed)[0]
    saved = json.loads((bundle / "series_results" / "DNAJA1" / "batch_results.json").read_text())
    rows = _fitted(saved["results"][REFERENCE.name]["integration_results"])
    assert ("intensity" in rows[0]["region_2d"]) is kept
    assert ("fitted_2d_surface" in rows[0]) is kept

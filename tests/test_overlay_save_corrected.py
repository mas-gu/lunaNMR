# ABOUTME: Tests the overlay viewer's refit + 'Save corrected' on real series output written by `lunaNMR series`.
# ABOUTME: A correction must land in its own point's columns of both CSVs, with the refit's quality and bound flags.

import contextlib
import io
import json
import math
import shutil
import sys
from pathlib import Path

import pandas as pd
import pytest

from PySide6.QtWidgets import QApplication

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "modules" / "dynamiXs_v2o0"))

DNAJA1 = Path(__file__).resolve().parents[2] / "data_test" / "data_binding" / "DNAJA1_HSPA8"


@pytest.fixture(scope="module")
def app():
    return QApplication.instance() or QApplication(sys.argv)


def _spectrum_index(viewer, name):
    return [s['name'] for s in viewer.spectra].index(name)


def _refit_in_place(viewer, spec_idx, assignment):
    """Refit one peak at the position the series left it, as Re-fit does for it."""
    cx, cy = viewer._peak_center(viewer._find_peak_obj(spec_idx, assignment))
    r2 = viewer._refit_peak_in_spectrum(spec_idx, assignment, cx, cy)
    assert r2 is not None
    viewer._edited.add((spec_idx, assignment))
    return r2


@pytest.fixture
def titration(app, dnaja1_titration_series, tmp_path, open_overlay_viewer):
    """A writable copy of the 3-point titration output, and the viewer open on it."""
    folder = tmp_path / "out"
    shutil.copytree(dnaja1_titration_series, folder)
    viewer = open_overlay_viewer(folder)
    yield folder, viewer
    viewer.deleteLater()


@pytest.fixture(scope="module")
def repeat_series(tmp_path_factory):
    """Output of a time series with a repeat acquisition named by the trailing-letter
    convention (10ms, 10msb, 20ms), run on real spectra through `lunaNMR series`."""
    if not DNAJA1.exists():
        pytest.skip("data_test not present")
    from lunaNMR.cli import main
    root = tmp_path_factory.mktemp("repeat_series")
    (root / "spectra").mkdir()
    for point, name in (("0o0", "s_10ms.ft"), ("0o5", "s_10msb.ft"), ("3o0", "s_20ms.ft")):
        (root / "spectra" / name).symlink_to(next(DNAJA1.glob(f"2D_*_{point}.ft")))
    with contextlib.redirect_stdout(io.StringIO()):
        code = main(["series", "--spectra", str(root / "spectra"), "--peaks", str(DNAJA1 / "DNAJA1.txt"),
                     "--mode", "time", "--out", str(root / "out"), "--format", "json"])
    assert code == 0
    return root / "out"


class TestCorrectionLandsInItsOwnPoint:
    def test_a_repeat_acquisition_saves_to_the_column_the_series_gave_it(self, app, repeat_series, tmp_path,
                                                                         open_overlay_viewer):
        folder = tmp_path / "out"
        shutil.copytree(repeat_series, folder)
        columns = {c['spectrum']: c['column']
                   for c in json.loads((folder / "series_metadata.json").read_text())['columns']}
        assert columns['s_10msb.ft'] == '10_2'
        viewer = open_overlay_viewer(folder)
        try:
            idx = _spectrum_index(viewer, 's_10msb.ft')
            peak = viewer._find_peak_obj(idx, '3LysH')
            peak['height'] = 999.0
            viewer._edited.add((idx, '3LysH'))
            viewer._on_save_corrected()
        finally:
            viewer.deleteLater()

        tracking = pd.read_csv(folder / "comprehensive_peak_tracking.csv")
        row = tracking[tracking['Assignment'] == '3LysH'].iloc[0]
        assert [c for c in tracking.columns if c.endswith('_Height') and row[c] == 999.0] == ['10_2_Height']
        tidy = pd.read_csv(folder / "series_analysis_tidy.csv", dtype={'spectrum_name': str})
        assert tidy.loc[tidy['height'] == 999.0, 'spectrum_name'].tolist() == ['10_2']

    def test_a_series_without_its_column_map_is_not_saved(self, titration):
        folder, viewer = titration
        (folder / "series_metadata.json").unlink()
        before = (folder / "comprehensive_peak_tracking.csv").read_text()
        idx = _spectrum_index(viewer, '2D_DNAJA1_HSPA8_3o0.ft')
        viewer._find_peak_obj(idx, '3LysH')['height'] = 999.0
        viewer._edited.add((idx, '3LysH'))
        viewer._on_save_corrected()
        assert (folder / "comprehensive_peak_tracking.csv").read_text() == before
        assert 'series_metadata.json' in viewer.edit_status_label.text()


class TestTitrationTidyRows:
    def test_a_titration_correction_reaches_the_tidy_csv(self, titration):
        folder, viewer = titration
        idx = _spectrum_index(viewer, '2D_DNAJA1_HSPA8_3o0.ft')
        _refit_in_place(viewer, idx, '42GlyH')
        viewer._on_save_corrected()

        tidy = pd.read_csv(folder / "series_analysis_tidy.csv", dtype={'spectrum_name': str})
        row = tidy[(tidy['spectrum_name'] == '3') & (tidy['assignment'] == '42GlyH')].iloc[0]
        assert math.isclose(row['height'], viewer._find_peak_obj(idx, '42GlyH')['height'], rel_tol=1e-9)

    def test_the_point_labels_survive_the_rewrite(self, titration):
        folder, viewer = titration
        idx = _spectrum_index(viewer, '2D_DNAJA1_HSPA8_3o0.ft')
        _refit_in_place(viewer, idx, '42GlyH')
        viewer._on_save_corrected()
        tidy = pd.read_csv(folder / "series_analysis_tidy.csv", dtype={'spectrum_name': str})
        assert sorted(set(tidy['spectrum_name'])) == ['0', '0.5', '3']


class TestRefitQualityRule:
    def test_a_refit_below_the_failed_threshold_is_saved_failed(self, titration):
        from dynamiXs_Kd.kd_input import load_titration_tracking
        folder, viewer = titration
        idx = _spectrum_index(viewer, '2D_DNAJA1_HSPA8_3o0.ft')
        r2 = _refit_in_place(viewer, idx, '28LysH')
        assert r2 < 0.2
        viewer._on_save_corrected()

        tracking = pd.read_csv(folder / "comprehensive_peak_tracking.csv")
        assert tracking.loc[tracking['Assignment'] == '28LysH', '3_Quality'].item() == 'Failed'
        points, residues = load_titration_tracking(str(folder / "comprehensive_peak_tracking.csv"))
        assert math.isnan(residues['28LysH']['height'][points.index(3.0)])


class TestRefitFlags:
    def test_a_point_refitted_off_its_intensity_floor_is_no_longer_pinned(self, titration):
        """48IleH on the reference point: the series leaves it on its intensity floor at
        1.1 sigma; refitted in place it comes back off the floor at ~3.9 sigma with a
        positive height. The reference point is fitted before any linewidth cap exists,
        so this case does not move with the cap."""
        from dynamiXs_Kd.kd_input import load_titration_tidy, load_titration_tracking
        folder, viewer = titration
        idx = _spectrum_index(viewer, '2D_DNAJA1_0o0.ft')
        assert viewer._find_peak_obj(idx, '48IleH')['intensity_pinned'] is True
        _refit_in_place(viewer, idx, '48IleH')
        peak = viewer._find_peak_obj(idx, '48IleH')
        assert peak['intensity_pinned'] is False
        assert peak['height'] > 0
        viewer._on_save_corrected()

        tracking = pd.read_csv(folder / "comprehensive_peak_tracking.csv")
        assert tracking.loc[tracking['Assignment'] == '48IleH', '0_Intensity_Pinned'].item() == 0.0
        tidy = pd.read_csv(folder / "series_analysis_tidy.csv", dtype={'spectrum_name': str})
        row = tidy[(tidy['spectrum_name'] == '0') & (tidy['assignment'] == '48IleH')].iloc[0]
        assert not row['intensity_pinned']
        assert math.isclose(row['detection_sigma'], peak['detection_sigma'], rel_tol=1e-9)
        for loader, name in ((load_titration_tracking, "comprehensive_peak_tracking.csv"),
                             (load_titration_tidy, "series_analysis_tidy.csv")):
            points, residues = loader(str(folder / name))
            assert math.isclose(residues['48IleH']['height'][points.index(0.0)], peak['height'],
                                rel_tol=1e-6)

    def test_the_detection_of_a_refit_is_its_own(self, titration):
        folder, viewer = titration
        idx = _spectrum_index(viewer, '2D_DNAJA1_HSPA8_3o0.ft')
        _refit_in_place(viewer, idx, '48IleH')
        sigma = viewer._find_peak_obj(idx, '48IleH')['detection_sigma']
        assert sigma < 3.0
        viewer._on_save_corrected()

        tracking = pd.read_csv(folder / "comprehensive_peak_tracking.csv")
        row = tracking[tracking['Assignment'] == '48IleH'].iloc[0]
        assert row['3_Detected'] == 0.0
        assert math.isclose(row['3_SNR'], sigma, rel_tol=1e-9)


def _metadata(folder):
    return json.loads((folder / "series_metadata.json").read_text())


def _rewrite_metadata(folder, edit):
    meta = _metadata(folder)
    edit(meta)
    (folder / "series_metadata.json").write_text(json.dumps(meta))


@contextlib.contextmanager
def _global_ps2d_config(nucleus):
    """The global PS2D config holds `nucleus`'s defaults - as after the main window
    loaded such a spectrum - and is put back afterwards."""
    from lunaNMR.core.ps2d_config import PS2DConfig, get_ps2d_config
    config = get_ps2d_config()
    saved = dict(vars(config))
    vars(config).update(vars(PS2DConfig(nucleus)))
    try:
        yield config
    finally:
        vars(config).clear()
        vars(config).update(saved)


class TestRefitFitsAsTheSeriesDid:
    def test_a_time_series_records_no_territories(self, repeat_series):
        meta = _metadata(repeat_series)
        assert meta['ps2d_config']['nucleus_type'] == '15N'
        assert not any('territories' in c for c in meta['columns'])

    def test_a_refit_uses_the_recorded_config_and_leaves_the_global_one(self, titration,
                                                                         open_overlay_viewer):
        folder, viewer = titration
        idx = _spectrum_index(viewer, '2D_DNAJA1_HSPA8_3o0.ft')
        _refit_in_place(viewer, idx, '42GlyH')
        expected = viewer._find_peak_obj(idx, '42GlyH')['height']

        other = open_overlay_viewer(folder)
        try:
            with _global_ps2d_config('13C') as config:
                carbon = dict(vars(config))
                _refit_in_place(other, idx, '42GlyH')
                assert dict(vars(config)) == carbon
            assert math.isclose(other._find_peak_obj(idx, '42GlyH')['height'], expected,
                                rel_tol=1e-9)
        finally:
            other.deleteLater()

    def test_a_refit_keeps_the_linewidth_caps_of_its_point(self, titration):
        """Caps are (max total F1, max total F2); lw_y/lw_x are the Gaussian parts,
        which a capped total bounds. Uncapped, 42GlyH at 3 eq fits 0.135/0.031 ppm."""
        folder, viewer = titration
        caps = (0.12, 0.028)

        def cap(meta):
            entry = next(c for c in meta['columns'] if c['column'] == '3')
            entry['linewidth_caps'] = {'42GlyH': list(caps)}
        _rewrite_metadata(folder, cap)
        idx = _spectrum_index(viewer, '2D_DNAJA1_HSPA8_3o0.ft')
        _refit_in_place(viewer, idx, '42GlyH')
        peak = viewer._find_peak_obj(idx, '42GlyH')
        assert peak['lw_y'] <= caps[0] and peak['lw_x'] <= caps[1]
        assert viewer.spectra[idx]['integrator'].titration_lw_margins is None

    def test_a_refit_without_a_recorded_config_says_which_it_used(self, titration):
        folder, viewer = titration
        _rewrite_metadata(folder, lambda meta: meta.pop('ps2d_config'))
        idx = _spectrum_index(viewer, '2D_DNAJA1_HSPA8_3o0.ft')
        viewer._pending_corrections[(idx, '42GlyH')] = viewer._peak_center(
            viewer._find_peak_obj(idx, '42GlyH'))
        viewer._on_refit_corrected()
        status = viewer.edit_status_label.text()
        assert "no PS2D config recorded" in status and "15N" in status

    def test_a_refit_with_the_recorded_config_says_nothing_of_it(self, titration):
        folder, viewer = titration
        idx = _spectrum_index(viewer, '2D_DNAJA1_HSPA8_3o0.ft')
        viewer._pending_corrections[(idx, '42GlyH')] = viewer._peak_center(
            viewer._find_peak_obj(idx, '42GlyH'))
        viewer._on_refit_corrected()
        assert "PS2D" not in viewer.edit_status_label.text()


METHYL = Path(__file__).resolve().parents[2] / "data_test" / "methyl" / "sampleNBD"


@pytest.fixture(scope="module")
def methyl_series(tmp_path_factory):
    """A real two-point 13C methyl T2 series, the spectra linked under neutral names."""
    if not METHYL.exists():
        pytest.skip("data_test not present")
    from lunaNMR.cli import main
    root = tmp_path_factory.mktemp("methyl_series")
    (root / "spectra").mkdir()
    for delay in ("3ms", "18ms"):
        (root / "spectra" / f"T2_{delay}.ft").symlink_to(METHYL / f"001_T2_NBD_{delay}.ft")
    with contextlib.redirect_stdout(io.StringIO()):
        assert main(["series", "--spectra", str(root / "spectra"),
                     "--peaks", str(METHYL / "001_NBDADP_NBD.txt"), "--mode", "time",
                     "--peak-source", "reference", "--out", str(root / "out")]) == 0
    return root / "out"


class TestARefitOfA13CSeriesUnderA15NWindow:
    """The main window on a 15N spectrum while a 13C series is refitted in the overlay:
    the refit must fit with the 13C config the series recorded, not the window's."""

    def test_the_refit_matches_a_refit_under_the_series_own_nucleus(self, app, methyl_series,
                                                                     open_overlay_viewer):
        assert _metadata(methyl_series)['ps2d_config']['nucleus_type'] == '13C'
        heights = {}
        for nucleus in ('13C', '15N'):
            viewer = open_overlay_viewer(methyl_series)
            try:
                with _global_ps2d_config(nucleus) as config:
                    before = dict(vars(config))
                    idx = _spectrum_index(viewer, 'T2_18ms.ft')
                    _refit_in_place(viewer, idx, '27_cg2_D')
                    assert dict(vars(config)) == before
                heights[nucleus] = viewer._find_peak_obj(idx, '27_cg2_D')['height']
            finally:
                viewer.deleteLater()
        assert math.isclose(heights['15N'], heights['13C'], rel_tol=1e-9)

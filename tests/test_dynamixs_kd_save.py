# ABOUTME: Tests project-save coverage for the methyl-T2 dynamiXs page and the Kd module.
# ABOUTME: Page session-state round-trips plus ProjectManager Kd state save/load.

import os
import sys
import types
import json
from pathlib import Path

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

REPO_ROOT = Path(__file__).resolve().parents[1]
_DYNAMIXS_DIR = REPO_ROOT / "modules" / "dynamiXs_v2o0"
if str(_DYNAMIXS_DIR) not in sys.path:
    sys.path.insert(0, str(_DYNAMIXS_DIR))

from PySide6.QtWidgets import QApplication  # noqa: E402

from lunaNMR.utils.project_manager import ProjectManager  # noqa: E402


@pytest.fixture(scope="module")
def app():
    return QApplication.instance() or QApplication(sys.argv)


def _fake_main_window():
    return types.SimpleNamespace(show_main_menu=lambda: None, current_dir="")


def test_methyl_t2_page_session_roundtrip(app):
    from methyl_t2_page import MethylT2FittingPage

    page = MethylT2FittingPage(_fake_main_window())
    page.input_file = "/data/methyl.csv"
    page.output_dir = "/data/out"
    page.last_json_folder = "/data/out/json"
    page.last_results_file = "/data/out/field1_methylT2_results.csv"
    page.t2a_spin.setValue(150.0)
    page.t2b_spin.setValue(25.0)
    page.field_freq_spin.setValue(800.0)
    page.boot_spin.setValue(2000)
    page.error_combo.setCurrentIndex(1)

    state = page.get_session_state()

    fresh = MethylT2FittingPage(_fake_main_window())
    fresh.restore_session_state(state)

    assert fresh.input_file == "/data/methyl.csv"
    assert fresh.output_dir == "/data/out"
    assert fresh.last_json_folder == "/data/out/json"
    assert fresh.last_results_file == "/data/out/field1_methylT2_results.csv"
    assert fresh.t2a_spin.value() == 150.0
    assert fresh.t2b_spin.value() == 25.0
    assert fresh.field_freq_spin.value() == 800.0
    assert fresh.boot_spin.value() == 2000
    assert fresh.error_combo.currentIndex() == 1


def test_kd_page_session_roundtrip(app):
    from kd_titration_page import KdTitrationPage

    page = KdTitrationPage(_fake_main_window())
    page.input_file = "/data/series_analysis_tidy.csv"
    page.output_dir = "/data/kd_out"
    page.detected_points = [0.0, 1.0, 2.0]
    page.last_json_file = "/data/kd_out/kd_fit_data.json"
    page.last_json_folder = "/data/kd_out"
    page._build_conc_rows(page.detected_points)
    for spin, val in zip(page.conc_spins, [0.0, 50.0, 100.0]):
        spin.setValue(val)
    for spin, val in zip(page.scale_spins, [1.0, 2.0, 0.5]):
        spin.setValue(val)
    page.p0_spin.setValue(75.0)
    page.alpha_spin.setValue(0.10)
    page.obs_combo.setCurrentIndex(1)
    page.intvalue_combo.setCurrentText("volume")
    page.boot_spin.setValue(500)

    state = page.get_session_state()

    fresh = KdTitrationPage(_fake_main_window())
    fresh.restore_session_state(state)

    assert fresh.input_file == "/data/series_analysis_tidy.csv"
    assert fresh.output_dir == "/data/kd_out"
    assert fresh.detected_points == [0.0, 1.0, 2.0]
    assert [s.value() for s in fresh.conc_spins] == [0.0, 50.0, 100.0]
    assert [s.value() for s in fresh.scale_spins] == [1.0, 2.0, 0.5]
    assert fresh.p0_spin.value() == 75.0
    assert fresh.alpha_spin.value() == 0.10
    assert fresh.obs_combo.currentIndex() == 1
    assert fresh.intvalue_combo.currentText() == "volume"
    assert fresh.boot_spin.value() == 500


def _kd_analyses():
    return {
        'A1_assi': {
            'fit_data': {'metadata': {'protein_conc': 75.0, 'alpha': 0.10,
                                      'points': [0.0, 50.0, 100.0]},
                         'fits': [{'residue': 'K14', 'series': {'ppm_x': [8.0, 8.1, 8.3]}}]},
            'meta': {'p0': 75.0, 'observable': 1, 'intensity_value': 'volume',
                     'input_basename': 'A1_assi.csv',
                     'comparison': {'ref': 0, 'cmp': 1}},
        },
        'B4_assi': {
            'fit_data': {'metadata': {'protein_conc': 50.0}, 'fits': []},
            'meta': {'p0': 50.0, 'input_basename': 'B4_assi.csv'},
        },
    }


def test_kd_multi_analysis_save_load_roundtrip(tmp_path):
    mw = types.SimpleNamespace(kd_titration_dialog=None, kd_analyses=_kd_analyses())
    pm = ProjectManager(mw)
    bundle = tmp_path / "proj.lunaNMR"
    bundle.mkdir()
    assert pm._save_kd_state(bundle) is True

    a1 = bundle / "kd" / "analyses" / "A1_assi"
    assert (a1 / "fit_data.json").exists() and (a1 / "meta.json").exists()
    fit = json.loads((a1 / "fit_data.json").read_text())
    assert fit['metadata']['protein_conc'] == 75.0
    assert fit['fits'][0]['series']['ppm_x'] == [8.0, 8.1, 8.3]
    meta = json.loads((a1 / "meta.json").read_text())
    assert meta['comparison'] == {'ref': 0, 'cmp': 1}

    mw2 = types.SimpleNamespace(kd_titration_dialog=None, kd_analyses=None)
    pm2 = ProjectManager(mw2)
    assert pm2._load_kd_state(bundle) is True
    assert set(mw2.kd_analyses) == {'A1_assi', 'B4_assi'}
    assert mw2.kd_analyses['A1_assi']['meta']['p0'] == 75.0
    # the bundled fit JSON path is exposed so the viewer can open it
    assert mw2.kd_analyses['A1_assi']['fit_data_path'].endswith('A1_assi/fit_data.json')


def test_kd_no_analyses_saves_nothing(tmp_path):
    mw = types.SimpleNamespace(kd_titration_dialog=None, kd_analyses={})
    pm = ProjectManager(mw)
    bundle = tmp_path / "proj.lunaNMR"
    bundle.mkdir()
    assert pm._save_kd_state(bundle) is False
    assert not (bundle / "kd" / "analyses").exists()


def test_kd_analyses_listed_in_inventory(tmp_path):
    mw = types.SimpleNamespace(kd_titration_dialog=None, kd_analyses=_kd_analyses())
    pm = ProjectManager(mw)
    bundle = tmp_path / "proj.lunaNMR"
    bundle.mkdir()
    pm._save_kd_state(bundle)
    cats = pm.inventory(bundle)
    kd = next(c for c in cats if c['id'] == 'kd')
    labels = {it['label'] for it in kd['items']}
    assert {'A1_assi', 'B4_assi'} <= labels
    assert all(it['removable'] for it in kd['items'])


def test_kd_dialog_state_roundtrip_and_close_transfer(app):
    from lunaNMR.gui.dialogs.kd_titration_dialog import KdTitrationDialog

    mw = types.SimpleNamespace(current_nmr_folder=None, kd_state=None, kd_file_refs=None)
    dlg = KdTitrationDialog(parent=None, main_window=mw)
    dlg.kd_page.input_file = "/d/x.csv"
    dlg.kd_page.detected_points = [0.0, 1.0]
    dlg.kd_page._build_conc_rows([0.0, 1.0])
    dlg.kd_page.p0_spin.setValue(60.0)

    assert dlg.get_state()['p0'] == 60.0

    dlg.close()  # closeEvent should transfer state to the main window
    assert mw.kd_state['p0'] == 60.0
    assert mw.kd_file_refs['input_file'] == "/d/x.csv"
    dlg.deleteLater()


def test_dynamixs_analyses_roundtrip_bundles_and_repoints_results(tmp_path):
    # A methyl-T2 run (results in its output_dir) is saved under
    # dynamixs/analyses/<name>/results/methyl_t2/, and load repoints into the bundle.
    out_dir = tmp_path / "methyl_out"
    (out_dir / "json").mkdir(parents=True)
    (out_dir / "field1_methylT2_results.csv").write_text("residue,t2a\nA,100\n")
    (out_dir / "json" / "A_methylT2_fit_data.json").write_text("{}")
    state = {'methyl_t2': {
        'output_dir': str(out_dir),
        'last_json_folder': str(out_dir / "json"),
        'last_results_file': str(out_dir / "field1_methylT2_results.csv"),
    }}
    mw = types.SimpleNamespace(dynamixs_dialog=None, dynamixs_analyses={
        'HSPA8_methylT2': {'state': state, 'file_refs': {}}})
    pm = ProjectManager(mw)
    bundle = tmp_path / "proj.lunaNMR"
    bundle.mkdir()
    assert pm._save_dynamixs_state(bundle) is True
    adir = bundle / "dynamixs" / "analyses" / "HSPA8_methylT2"
    assert (adir / "results" / "methyl_t2" / "field1_methylT2_results.csv").exists()
    assert (adir / "results" / "methyl_t2" / "json").is_dir()

    mw2 = types.SimpleNamespace(dynamixs_dialog=None, dynamixs_analyses=None)
    pm2 = ProjectManager(mw2)
    assert pm2._load_dynamixs_state(bundle) is True
    loaded = mw2.dynamixs_analyses['HSPA8_methylT2']['state']['methyl_t2']
    assert str(out_dir) not in loaded['last_json_folder']      # repointed into the bundle
    assert Path(loaded['last_json_folder']).is_dir()
    assert Path(loaded['output_dir']) == adir / "results" / "methyl_t2"


def test_dynamixs_save_prunes_removed_analyses(tmp_path):
    mw = types.SimpleNamespace(dynamixs_dialog=None, dynamixs_analyses={
        'keep': {'state': {}, 'file_refs': {}},
        'drop': {'state': {}, 'file_refs': {}}})
    pm = ProjectManager(mw)
    bundle = tmp_path / "proj.lunaNMR"
    bundle.mkdir()
    assert pm._save_dynamixs_state(bundle) is True
    adir = bundle / "dynamixs" / "analyses"
    assert {p.name for p in adir.iterdir()} == {'keep', 'drop'}
    del mw.dynamixs_analyses['drop']
    pm._save_dynamixs_state(bundle)
    assert {p.name for p in adir.iterdir()} == {'keep'}


def test_kd_fit_json_bundled_and_repointed(tmp_path):
    # A saved analysis carries its fit JSON; on load it must repoint into the
    # bundle (the original output folder may be gone).
    mw = types.SimpleNamespace(kd_titration_dialog=None, kd_analyses={
        'A1_assi': {'fit_data': {'metadata': {'protein_conc': 50.0}, 'fits': []},
                    'meta': {'p0': 50.0}}})
    pm = ProjectManager(mw)
    bundle = tmp_path / "proj.lunaNMR"
    bundle.mkdir()
    assert pm._save_kd_state(bundle) is True

    bundled = bundle / "kd" / "analyses" / "A1_assi" / "fit_data.json"
    assert bundled.exists(), "fit JSON should be bundled"

    mw2 = types.SimpleNamespace(kd_titration_dialog=None, kd_analyses=None)
    pm2 = ProjectManager(mw2)
    assert pm2._load_kd_state(bundle) is True
    assert mw2.kd_analyses['A1_assi']['fit_data_path'] == str(bundled)


def test_load_kd_state_absent_is_noop(tmp_path):
    bundle = tmp_path / "proj.lunaNMR"
    bundle.mkdir()
    mw = types.SimpleNamespace(kd_titration_dialog=None, kd_analyses=None)
    pm = ProjectManager(mw)
    assert pm._load_kd_state(bundle) is False
    assert mw.kd_analyses is None


def test_dynamixs_summary_reports_analysis_names(tmp_path):
    mw = types.SimpleNamespace(dynamixs_dialog=None, dynamixs_analyses={
        'HSPA1A_t1t2': {'state': {}, 'file_refs': {}}})
    pm = ProjectManager(mw)
    bundle = tmp_path / "proj.lunaNMR"
    bundle.mkdir()
    assert pm._save_dynamixs_state_with_summary(bundle) == {'analyses': ['HSPA1A_t1t2']}

    mw2 = types.SimpleNamespace(dynamixs_dialog=None, dynamixs_analyses=None)
    pm2 = ProjectManager(mw2)
    assert pm2._load_dynamixs_state_with_summary(bundle) == {'analyses': ['HSPA1A_t1t2']}


def test_page_store_analysis_autonames_and_stores(app, tmp_path):
    from kd_titration_page import KdTitrationPage
    page = KdTitrationPage(_fake_main_window())
    page.input_file = str(tmp_path / "A1_assi.csv")
    jf = tmp_path / "kd_fit_data.json"
    jf.write_text('{"metadata": {"protein_conc": 50.0}, "fits": []}')
    page.last_json_file = str(jf)

    name1 = page._store_current_analysis()
    assert name1 == "A1_assi"
    entry = page.main_window.kd_analyses["A1_assi"]
    assert entry["fit_data"]["metadata"]["protein_conc"] == 50.0
    assert entry["meta"]["input_basename"] == "A1_assi.csv"

    # saving the same input again gets a unique suffix
    assert page._store_current_analysis() == "A1_assi_2"
    assert set(page.main_window.kd_analyses) == {"A1_assi", "A1_assi_2"}


def test_page_store_analysis_none_without_fit(app, tmp_path):
    from kd_titration_page import KdTitrationPage
    page = KdTitrationPage(_fake_main_window())
    page.last_json_file = None
    assert page._store_current_analysis() is None


_FIT_JSON = {
    "metadata": {"protein_conc": 50.0, "alpha": 0.14,
                 "concentrations": [0.0, 10.0, 25.0], "points": [0, 10, 25],
                 "intensity_value": "height"},
    "fits": [{"residue": "K14",
              "series": {"ppm_x": [8.0, 8.1, 8.3], "ppm_y": [120.0, 120.5, 121.0],
                         "height": [1000.0, 500.0, 250.0], "volume": [2000.0, 1000.0, 500.0]},
              "csp": {"success": False}, "intensity": {"success": False}}],
    "global": {},
}


def test_kd_end_to_end_store_save_load_open(app, tmp_path):
    """Store a fit on the page → save project → load into a fresh window →
    open the saved analysis back into a new page (the new multi-analysis flow)."""
    from kd_titration_page import KdTitrationPage

    jf = tmp_path / "kd_fit_data.json"
    jf.write_text(json.dumps(_FIT_JSON))

    page = KdTitrationPage(_fake_main_window())
    page.input_file = str(tmp_path / "A1_assi.csv")
    page.last_json_file = str(jf)
    name = page._store_current_analysis()
    assert name == "A1_assi"

    pm = ProjectManager(page.main_window)
    bundle = tmp_path / "proj.lunaNMR"
    bundle.mkdir()
    assert pm._save_kd_state(bundle) is True

    mw2 = _fake_main_window()
    mw2.kd_analyses = None
    pm2 = ProjectManager(mw2)
    assert pm2._load_kd_state(bundle) is True

    page2 = KdTitrationPage(_fake_main_window())
    viewer = page2.open_saved_analysis(mw2.kd_analyses["A1_assi"])
    try:
        assert viewer is not None
        # the page now points at the bundled fit JSON, inside the project
        assert page2.last_json_file.endswith("A1_assi/fit_data.json")
        assert "proj.lunaNMR" in page2.last_json_file
    finally:
        viewer.deleteLater()


def test_kd_analysis_stored_on_real_main_window_through_dialog(app, tmp_path):
    """Regression: the page's main_window is the hosting dialog, but analyses must
    land on the REAL main window so ProjectManager (which holds the real one) saves
    them. Exercises the dialog indirection that a direct-page test bypasses."""
    from lunaNMR.gui.dialogs.kd_titration_dialog import KdTitrationDialog

    real_mw = types.SimpleNamespace(current_nmr_folder=None, kd_state=None,
                                    kd_file_refs=None, kd_titration_dialog=None,
                                    kd_analyses=None)
    dlg = KdTitrationDialog(parent=None, main_window=real_mw)
    jf = tmp_path / "kd_fit_data.json"
    jf.write_text(json.dumps(_FIT_JSON))
    dlg.kd_page.input_file = str(tmp_path / "B1_assi.csv")
    dlg.kd_page.last_json_file = str(jf)

    name = dlg.kd_page._store_current_analysis()
    assert name == "B1_assi"
    # MUST be on the real main window, not the dialog
    assert getattr(real_mw, 'kd_analyses', None) and "B1_assi" in real_mw.kd_analyses

    pm = ProjectManager(real_mw)
    bundle = tmp_path / "proj.lunaNMR"
    bundle.mkdir()
    assert pm._save_kd_state(bundle) is True
    assert (bundle / "kd" / "analyses" / "B1_assi" / "fit_data.json").exists()
    dlg.deleteLater()


def test_kd_load_skips_corrupt_analysis(tmp_path):
    """A corrupt analysis JSON must not abort loading the whole project; it is
    skipped and the rest load."""
    bundle = tmp_path / "proj.lunaNMR"
    good = bundle / "kd" / "analyses" / "good"
    bad = bundle / "kd" / "analyses" / "bad"
    good.mkdir(parents=True)
    bad.mkdir(parents=True)
    (good / "fit_data.json").write_text('{"metadata": {}, "fits": []}')
    (good / "meta.json").write_text('{"p0": 50.0}')
    (bad / "fit_data.json").write_text('{ this is not valid json')

    mw = types.SimpleNamespace(kd_titration_dialog=None, kd_analyses=None)
    pm = ProjectManager(mw)
    assert pm._load_kd_state(bundle) is True
    assert set(mw.kd_analyses) == {'good'}        # 'bad' skipped, not fatal


def test_kd_series_picker_lists_only_titration_runs(app):
    from kd_titration_page import KdTitrationPage
    mw = _fake_main_window()
    mw.saved_series = {
        'titr_A': types.SimpleNamespace(
            series_mode='titration',
            metadata={'csv_path': '/x/titr_A/series_analysis_tidy.csv'}),
        'relax_B': types.SimpleNamespace(
            series_mode='time', metadata={'csv_path': '/x/relax_B/tidy.csv'}),
        'legacy_C': types.SimpleNamespace(metadata={}),   # no series_mode -> 'time' -> excluded
    }
    mw.current_project_path = None
    page = KdTitrationPage(mw)
    got = page._get_available_series()
    assert [s['name'] for s in got] == ['titr_A']
    assert got[0]['csv_path'] == '/x/titr_A/series_analysis_tidy.csv'
    # the picker widget shows the run; the "no series" hint is hidden
    assert page.series_list_widget.count() == 1
    assert page.no_series_label.isVisible() is False


def test_kd_autosave_upserts_by_series_name(app, tmp_path):
    # Auto-save uses the source-series name and OVERWRITES (no _2 duplicate) so
    # re-saving the project keeps one fit per series.
    from kd_titration_page import KdTitrationPage
    page = KdTitrationPage(_fake_main_window())
    page.input_file = str(tmp_path / "A1_assi.csv")
    jf = tmp_path / "kd_fit_data.json"
    jf.write_text('{"metadata": {}, "fits": []}')
    page.last_json_file = str(jf)

    assert page.ensure_current_saved() == "A1_assi"
    assert page.ensure_current_saved() == "A1_assi"     # upsert, not A1_assi_2
    assert list(page.main_window.kd_analyses) == ["A1_assi"]


def test_kd_save_prunes_deleted_analyses(tmp_path):
    mw = types.SimpleNamespace(kd_titration_dialog=None, kd_analyses={
        'keep_me': {'fit_data': {}, 'meta': {}},
        'delete_me': {'fit_data': {}, 'meta': {}}})
    pm = ProjectManager(mw)
    bundle = tmp_path / "proj.lunaNMR"
    bundle.mkdir()
    assert pm._save_kd_state(bundle) is True
    adir = bundle / "kd" / "analyses"
    assert {p.name for p in adir.iterdir()} == {'keep_me', 'delete_me'}

    # user deletes one, re-saves → its folder is pruned
    del mw.kd_analyses['delete_me']
    assert pm._save_kd_state(bundle) is True
    assert {p.name for p in adir.iterdir()} == {'keep_me'}


def test_kd_page_lists_saved_fits_from_project(app):
    from kd_titration_page import KdTitrationPage
    mw = _fake_main_window()
    mw.kd_analyses = {'B4_assi': {'fit_data': {}, 'meta': {}},
                      'A1_assi': {'fit_data': {}, 'meta': {}}}
    page = KdTitrationPage(mw)                           # populates on setup
    shown = [page.saved_fits_list.item(i).text()
             for i in range(page.saved_fits_list.count())]
    assert shown == ['A1_assi', 'B4_assi']               # sorted
    assert page.no_saved_fits_label.isVisible() is False


def test_loading_input_autofills_output_dir(app, tmp_path):
    import pandas as pd
    from kd_titration_page import KdTitrationPage
    sub = tmp_path / "series_results_X"
    sub.mkdir()
    csv = sub / "series_analysis_tidy.csv"
    pd.DataFrame({"spectrum_name": ["0", "1"], "assignment": ["R1", "R1"],
                  "ppm_x": [8.0, 8.1], "ppm_y": [120.0, 120.0],
                  "height": [1000.0, 500.0], "volume": [2000.0, 1000.0]}).to_csv(csv, index=False)
    page = KdTitrationPage(_fake_main_window())
    page._set_input_file(str(csv))
    assert page.input_file == str(csv)
    # auto-filled to a 'kd_analysis' subfolder of the series (peak-integration) folder
    assert page.output_dir == str(sub / "kd_analysis")
    assert page.outdir_label.text() == str(sub / "kd_analysis")


def test_dynamixs_dialog_autosaves_active_run_by_series_and_type(app, tmp_path):
    from lunaNMR.gui.dialogs.dynamixs_dialog import DynamiXsDialog
    mw = types.SimpleNamespace(current_nmr_folder=None, dynamixs_state=None,
                               dynamixs_file_refs=None, dynamixs_dialog=None,
                               dynamixs_analyses=None)
    dlg = DynamiXsDialog(parent=None, main_window=mw)
    try:
        # On the methyl-T2 page with a completed run + a known source series.
        dlg.stack.setCurrentWidget(dlg.methyl_t2_page)
        dlg.set_source_series("HSPA8")
        dlg.methyl_t2_page.output_dir = str(tmp_path)
        dlg.methyl_t2_page.last_results_file = str(tmp_path / "field1_methylT2_results.csv")

        name = dlg.ensure_current_saved()
        assert name == "HSPA8_methylT2"                     # <series>_<type>
        assert "HSPA8_methylT2" in mw.dynamixs_analyses
        # re-capture upserts (no duplicate)
        assert dlg.ensure_current_saved() == "HSPA8_methylT2"
        assert list(mw.dynamixs_analyses) == ["HSPA8_methylT2"]
    finally:
        dlg.deleteLater()


def test_dynamixs_no_results_not_saved(app, tmp_path):
    from lunaNMR.gui.dialogs.dynamixs_dialog import DynamiXsDialog
    mw = types.SimpleNamespace(current_nmr_folder=None, dynamixs_state=None,
                               dynamixs_file_refs=None, dynamixs_dialog=None,
                               dynamixs_analyses=None)
    dlg = DynamiXsDialog(parent=None, main_window=mw)
    try:
        dlg.stack.setCurrentWidget(dlg.methyl_t2_page)       # no run done
        assert dlg.ensure_current_saved() is None
    finally:
        dlg.deleteLater()

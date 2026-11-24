# ABOUTME: A T1/T2 or methyl fit JSON holds no NaN token: a number with no value is null, so a strict
# ABOUTME: parser reads the file, and the fit viewers, a refit and `t1rho` read null as "no value".

import json
import os
import sys
from pathlib import Path

import numpy as np
import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

REPO_ROOT = Path(__file__).resolve().parents[1]
DATA = REPO_ROOT.parent / "data_test"
T1_700 = (DATA / "data_JA1_relaxation" / "run_20260825_1405" / "series" / "700_T1"
          / "peak_intensity_matrix.csv")
METHYL = (DATA / "methyl" / "sampleNBD" / "series_results_sampleNBD_20260818_124619"
          / "peak_intensity_matrix.csv")
# With the baseline fitted, the fit of this residue of the real 700 MHz T1 series does
# not converge: it is written with no value at all.
NO_FIT = "10ValH"
KEY = "field1_T1"

pytestmark = pytest.mark.skipif(not (T1_700.exists() and METHYL.exists()),
                                reason="data_test not present")

for folder in ("", "dynamiXs_T1_T2"):
    sys.path.insert(0, str(REPO_ROOT / "modules" / "dynamiXs_v2o0" / folder))


def _refuse(token):
    raise ValueError(f"{token} is not JSON")


def strict(text):
    """`text` parsed as a strict JSON parser does: NaN and Infinity are refused."""
    return json.loads(text, parse_constant=_refuse)


def _t1_run(module_name, folder, fit_baseline):
    import importlib
    import matplotlib
    matplotlib.use("Agg")
    folder.mkdir(parents=True)
    summary = importlib.import_module(module_name).run_analysis_with_params({
        'input_csv_file': str(T1_700), 'output_prefix': str(folder / 'P'),
        'results_txt_file': str(folder / 'P_fit_results.txt'), 'experiment_type': 'T1',
        'time_units': 'ms', 'error_method': 'analytical', 'field_name': 'P_field1',
        'field_freq': 700.0, 'json_folder': str(folder), 'fit_baseline': fit_baseline})
    return Path(summary['json_file'])


@pytest.fixture(scope="module")
def unfitted_run(tmp_path_factory):
    """A single-core T1 run with the baseline fitted: (fit JSON, results file)."""
    folder = tmp_path_factory.mktemp("t1") / "run"
    return _t1_run("fit_Tx_NMRRE", folder, True), folder / "P_fit_results.txt"


@pytest.mark.parametrize("module", ["fit_Tx_NMRRE", "fitmulti__Tx_NMRRE"])
@pytest.mark.parametrize("fit_baseline", [False, True], ids=["fixed", "free"])
def test_a_t1_fit_json_is_strict_json(module, fit_baseline, tmp_path):
    payload = strict(_t1_run(module, tmp_path / "run", fit_baseline).read_text())
    assert len(payload['fits']) == 60
    if fit_baseline:
        entry = next(f for f in payload['fits'] if f['residue'] == NO_FIT)
        assert entry['success'] is False
        assert entry['t2'] is None and entry['A'] is None and entry['window_ratio'] is None


def test_a_methyl_fit_json_is_strict_json(tmp_path):
    from lunaNMR.cli import main
    assert main(["dynamixs", "methyl-t2", "--input", str(METHYL), "--out", str(tmp_path),
                 "--prefix", "M"]) == 0
    text = (tmp_path / "M_field1_methylT2_fit_data.json").read_text()
    payload = strict(text)
    assert "null" in text, "this real run holds errors with no value"
    assert all(set(fit) == set(payload['fits'][0]) for fit in payload['fits'])


def test_reading_a_fit_json_gives_nan_for_a_number_with_no_value(unfitted_run):
    from fit_Tx_NMRRE import read_fit_json
    payload = read_fit_json(unfitted_run[0])
    entry = next(f for f in payload['fits'] if f['residue'] == NO_FIT)
    assert np.isnan(entry['t2']) and np.isnan(entry['A'])
    assert np.isnan(entry['fit_curve']['intensity']).all()
    assert entry['success'] is False
    assert payload['metadata']['results_file'] == "P_fit_results.txt"


@pytest.fixture(scope="module")
def app():
    from PySide6.QtWidgets import QApplication
    return QApplication.instance() or QApplication(sys.argv)


@pytest.fixture
def open_viewer(app, monkeypatch):
    import lunaNMR.gui.dialogs.dynamixs_dialog  # noqa: F401 - puts the module on sys.path
    import visualization.fit_viewer as fit_viewer
    for box in ("show_info", "show_warning", "show_error"):
        monkeypatch.setattr(fit_viewer, box, lambda *a, **k: None)
    opened = []

    def open_(json_file):
        opened.append(fit_viewer.FitViewer(json_files=[str(json_file)]))
        return opened[-1]
    yield open_
    for viewer in opened:
        viewer.close()
        viewer.deleteLater()


def _show(viewer, residue):
    from PySide6.QtCore import Qt
    for row in range(viewer.residue_list.count()):
        if viewer.residue_list.item(row).data(Qt.UserRole) == residue:
            viewer.residue_list.setCurrentRow(row)
            return viewer.residue_list.item(row), viewer.axes[0]
    raise KeyError(residue)


def _assert_shown_as_not_fitted(viewer):
    item, ax = _show(viewer, NO_FIT)
    assert "not fitted" in item.toolTip()
    assert len(ax.lines) == 0
    assert [t.get_text() for t in ax.texts] == ["T1 (700.0 MHz): not fitted"]
    assert sum(len(points.get_offsets()) for points in ax.collections) == 10


def test_the_fit_viewer_shows_a_residue_with_no_value_as_not_fitted(unfitted_run, open_viewer):
    _assert_shown_as_not_fitted(open_viewer(unfitted_run[0]))


def test_the_fit_viewer_opens_a_fit_json_that_holds_nan_tokens(unfitted_run, open_viewer, tmp_path):
    """The file as it was written before a number with no value became null."""
    from fit_Tx_NMRRE import read_fit_json
    with_tokens = tmp_path / unfitted_run[0].name
    with_tokens.write_text(json.dumps(read_fit_json(unfitted_run[0]), indent=2))
    assert "NaN" in with_tokens.read_text()
    _assert_shown_as_not_fitted(open_viewer(with_tokens))


def test_a_refit_of_a_residue_with_no_value_keeps_the_file_strict(unfitted_run, open_viewer,
                                                                 tmp_path):
    import shutil
    folder = tmp_path / "copy"
    shutil.copytree(unfitted_run[0].parent, folder)
    json_file = folder / unfitted_run[0].name
    viewer = open_viewer(json_file)
    _show(viewer, NO_FIT)
    time_points = viewer.data[KEY]["metadata"]["time_points"]
    viewer.exclusions.setdefault(KEY, {})[NO_FIT] = [time_points.index(max(time_points))]

    viewer._on_refit_clicked()

    payload = strict(json_file.read_text())
    entry = next(f for f in payload['fits'] if f['residue'] == NO_FIT)
    assert isinstance(entry['success'], bool)
    assert set(entry) == set(payload['fits'][0])


def test_the_methyl_viewer_reads_and_rewrites_strict_json(app, tmp_path, monkeypatch):
    from lunaNMR.cli import main
    import lunaNMR.gui.dialogs.dynamixs_dialog  # noqa: F401 - puts the module on sys.path
    import visualization.methyl_t2_fit_viewer as methyl_viewer
    for box in ("show_info", "show_warning", "show_error"):
        monkeypatch.setattr(methyl_viewer, box, lambda *a, **k: None)
    assert main(["dynamixs", "methyl-t2", "--input", str(METHYL), "--out", str(tmp_path),
                 "--prefix", "M"]) == 0
    json_file = tmp_path / "M_field1_methylT2_fit_data.json"
    without_value = [f['residue'] for f in strict(json_file.read_text())['fits']
                     if f['t2_a_err'] is None]
    assert without_value
    viewer = methyl_viewer.MethylT2FitViewer(json_file=str(json_file))
    try:
        (key,) = viewer.data
        entry = next(f for f in viewer.data[key]['fits'] if f['residue'] == without_value[0])
        assert np.isnan(entry['t2_a_err'])
        viewer._update_plot(without_value[0])
        viewer._refit_one(key, without_value[0])
    finally:
        viewer.close()
        viewer.deleteLater()
    strict(json_file.read_text())

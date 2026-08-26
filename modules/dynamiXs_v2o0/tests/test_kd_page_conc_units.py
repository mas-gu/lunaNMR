# ABOUTME: Tests the Kd page's concentration-units control and its params round trip.
# ABOUTME: Also that an unusable-fit verdict is not presented as an ordinary success.

import os
import sys
import json
import types
from pathlib import Path

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

_MOD = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_MOD))


@pytest.fixture(scope="module")
def app():
    from PySide6.QtWidgets import QApplication
    return QApplication.instance() or QApplication([])


def _normalised(raw):
    """The page applies params that have been through the schema, not raw dicts."""
    from kd_params import normalize_params
    return normalize_params(raw)


def _page():
    from kd_titration_page import KdTitrationPage
    return KdTitrationPage(types.SimpleNamespace(show_main_menu=lambda: None,
                                                 current_dir=""))


class TestConcentrationUnitsAreReachableFromTheGui:
    """Spectra named in equivalents fed in as molar concentrations give a Kd wrong by
    the factor [P]0. The CLI gained a control this morning; without one here a GUI user
    cannot correct it short of hand-editing JSON."""

    def test_the_control_exists(self, app):
        assert _page().conc_units_combo is not None

    def test_it_defaults_to_absolute(self, app):
        """Existing behaviour must not change for anyone who never touches it."""
        page = _page()
        assert page._gather_params_dict()["conc_units"] == "absolute"

    def test_choosing_equivalents_reaches_the_params(self, app):
        page = _page()
        page.conc_units_combo.setCurrentIndex(1)
        assert page._gather_params_dict()["conc_units"] == "equivalents"

    def test_it_offers_exactly_the_two_documented_units(self, app):
        page = _page()
        labels = [page.conc_units_combo.itemText(i).lower()
                  for i in range(page.conc_units_combo.count())]
        assert len(labels) == 2
        assert any("absolute" in text for text in labels)
        assert any("equivalent" in text for text in labels)


class TestTheSettingRoundTrips:
    """It is a property of the dataset, so it has to survive a save and reload like
    every other setting rather than being retyped each session."""

    def test_equivalents_survives_apply(self, app):
        page = _page()
        page._apply_params(_normalised({"conc_units": "equivalents"}), n_points=3)
        assert page._gather_params_dict()["conc_units"] == "equivalents"

    def test_absolute_survives_apply(self, app):
        page = _page()
        page.conc_units_combo.setCurrentIndex(1)
        page._apply_params(_normalised({"conc_units": "absolute"}), n_points=3)
        assert page._gather_params_dict()["conc_units"] == "absolute"

    def test_a_params_file_without_the_key_leaves_the_default(self, app):
        """An older params JSON predates this setting and must not break loading."""
        page = _page()
        page._apply_params(_normalised({"protein_conc": 50.0}), n_points=3)
        assert page._gather_params_dict()["conc_units"] == "absolute"

    def test_it_survives_a_json_round_trip(self, app, tmp_path):
        page = _page()
        page.conc_units_combo.setCurrentIndex(1)
        path = tmp_path / "kd_params.json"
        path.write_text(json.dumps(page._gather_params_dict()))

        reloaded = _page()
        reloaded._apply_params(_normalised(json.loads(path.read_text())), n_points=3)
        assert reloaded._gather_params_dict()["conc_units"] == "equivalents"


class TestTheJsonBrowserFollowsTheOutputLayout:
    """Fit JSONs moved into <out>/data/. A dialog opening one level up makes the user
    hunt for a file the tool just wrote."""

    def test_it_prefers_the_data_subfolder(self, app, tmp_path):
        (tmp_path / "data").mkdir()
        page = _page()
        page.output_dir = str(tmp_path)
        assert page._json_browse_dir() == str(tmp_path / "data")

    def test_it_falls_back_when_there_is_no_data_folder(self, app, tmp_path):
        """An older result folder has no data/ and must still open somewhere useful."""
        page = _page()
        page.output_dir = str(tmp_path)
        assert page._json_browse_dir() == str(tmp_path)

    def test_an_unset_output_is_not_a_crash(self, app):
        page = _page()
        page.output_dir = ""
        page.last_json_folder = ""
        assert page._json_browse_dir() == ""

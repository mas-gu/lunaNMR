# ABOUTME: Tests the Kd page's Save-analysis-to-project name-collision handling.
# ABOUTME: On a duplicate name the user is asked to replace, rename, or cancel.

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


def _page():
    from kd_titration_page import KdTitrationPage
    return KdTitrationPage(types.SimpleNamespace(show_main_menu=lambda: None, current_dir=""))


def _page_with_fit(tmp_path, monkeypatch):
    """A page with a stored fit JSON and a fixed source-series name, plus the modal
    info/warning popups silenced so a save doesn't block the test."""
    import kd_titration_page as mod
    monkeypatch.setattr(mod, "show_info", lambda *a, **k: None)
    monkeypatch.setattr(mod, "show_warning", lambda *a, **k: None)
    jf = tmp_path / "SER_kd_fit_data.json"
    jf.write_text(json.dumps({"metadata": {"name": "SER"}, "fits": [], "global": {}}))
    page = _page()
    page.last_json_file = str(jf)
    page._current_analysis_name = "SER"
    return page


def test_store_under_explicit_name(app, tmp_path, monkeypatch):
    page = _page_with_fit(tmp_path, monkeypatch)
    assert page._store_current_analysis(name="custom") == "custom"
    assert "custom" in page._lunanmr_main_window().kd_analyses


def test_first_save_uses_base_name_no_prompt(app, tmp_path, monkeypatch):
    page = _page_with_fit(tmp_path, monkeypatch)
    monkeypatch.setattr(page, "_prompt_name_collision",
                        lambda base: pytest.fail("should not prompt without a collision"))
    page._save_analysis_to_project()
    assert list(page._lunanmr_main_window().kd_analyses) == ["SER"]


def test_collision_replace_overwrites(app, tmp_path, monkeypatch):
    page = _page_with_fit(tmp_path, monkeypatch)
    page._save_analysis_to_project()                    # first save -> SER
    monkeypatch.setattr(page, "_prompt_name_collision", lambda base: base)   # Replace
    page._save_analysis_to_project()
    assert list(page._lunanmr_main_window().kd_analyses) == ["SER"]   # replaced, not suffixed


def test_collision_rename_adds_second(app, tmp_path, monkeypatch):
    page = _page_with_fit(tmp_path, monkeypatch)
    page._save_analysis_to_project()                    # SER
    monkeypatch.setattr(page, "_prompt_name_collision", lambda base: "SER_v2")
    page._save_analysis_to_project()
    assert set(page._lunanmr_main_window().kd_analyses) == {"SER", "SER_v2"}


def test_collision_cancel_leaves_unchanged(app, tmp_path, monkeypatch):
    page = _page_with_fit(tmp_path, monkeypatch)
    page._save_analysis_to_project()                    # SER
    stored = dict(page._lunanmr_main_window().kd_analyses)
    monkeypatch.setattr(page, "_prompt_name_collision", lambda base: None)   # Cancel
    page._save_analysis_to_project()
    assert page._lunanmr_main_window().kd_analyses == stored


def test_save_without_fit_warns_and_stores_nothing(app, tmp_path, monkeypatch):
    import kd_titration_page as mod
    warned = {"n": 0}
    monkeypatch.setattr(mod, "show_warning", lambda *a, **k: warned.__setitem__("n", warned["n"] + 1))
    monkeypatch.setattr(mod, "show_info", lambda *a, **k: None)
    page = _page()
    page.last_json_file = None
    page._save_analysis_to_project()
    assert warned["n"] == 1
    assert not (getattr(page._lunanmr_main_window(), "kd_analyses", None) or {})

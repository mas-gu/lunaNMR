# ABOUTME: Tests what a project save and load report about themselves, through the real ProjectManager:
# ABOUTME: the number of peaks in the bundled peak list, and the lunaNMR version the manifest records.

import json
import types
from pathlib import Path

import pytest

import lunaNMR

PEAK_LIST = (Path(__file__).resolve().parents[2] / "data_test" / "data_binding"
             / "DNAJA1_HSPA8" / "DNAJA1.txt")

pytestmark = pytest.mark.skipif(not PEAK_LIST.exists(), reason="data_test not present")


def _window(peak_list=None):
    return types.SimpleNamespace(integrator=types.SimpleNamespace(peak_list=peak_list))


@pytest.fixture
def peak_list():
    from lunaNMR.utils.file_manager import NMRFileManager
    return NMRFileManager().load_peak_list(str(PEAK_LIST))


@pytest.fixture
def saved(peak_list, tmp_path):
    """A project saved from a session holding the peak list: (bundle, save summary)."""
    from lunaNMR.utils.project_manager import ProjectManager
    bundle = tmp_path / "P.lunaNMR"
    ok, summary = ProjectManager(_window(peak_list)).save_project(bundle)
    assert ok
    return bundle, summary


def test_the_save_summary_counts_the_peaks_it_saved(saved, peak_list):
    _bundle, summary = saved
    assert len(peak_list) > 0
    assert f"Peak list ({len(peak_list)} peaks)" in summary['saved_items']


def test_the_load_summary_counts_the_peaks_it_loaded(saved, peak_list):
    from lunaNMR.utils.project_manager import ProjectManager
    bundle, _summary = saved
    window = _window()
    ok, _missing, summary = ProjectManager(window).load_project(bundle)
    assert ok
    assert len(window.integrator.peak_list) == len(peak_list)
    assert f"Peak list ({len(peak_list)} peaks)" in summary['loaded_items']


def test_a_session_without_peaks_reports_no_peak_list(tmp_path):
    from lunaNMR.utils.project_manager import ProjectManager
    ok, summary = ProjectManager(_window()).save_project(tmp_path / "P.lunaNMR")
    assert ok
    assert not [item for item in summary['saved_items'] if item.startswith("Peak list")]


def test_the_manifest_records_the_release_version(saved):
    bundle, _summary = saved
    manifest = json.loads((bundle / "project.json").read_text())
    assert manifest["lunaNMR_version"] == lunaNMR.__version__

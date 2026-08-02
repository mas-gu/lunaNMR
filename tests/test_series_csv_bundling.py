# ABOUTME: Tests that a series' titration CSV is bundled on save and resolved from the bundle.
# ABOUTME: Makes a saved project self-contained so a moved/renamed run folder still populates the Kd box.

import os
import sys
import types
from pathlib import Path

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

REPO_ROOT = Path(__file__).resolve().parents[1]
_DYNAMIXS_DIR = REPO_ROOT / "modules" / "dynamiXs_v2o0"
if str(_DYNAMIXS_DIR) not in sys.path:
    sys.path.insert(0, str(_DYNAMIXS_DIR))

from lunaNMR.processors.multi_spectrum_processor import BatchResults
from lunaNMR.utils.project_manager import ProjectManager


def _make_series(tmp_path):
    run = tmp_path / "run" / "series_results_S1"
    run.mkdir(parents=True)
    tidy = run / "series_analysis_tidy.csv"
    tidy.write_text("spectrum_name,ppm_x\nA,8.0\n")
    batch = BatchResults()
    batch.series_mode = 'titration'
    batch.metadata['csv_path'] = str(tidy)
    batch.metadata['output_folder'] = str(run)
    batch.add_result('A', {'status': 'success'})
    return batch, tidy


# ---------- save side: CSV copied into the bundle ----------

def test_save_copies_tidy_csv_into_bundle(tmp_path):
    batch, tidy = _make_series(tmp_path)
    mw = types.SimpleNamespace(saved_series={'S1': batch})
    bundle = tmp_path / "proj.lunaNMR"
    bundle.mkdir()
    ProjectManager(mw)._save_series_results(bundle)

    copied = bundle / "series_results" / "S1" / "series_analysis_tidy.csv"
    assert copied.exists()
    assert copied.read_text() == tidy.read_text()


# ---------- resolve side: bundled copy preferred, falls back correctly ----------

def _resolve():
    from kd_titration_page import resolve_series_csv
    return resolve_series_csv


def test_resolve_prefers_bundled_copy(tmp_path):
    resolve_series_csv = _resolve()
    bundle = tmp_path / "proj.lunaNMR"
    bundled = bundle / "series_results" / "S1" / "series_analysis_tidy.csv"
    bundled.parent.mkdir(parents=True)
    bundled.write_text("bundled\n")
    orig = tmp_path / "elsewhere" / "series_analysis_tidy.csv"
    orig.parent.mkdir(parents=True)
    orig.write_text("original\n")

    meta = {'csv_path': str(orig), 'output_folder': str(orig.parent)}
    got = resolve_series_csv('S1', meta, str(bundle))
    assert got == str(bundled)


def test_resolve_falls_back_to_recorded_path(tmp_path):
    resolve_series_csv = _resolve()
    bundle = tmp_path / "proj.lunaNMR"          # no bundled copy present
    orig = tmp_path / "elsewhere" / "series_analysis_tidy.csv"
    orig.parent.mkdir(parents=True)
    orig.write_text("original\n")

    meta = {'csv_path': str(orig), 'output_folder': str(orig.parent)}
    got = resolve_series_csv('S1', meta, str(bundle))
    assert got == str(orig)


def test_resolve_no_extra_lunaNMR_component(tmp_path):
    """The bundle dir IS the .lunaNMR path — resolution must not insert a second one."""
    resolve_series_csv = _resolve()
    bundle = tmp_path / "proj.lunaNMR"
    bundled = bundle / "series_results" / "S1" / "series_analysis_tidy.csv"
    bundled.parent.mkdir(parents=True)
    bundled.write_text("bundled\n")

    got = resolve_series_csv('S1', {}, str(bundle))
    assert got == str(bundled)
    assert ".lunaNMR/.lunaNMR" not in got.replace(os.sep, "/")

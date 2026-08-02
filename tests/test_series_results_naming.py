# ABOUTME: Tests that the series-run output folder embeds the user-provided series name.
# ABOUTME: series_results_<name>_<timestamp>, falling back to timestamp-only when unnamed.

import os
import sys
from pathlib import Path

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))


def test_folder_includes_series_name_and_timestamp():
    from lunaNMR.gui.dialogs.series_integration_dialog import _series_results_folder
    out = _series_results_folder("/data/spectra", "HSPA8", "20260722_153045")
    assert out == os.path.join("/data/spectra", "series_results_HSPA8_20260722_153045")


def test_folder_falls_back_to_timestamp_only_when_unnamed():
    from lunaNMR.gui.dialogs.series_integration_dialog import _series_results_folder
    for empty in (None, "", "   "):
        out = _series_results_folder("/data", empty, "20260722_153045")
        assert out == os.path.join("/data", "series_results_20260722_153045")

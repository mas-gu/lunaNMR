# ABOUTME: A titration series cannot run with the 'independent' peak source, in the CLI or GUI.
# ABOUTME: Titration tracks by cascade; independent re-anchors every point to the free state.

"""
Titration mode forces cascade tracking, but 'independent' also switches on
use_original_reference_for_detection, which titration left in place: every point was then
detected against the 0-eq peak list and the n-1 seed never used, while the log said
"cascade tracking". The combination is refused instead.
"""

import sys
from pathlib import Path

import pytest

DNAJA1 = Path(__file__).resolve().parents[2] / "data_test" / "data_binding" / "DNAJA1_HSPA8"


@pytest.mark.skipif(not DNAJA1.exists(), reason="data_test not present")
class TestTheCliRefusesIt:

    def _run(self, tmp_path, *extra):
        from lunaNMR.cli import main
        return main(["series", "--spectra", str(DNAJA1), "--peaks", str(DNAJA1 / "DNAJA1.txt"),
                     "--out", str(tmp_path / "o"), "--mode", "titration",
                     "--peak-source", "independent", *extra])

    def test_a_run_is_refused_and_says_titration_tracks_by_cascade(self, tmp_path, capsys):
        assert self._run(tmp_path) != 0
        err = capsys.readouterr().err
        assert "independent" in err and "cascade" in err
        assert not (tmp_path / "o").exists()

    def test_the_dry_run_refuses_it_too(self, tmp_path):
        assert self._run(tmp_path, "--dry-run") != 0

    def test_time_mode_still_accepts_independent(self, tmp_path):
        from lunaNMR.cli import main
        assert main(["series", "--spectra", str(DNAJA1), "--peaks", str(DNAJA1 / "DNAJA1.txt"),
                     "--out", str(tmp_path / "o"), "--mode", "time",
                     "--peak-source", "independent", "--dry-run"]) == 0


@pytest.fixture(scope="module")
def app():
    from PySide6.QtWidgets import QApplication
    return QApplication.instance() or QApplication(sys.argv)


@pytest.fixture
def dialog(app):
    from lunaNMR.gui.dialogs.series_integration_dialog import SeriesIntegrationDialog
    dlg = SeriesIntegrationDialog(main_window=None)
    yield dlg
    dlg.deleteLater()


TITRATION = 2


class TestTheDialogMakesItUnavailable:

    def test_independent_cannot_be_chosen_for_a_titration(self, dialog):
        dialog.xaxis_mode_combo.setCurrentIndex(TITRATION)
        assert not dialog.independent_radio.isEnabled()

    def test_choosing_titration_moves_an_independent_choice_to_cascade(self, dialog):
        dialog.independent_radio.setChecked(True)
        dialog.xaxis_mode_combo.setCurrentIndex(TITRATION)
        assert not dialog.independent_radio.isChecked()
        assert dialog.cascade_radio.isChecked()

    @pytest.mark.parametrize("index", [0, 1])
    def test_it_is_available_again_outside_titration(self, dialog, index):
        dialog.xaxis_mode_combo.setCurrentIndex(TITRATION)
        dialog.xaxis_mode_combo.setCurrentIndex(index)
        assert dialog.independent_radio.isEnabled()

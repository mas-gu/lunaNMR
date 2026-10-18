# ABOUTME: `diagnose` on real spectra whose peak list or folder holds something it cannot measure.
# ABOUTME: A peak outside the spectral window, or a file it cannot read, must not cost the rest.

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

REPO = Path(__file__).resolve().parents[2]
TITRATION = REPO / "data_test" / "data_binding" / "DNAJA1_HSPA8"
ONED = Path(__file__).resolve().parents[1] / "data_example" / "1D" / "1D_series_001.ft1"
SPECTRA = ("2D_DNAJA1_0o0.ft", "2D_DNAJA1_HSPA8_3o0.ft")

pytestmark = pytest.mark.skipif(not TITRATION.exists(), reason="data_test not present")


def _experiment(folder, extra_rows=()):
    """The titration's first and last spectra beside its peak list, plus `extra_rows`."""
    folder.mkdir(parents=True)
    for name in SPECTRA:
        (folder / name).symlink_to(TITRATION / name)
    text = (TITRATION / "DNAJA1.txt").read_text().rstrip("\n")
    rows = "".join(f"\n{name}, {x}, {y}" for name, x, y in extra_rows)
    (folder / "DNAJA1.txt").write_text(text + rows + "\n")
    return folder


class TestAPeakOutsideTheWindow:
    """The 15N window is about 100.5-132.5 ppm. Beyond its top a box was an empty slice
    and the whole run raised; below its bottom the index went negative and the box read
    a region nowhere near the peak, silently."""

    @pytest.fixture(scope="class")
    def clean(self, tmp_path_factory):
        from lunaNMR.validation.spectra_check import check_experiment
        d = _experiment(tmp_path_factory.mktemp("clean") / "exp")
        return check_experiment(d, d / "DNAJA1.txt", quick=True)

    @pytest.fixture(scope="class")
    def outside(self, tmp_path_factory):
        from lunaNMR.validation.spectra_check import check_experiment
        d = _experiment(tmp_path_factory.mktemp("outside") / "exp",
                        [("below", 8.3, 95.0), ("above", 8.3, 140.0)])
        return check_experiment(d, d / "DNAJA1.txt", quick=True)

    def test_every_spectrum_is_still_assessed(self, outside):
        assert [s["spectrum"] for s in outside["spectra"]] == list(SPECTRA)

    def test_the_peaks_outside_are_named(self, outside):
        messages = [f.message for f in outside["findings"] if f.check == "window"]
        assert len(messages) == 1
        assert "below" in messages[0] and "above" in messages[0]

    def test_they_do_not_move_what_is_measured_for_the_others(self, clean, outside):
        for c, o in zip(clean["spectra"], outside["spectra"]):
            assert o["capture"] == c["capture"]
            assert o["weak"] == c["weak"]
            assert o["median_snr"] == pytest.approx(c["median_snr"], rel=0.05)

    def test_a_box_outside_reads_nothing(self, tmp_path):
        from lunaNMR.validation.spectra_check import peak_heights
        heights, _ = peak_heights(TITRATION / SPECTRA[0],
                                  [("below", 8.3, 95.0), ("above", 8.3, 140.0),
                                   ("3LysH", 8.224692, 126.270531)])
        assert np.isnan(heights[0]) and np.isnan(heights[1])
        assert np.isfinite(heights[2])


@pytest.mark.skipif(not ONED.exists(), reason="data_example not present")
class TestAFolderThatCannotBeChecked:
    """One folder whose spectrum the 2D check cannot read used to abort the dataset, and
    every other experiment's result went with it."""

    @pytest.fixture(scope="class")
    def result(self, tmp_path_factory):
        from lunaNMR.validation.spectra_check import check_dataset
        root = tmp_path_factory.mktemp("root")
        _experiment(root / "good")
        (root / "oned").mkdir()
        (root / "oned" / ONED.name).symlink_to(ONED)
        (root / "oned" / "DNAJA1.txt").write_text((TITRATION / "DNAJA1.txt").read_text())
        return check_dataset(root, quick=True)

    def test_the_readable_experiment_is_still_reported(self, result):
        assert [e["experiment"] for e in result["experiments"]] == ["good"]

    def test_the_other_is_a_fail_naming_it(self, result):
        fails = [f for f in result["findings"] if f.severity == "FAIL"]
        assert len(fails) == 1
        assert fails[0].check == "spectra" and "oned" in fails[0].message

# ABOUTME: A series detects its first spectrum with the search window it was given.
# ABOUTME: Real DNAJA1 spectrum, peak list shifted 0.03 ppm in 1H: inside 0.070, outside the 0.01 default.

"""
The series parameters reached the integrator in process_peak_list, which runs after
detection, so spectrum 1 was detected with the integrator's own defaults (0.01 / 0.05 ppm)
whatever window was set - and cascade and titration never detect again, so there the
setting had no effect at all.
"""

import json
from pathlib import Path

import pandas as pd
import pytest

DNAJA1 = Path(__file__).resolve().parents[2] / "data_test" / "data_binding" / "DNAJA1_HSPA8"
SHIFT = 0.03

pytestmark = pytest.mark.skipif(not DNAJA1.exists(), reason="data_test not present")


def test_spectrum_one_is_detected_with_the_series_window(tmp_path):
    from lunaNMR.cli import main
    from lunaNMR.utils.file_manager import NMRFileManager
    peaks = NMRFileManager().load_peak_list(str(DNAJA1 / "DNAJA1.txt"))
    peaks["Position_X"] = peaks["Position_X"] + SHIFT
    shifted = tmp_path / "shifted.txt"
    peaks[["Assignment", "Position_X", "Position_Y"]].to_csv(shifted, sep="\t", index=False)
    (tmp_path / "spectra").mkdir()
    (tmp_path / "spectra" / "2D_DNAJA1_0o0.ft").symlink_to(DNAJA1 / "2D_DNAJA1_0o0.ft")
    assert main(["series", "--spectra", str(tmp_path / "spectra"), "--peaks", str(shifted),
                 "--mode", "time", "--peak-source", "reference",
                 "--out", str(tmp_path / "out")]) == 0
    batch = json.loads((tmp_path / "out" / "batch_results.json").read_text())
    fitted = pd.DataFrame(next(iter(batch["results"].values()))["fitted_results"])
    matched = fitted[~fitted["reference_retained"].astype(bool)]
    assert len(matched) > len(fitted) / 2
    assert matched["distance_from_reference_x"].max() > 0.01

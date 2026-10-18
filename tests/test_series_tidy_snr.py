# ABOUTME: The tidy series CSV reports each peak's S/N, negative peaks included.
# ABOUTME: Real DNAJA1 reference spectrum run through `lunaNMR series`.

"""
The `snr` column was 0 on every row: the converter read snr / peak_snr / noise_level off
the fit, which reports none of them, and its last fallback skipped height <= 0. The fit
reports detection_sigma, |height| / noise, which is the same quantity.
"""

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

DNAJA1 = Path(__file__).resolve().parents[2] / "data_test" / "data_binding" / "DNAJA1_HSPA8"

pytestmark = pytest.mark.skipif(not DNAJA1.exists(), reason="data_test not present")


@pytest.fixture(scope="module")
def tidy(tmp_path_factory):
    from lunaNMR.cli import main
    root = tmp_path_factory.mktemp("snr")
    (root / "spectra").mkdir()
    (root / "spectra" / "2D_DNAJA1_0o0.ft").symlink_to(DNAJA1 / "2D_DNAJA1_0o0.ft")
    assert main(["series", "--spectra", str(root / "spectra"), "--peaks", str(DNAJA1 / "DNAJA1.txt"),
                 "--mode", "time", "--peak-source", "reference", "--out", str(root / "out")]) == 0
    return pd.read_csv(root / "out" / "series_analysis_tidy.csv")


def test_every_fitted_peak_reports_its_snr(tidy):
    reported = tidy.dropna(subset=["detection_sigma"])
    assert len(reported) == len(tidy)
    np.testing.assert_allclose(reported.snr, reported.detection_sigma, rtol=1e-4)


def test_a_negative_peak_has_a_positive_snr(tidy):
    negative = tidy[tidy.height < 0]
    assert len(negative)
    assert (negative.snr > 0).all()

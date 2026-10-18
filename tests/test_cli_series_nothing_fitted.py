# ABOUTME: A series that fits nothing is reported as NothingFitted, as the other fit commands are.
# ABOUTME: Real JA1 T1 spectrum, truncated, so it cannot be read.

import contextlib
import io
import json
from pathlib import Path

import pytest

T1 = (Path(__file__).resolve().parents[2] / "data_test" / "data_JA1_relaxation"
      / "600_WT" / "T1")

pytestmark = pytest.mark.skipif(not T1.exists(), reason="data_test not present")


def test_a_series_that_fits_nothing_says_so(tmp_path):
    # Imported before pytest.warns, which re-emits every warning it does not match -
    # nmrglue's own import-time deprecations among them
    import nmrglue  # noqa: F401
    from lunaNMR.cli import main
    (tmp_path / "spectra").mkdir()
    (tmp_path / "spectra" / "T1_150ms.ft").write_bytes(
        next(T1.glob("*_150ms.ft")).read_bytes()[:4096])
    out = io.StringIO()
    # nmrglue warns on the truncated file - the condition under test
    with pytest.warns(UserWarning, match="cannot be shaped"), \
            contextlib.redirect_stdout(out), contextlib.redirect_stderr(io.StringIO()):
        code = main(["series", "--spectra", str(tmp_path / "spectra"),
                     "--peaks", str(next(T1.glob("*_assignment.txt"))), "--mode", "time",
                     "--peak-source", "reference", "--out", str(tmp_path / "o"),
                     "--format", "json"])
    envelope = json.loads(out.getvalue())
    assert code == 1
    assert envelope["ok"] is False
    assert envelope["error"]["type"] == "NothingFitted"

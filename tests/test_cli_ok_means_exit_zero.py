# ABOUTME: Under --format json, `ok` is true exactly when the exit code is 0.
# ABOUTME: Runs that emitted their summary and then returned 1 still said ok: true.

import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
DATA = Path(__file__).resolve().parents[2] / "data_test"
RUN = DATA / "data_JA1_relaxation" / "run_20260825_1405" / "series"
JA1_PEAKS = DATA / "data_JA1_relaxation" / "600_WT" / "T1" / "600_T1_assignment.txt"
METHYL = (DATA / "methyl" / "sampleNBD" / "series_results_sampleNBD_20260818_124619"
          / "peak_intensity_matrix.csv")
DNAJA1 = DATA / "data_binding" / "DNAJA1_HSPA8"

pytestmark = pytest.mark.skipif(not DATA.exists(), reason="data_test not present")

_ID_COLUMNS = ('Peak_Number', 'Assignment', 'Reference_X', 'Reference_Y')


def _run(*argv):
    return subprocess.run([sys.executable, "-m", "lunaNMR", *argv, "--format", "json"],
                          cwd=str(REPO_ROOT), capture_output=True, text=True)


def _failure(p, error_type):
    """The payload of a run that exited 1, checked to say so."""
    assert p.returncode == 1, p.stderr[-800:]
    payload = json.loads(p.stdout)
    assert payload["ok"] is False
    assert payload["error"]["type"] == error_type
    assert payload["error"]["message"]
    return payload


def _delay_columns(df):
    return [c for c in df.columns if c not in _ID_COLUMNS]


class TestADryRunWithAMissingInput:

    def test_kd(self, tmp_path):
        missing = str(tmp_path / "nope.csv")
        payload = _failure(_run("kd", "--input", missing, "--out", str(tmp_path / "o"),
                                "--p0", "50", "--dry-run"), "MissingInput")
        assert payload["missing_inputs"] == [missing]
        assert missing in payload["error"]["message"]

    def test_series(self, tmp_path):
        missing = str(tmp_path / "nope.txt")
        payload = _failure(_run("series", "--spectra", str(DNAJA1), "--peaks", missing,
                                "--out", str(tmp_path / "o"), "--dry-run"), "MissingInput")
        assert payload["missing_inputs"] == [missing]
        assert payload["dry_run"] is True


class TestNothingFitted:

    def test_t1t2(self, tmp_path):
        """A real T2 matrix held at its first delay: no decay, so every fit is excluded."""
        df = pd.read_csv(RUN / "600_T2" / "peak_intensity_matrix.csv")
        delays = _delay_columns(df)
        for c in delays:
            df[c] = df[delays[0]]
        csv = tmp_path / "flat.csv"
        df.to_csv(csv, index=False)
        payload = _failure(_run("dynamixs", "t1t2", "--exp", "T2", "--input", str(csv),
                                "--out", str(tmp_path / "o")), "NothingFitted")
        assert payload["n_fitted"] == 0 and payload["n_total"] == len(df)

    def test_methyl_t2(self, tmp_path):
        """A real methyl matrix with only its first delay measured: no fit can run."""
        df = pd.read_csv(METHYL)
        for c in _delay_columns(df)[1:]:
            df[c] = np.nan
        csv = tmp_path / "one_point.csv"
        df.to_csv(csv, index=False)
        payload = _failure(_run("dynamixs", "methyl-t2", "--input", str(csv),
                                "--out", str(tmp_path / "o")), "NothingFitted")
        assert payload["n_fitted"] == 0 and payload["n_total"] == len(df)

    def test_t1rho(self, tmp_path):
        """The T1 series given as the T1rho one and the T2 series as T1: R1rho falls
        below R1 cos^2(theta) for every residue, so no R2 comes out positive."""
        payload = _failure(_run(
            "dynamixs", "t1rho", "--input", str(RUN / "600_T1" / "peak_intensity_matrix.csv"),
            "--t1", str(RUN / "600_T2" / "peak_intensity_matrix.csv"),
            "--peaks", str(JA1_PEAKS), "--omega1", "2000", "--carrier", "118",
            "--theta", "30", "--field-freq", "600", "--out", str(tmp_path / "o")),
            "NothingFitted")
        assert payload["n_residues"] == 0


def _misregistered(root):
    """A real spectrum beside its own peak list moved 0.06 ppm in 1H."""
    folder = root / "titration"
    folder.mkdir()
    os.symlink(DNAJA1 / "2D_DNAJA1_0o0.ft", folder / "2D_DNAJA1_0o0.ft")
    peaks = pd.read_csv(DNAJA1 / "DNAJA1.txt", skipinitialspace=True)
    peaks["Position_X"] = peaks["Position_X"] + 0.06
    peaks.to_csv(folder / "DNAJA1.txt", index=False)
    return folder


class TestAPreflightThatFails:

    def test_diagnose(self, tmp_path):
        _misregistered(tmp_path)
        payload = _failure(_run("diagnose", str(tmp_path), "--quick"), "PreflightFailed")
        assert any(f["severity"] == "FAIL" for f in payload["findings"])

    def test_diagnose_strict_on_a_warning(self, tmp_path):
        """A folder of spectra with no peak list is a WARN, which --strict gates on."""
        folder = tmp_path / "titration"
        folder.mkdir()
        os.symlink(DNAJA1 / "2D_DNAJA1_0o0.ft", folder / "2D_DNAJA1_0o0.ft")
        assert json.loads(_run("diagnose", str(tmp_path), "--quick").stdout)["ok"] is True
        payload = _failure(_run("diagnose", str(tmp_path), "--quick", "--strict"),
                           "PreflightFailed")
        assert {f["severity"] for f in payload["findings"]} == {"WARN"}

    def test_series_deep_dry_run(self, tmp_path):
        folder = _misregistered(tmp_path)
        payload = _failure(_run("series", "--spectra", str(folder),
                                "--peaks", str(folder / "DNAJA1.txt"),
                                "--out", str(tmp_path / "o"), "--dry-run", "--deep"),
                           "PreflightFailed")
        assert payload["missing_inputs"] == []
        assert any(f["severity"] == "FAIL" for f in payload["checks"]["findings"])

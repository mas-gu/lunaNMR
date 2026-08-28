# ABOUTME: The relaxation subcommands must read series_metadata.json beside their input.
# ABOUTME: It records the delay units and which spectra had none - the CLI asked the user instead.

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))


def _series_output(tmp_path, headers, spectra, value_units="ms", n_unparsed=0):
    """A series output folder: an intensity matrix plus the sidecar `series` writes."""
    out = tmp_path / "series_out"
    out.mkdir(exist_ok=True)
    decay = 1000.0 * np.exp(-np.array([8.0, 51.0, 102.0, 204.0]) / 120.0) + 5.0
    rows = [["Peak_Number", "Assignment", "Reference_X", "Reference_Y"] + list(headers)]
    for i, name in enumerate(("K3", "E4", "V5"), start=1):
        rows.append([i, name, 8.0, 120.0] + [float(v) * (1 + 0.01 * i) for v in decay])
    pd.DataFrame(rows).to_csv(out / "peak_intensity_matrix.csv", header=False, index=False)
    (out / "series_metadata.json").write_text(json.dumps({
        "series_mode": "time",
        "value_units": value_units,
        "n_spectra": len(headers),
        "n_value_unparsed": n_unparsed,
        "columns": [{"column": h, "spectrum": s, "value": None if i >= len(headers) - n_unparsed else 1.0}
                    for i, (h, s) in enumerate(zip(headers, spectra))],
        "repeat_scale": None,
    }))
    return out / "peak_intensity_matrix.csv"


def _run(argv):
    from lunaNMR.cli import main
    return main(argv)


class TestUnparsedDelaysAreRefused:
    """`series` records how many spectra had no parseable delay. Those columns cannot be
    fitted, so the run quietly measures fewer points than the user acquired — the sidecar
    is the only thing that knows, and nothing read it.
    """

    def test_a_run_with_unparsed_spectra_is_refused(self, tmp_path):
        csv = _series_output(
            tmp_path,
            headers=["8", "51", "102", "ref_noCa_004"],
            spectra=["a_8ms.ft", "a_51ms.ft", "a_102ms.ft", "a_ref_noCa_004.ft"],
            n_unparsed=1)
        assert _run(["dynamixs", "t1t2", "--exp", "T2", "--input", str(csv),
                     "--out", str(tmp_path / "o")]) == 1

    def test_the_refusal_names_the_offending_spectrum(self, tmp_path, capsys):
        csv = _series_output(
            tmp_path,
            headers=["8", "51", "102", "ref_noCa_004"],
            spectra=["a_8ms.ft", "a_51ms.ft", "a_102ms.ft", "a_ref_noCa_004.ft"],
            n_unparsed=1)
        _run(["dynamixs", "t1t2", "--exp", "T2", "--input", str(csv),
              "--out", str(tmp_path / "o")])
        assert "a_ref_noCa_004.ft" in capsys.readouterr().err

    def test_it_can_be_overridden_deliberately(self, tmp_path):
        """A hetNOE folder legitimately has no delays, and someone may want the
        surviving points anyway. The escape hatch must be explicit, not the default."""
        csv = _series_output(
            tmp_path,
            headers=["8", "51", "102", "ref_noCa_004"],
            spectra=["a_8ms.ft", "a_51ms.ft", "a_102ms.ft", "a_ref_noCa_004.ft"],
            n_unparsed=1)
        assert _run(["dynamixs", "t1t2", "--exp", "T2", "--input", str(csv),
                     "--out", str(tmp_path / "o"), "--allow-unparsed-delays"]) == 0

    def test_a_clean_series_runs(self, tmp_path):
        csv = _series_output(
            tmp_path,
            headers=["8", "51", "102", "204"],
            spectra=["a_8ms.ft", "a_51ms.ft", "a_102ms.ft", "a_204ms.ft"])
        assert _run(["dynamixs", "t1t2", "--exp", "T2", "--input", str(csv),
                     "--out", str(tmp_path / "o")]) == 0

    def test_no_sidecar_is_not_an_error(self, tmp_path):
        """Hand-built and DynamiXs-format matrices have no sidecar; they must still run."""
        csv = _series_output(
            tmp_path,
            headers=["8", "51", "102", "204"],
            spectra=["a_8ms.ft", "a_51ms.ft", "a_102ms.ft", "a_204ms.ft"])
        (csv.parent / "series_metadata.json").unlink()
        assert _run(["dynamixs", "t1t2", "--exp", "T2", "--input", str(csv),
                     "--out", str(tmp_path / "o")]) == 0


class TestTheUnitsComeFromTheSidecar:
    """`series` normalises every delay to ms and records that. Asking the user to
    re-assert it by hand is what produces the documented 1000x error on R1.
    """

    def _summary(self, tmp_path, extra, value_units="ms"):
        csv = _series_output(
            tmp_path,
            headers=["8", "51", "102", "204"],
            spectra=["a_8ms.ft", "a_51ms.ft", "a_102ms.ft", "a_204ms.ft"],
            value_units=value_units)
        out = tmp_path / "o"
        assert _run(["dynamixs", "t1t2", "--exp", "T2", "--input", str(csv),
                     "--out", str(out), "--format", "json", *extra]) == 0
        return csv, out

    def test_the_sidecar_supplies_the_units(self, tmp_path, capsys):
        """t1t2 defaults --time-units to 's', but a series-produced matrix is in ms."""
        self._summary(tmp_path, [])
        summary = json.loads(capsys.readouterr().out)
        assert summary["time_units"] == "ms"
        assert summary["series_metadata"]["value_units"] == "ms"

    def test_an_explicit_flag_that_disagrees_is_refused(self, tmp_path):
        """One of the two is wrong and guessing is the 1000x error. Refuse instead."""
        csv = _series_output(
            tmp_path,
            headers=["8", "51", "102", "204"],
            spectra=["a_8ms.ft", "a_51ms.ft", "a_102ms.ft", "a_204ms.ft"],
            value_units="ms")
        assert _run(["dynamixs", "t1t2", "--exp", "T2", "--input", str(csv),
                     "--out", str(tmp_path / "o"), "--time-units", "s"]) == 1

    def test_an_explicit_flag_that_agrees_is_fine(self, tmp_path):
        self._summary(tmp_path, ["--time-units", "ms"])

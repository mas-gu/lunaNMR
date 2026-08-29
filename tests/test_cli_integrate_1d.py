# ABOUTME: `integrate-1d` measures picked peaks across a 1D series, headlessly.
# ABOUTME: The region sum is the quantitative observable; the fitted area is diagnostic only.

import json
import sys
from pathlib import Path

import numpy as np
import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))
sys.path.insert(0, str(REPO_ROOT / "modules" / "integration_1d_v1o0"))
ng = pytest.importorskip("nmrglue")


def _write_series(tmp_path, n=4):
    """A 1:1 conversion: substrate at 5.0 ppm decays, product at 5.3 ppm grows."""
    from oned_voigt import voigt_1d
    ppm = np.linspace(6.0, 4.0, 4000)
    paths = []
    for k in range(n):
        f = k / (n - 1)
        data = (voigt_1d(ppm, 1.0 - 0.8 * f, 5.00, 0.0015, 0.0005)
                + voigt_1d(ppm, 0.2 + 0.8 * f, 5.30, 0.0015, 0.0005))
        data = data + np.random.default_rng(k).normal(0, 1e-4, data.shape)
        path = tmp_path / f"series_{k:03d}.ft1"
        _write_pipe_1d(path, data, ppm)
        paths.append(path)
    return paths


def _write_pipe_1d(path, data, ppm):
    udic = ng.fileiobase.create_blank_udic(1)
    sw = float(abs(ppm[0] - ppm[-1]) * 600.0)
    udic[0].update(dict(sw=sw, obs=600.0, car=float(np.mean(ppm)) * 600.0,
                        size=len(data), label='1H', complex=False, time=False, freq=True))
    dic = ng.pipe.create_dic(udic)
    ng.pipe.write(str(path), dic, np.asarray(data, 'float32'), overwrite=True)


def _run(argv):
    from lunaNMR.cli import main
    return main(argv)


class TestItMeasuresASeries:
    def test_a_csv_is_written_with_a_row_per_spectrum(self, tmp_path, capsys):
        paths = _write_series(tmp_path)
        out = tmp_path / "out"
        assert _run(["integrate-1d", "--spectra", str(tmp_path), "--out", str(out),
                     "--peaks", "5.00,5.30", "--format", "json"]) == 0
        payload = json.loads(capsys.readouterr().out)
        assert payload["n_spectra"] == len(paths)
        assert payload["n_peaks"] == 2
        assert Path(payload["outputs"]["series"]).exists()

    def test_the_long_table_carries_every_observable(self, tmp_path, capsys):
        import csv
        _write_series(tmp_path)
        out = tmp_path / "out"
        _run(["integrate-1d", "--spectra", str(tmp_path), "--out", str(out),
              "--peaks", "5.00,5.30", "--format", "json"])
        table = Path(json.loads(capsys.readouterr().out)["outputs"]["table"])
        cols = next(csv.reader(table.open()))
        for col in ("spectrum", "assignment", "ppm", "height", "area", "fit_area",
                    "r_squared", "fwhm"):
            assert col in cols

    def test_the_conversion_is_tracked(self, tmp_path, capsys):
        """Substrate falls, product rises. If the two peaks collapsed onto one maximum
        both would move together — the failure competitive matching exists to stop."""
        import csv
        _write_series(tmp_path)
        out = tmp_path / "out"
        _run(["integrate-1d", "--spectra", str(tmp_path), "--out", str(out),
              "--peaks", "5.00,5.30", "--format", "json"])
        rows = list(csv.DictReader(
            Path(json.loads(capsys.readouterr().out)["outputs"]["table"]).open()))
        by_peak = {}
        for r in rows:
            by_peak.setdefault(r["assignment"], []).append(float(r["area"]))
        first, second = by_peak.values()
        assert (first[-1] - first[0]) * (second[-1] - second[0]) < 0, by_peak


class TestTheObservableWarningReachesTheUser:
    """On the reference series a 1:1 conversion conserves to 6.4% for the region sum,
    11.5% for height and 38% for the FITTED area — because a close neighbour caps the
    window at ~1.4 linewidths and leaves the sigma/gamma split unconstrained, while R2
    stays above 0.97 throughout and does not validate the extrapolation. In the GUI you
    would see that; a CLI caller has no plot, so the output has to say it.
    """

    def test_the_summary_names_the_quantitative_observable(self, tmp_path, capsys):
        _write_series(tmp_path)
        _run(["integrate-1d", "--spectra", str(tmp_path), "--out", str(tmp_path / "o"),
              "--peaks", "5.00,5.30", "--format", "json"])
        payload = json.loads(capsys.readouterr().out)
        assert payload["quantitative_observable"] == "area"
        assert "fit_area" in payload["diagnostic_observables"]

    def test_the_human_output_says_which_number_to_use(self, tmp_path, capsys):
        _write_series(tmp_path)
        _run(["integrate-1d", "--spectra", str(tmp_path), "--out", str(tmp_path / "o"),
              "--peaks", "5.00,5.30"])
        out = capsys.readouterr().out
        assert "area" in out and "fit_area" in out


class TestInputHandling:
    def test_peaks_can_be_detected_instead_of_given(self, tmp_path, capsys):
        _write_series(tmp_path)
        assert _run(["integrate-1d", "--spectra", str(tmp_path),
                     "--out", str(tmp_path / "o"), "--detect", "--format", "json"]) == 0
        assert json.loads(capsys.readouterr().out)["n_peaks"] >= 2

    def test_it_needs_peaks_from_somewhere(self, tmp_path):
        _write_series(tmp_path)
        with pytest.raises(SystemExit) as exc:
            _run(["integrate-1d", "--spectra", str(tmp_path), "--out", str(tmp_path / "o")])
        assert exc.value.code == 2

    def test_no_spectra_is_an_error(self, tmp_path):
        assert _run(["integrate-1d", "--spectra", str(tmp_path / "empty"),
                     "--out", str(tmp_path / "o"), "--peaks", "5.0"]) == 1

    def test_dry_run_writes_nothing(self, tmp_path):
        _write_series(tmp_path)
        out = tmp_path / "o"
        assert _run(["integrate-1d", "--spectra", str(tmp_path), "--out", str(out),
                     "--peaks", "5.00", "--dry-run"]) == 0
        assert not out.exists()

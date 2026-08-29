# ABOUTME: `peaks shift` corrects a rigidly offset peak list, which diagnose only measured.
# ABOUTME: A mis-registered list is the most common recoverable error in the series workflow.

import json
import subprocess
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
ng = pytest.importorskip("nmrglue")

from test_spectra_check import _peaks, _spectrum, _write_list

ROOT = Path(__file__).resolve().parent.parent


def _run(*argv):
    return subprocess.run([sys.executable, "-m", "lunaNMR", *argv],
                          cwd=str(ROOT), capture_output=True, text=True)


def _dataset(tmp_path, shift=(0.06, 0.6)):
    """A spectrum whose peaks sit `shift` away from where the list says they are."""
    peaks = _peaks()
    spec = tmp_path / "T1_50ms.ft"
    _spectrum(spec, peaks, shift=shift)
    return spec, _write_list(tmp_path / "list.txt", peaks)


def _positions(path):
    from lunaNMR.utils.file_manager import NMRFileManager
    df = NMRFileManager().load_peak_list(str(path))
    return list(zip(df['Assignment'].astype(str),
                    df['Position_X'].astype(float), df['Position_Y'].astype(float)))


class TestAnExplicitShift:
    def test_every_peak_moves_by_the_given_offset(self, tmp_path):
        _, peaks = _dataset(tmp_path)
        out = tmp_path / "shifted.txt"
        assert _run("peaks", "shift", "--peaks", str(peaks), "--out", str(out),
                    "--dx", "0.05", "--dy", "-0.5").returncode == 0
        for (n0, x0, y0), (n1, x1, y1) in zip(_positions(peaks), _positions(out)):
            assert n0 == n1
            assert x1 == pytest.approx(x0 + 0.05, abs=1e-6)
            assert y1 == pytest.approx(y0 - 0.5, abs=1e-6)

    def test_assignments_and_order_survive(self, tmp_path):
        _, peaks = _dataset(tmp_path)
        out = tmp_path / "shifted.txt"
        _run("peaks", "shift", "--peaks", str(peaks), "--out", str(out), "--dx", "0.01")
        assert [n for n, _, _ in _positions(out)] == [n for n, _, _ in _positions(peaks)]

    def test_it_refuses_to_overwrite_its_input(self, tmp_path):
        """The original is the only copy of an assignment nobody wants to redo."""
        _, peaks = _dataset(tmp_path)
        p = _run("peaks", "shift", "--peaks", str(peaks), "--out", str(peaks), "--dx", "0.01")
        assert p.returncode == 1
        assert "same file" in (p.stdout + p.stderr)


class TestTheMeasuredShift:
    """diagnose already measures the offset and dead-ends in a report. --auto applies
    what it measured, which is the whole point: the error is rigid and correctable.
    """

    def test_auto_recovers_the_injected_offset(self, tmp_path):
        spec, peaks = _dataset(tmp_path, shift=(0.06, 0.6))
        out = tmp_path / "shifted.txt"
        p = _run("peaks", "shift", "--peaks", str(peaks), "--out", str(out),
                 "--auto", "--spectrum", str(spec), "--format", "json")
        assert p.returncode == 0, p.stderr
        payload = json.loads(p.stdout)
        assert payload["dx"] == pytest.approx(0.06, abs=0.01)
        assert payload["dy"] == pytest.approx(0.6, abs=0.1)

    def test_the_corrected_list_registers_against_the_spectrum(self, tmp_path):
        """The round-trip, and the reason the sign matters: applying the offset with the
        wrong sign doubles the error instead of removing it, and both look like 'it ran'.
        """
        from lunaNMR.validation.spectra_check import assess_registration
        spec, peaks = _dataset(tmp_path, shift=(0.06, 0.6))
        out = tmp_path / "shifted.txt"
        assert _run("peaks", "shift", "--peaks", str(peaks), "--out", str(out),
                    "--auto", "--spectrum", str(spec)).returncode == 0
        before = assess_registration(str(spec), _positions(peaks))
        after = assess_registration(str(spec), _positions(out))
        assert abs(after[0]) < abs(before[0]) / 2
        assert abs(after[1]) < abs(before[1]) / 2
        assert abs(after[0]) <= 0.02 and abs(after[1]) <= 0.2

    def test_auto_needs_a_spectrum(self, tmp_path):
        _, peaks = _dataset(tmp_path)
        p = _run("peaks", "shift", "--peaks", str(peaks),
                 "--out", str(tmp_path / "o.txt"), "--auto")
        assert p.returncode == 2
        assert "--spectrum" in p.stderr

    def test_a_dry_run_reports_without_writing(self, tmp_path):
        spec, peaks = _dataset(tmp_path)
        out = tmp_path / "shifted.txt"
        assert _run("peaks", "shift", "--peaks", str(peaks), "--out", str(out),
                    "--auto", "--spectrum", str(spec), "--dry-run").returncode == 0
        assert not out.exists()


class TestTheMeasurementGridIsTheResolutionLimit:
    """--quick quantises the 15N offset to 0.03 ppm instead of 0.015. A smaller real
    offset cannot be represented, so it is misreported in whichever direction the grid
    lands — on the real shifted list in data_test it came back as 0.000 (missing a real
    0.015 error), and on this synthetic one as 0.030 (twice the truth). Both look like
    a confident measurement.
    """

    def test_quick_quantises_a_sub_step_offset(self, tmp_path):
        spec, peaks = _dataset(tmp_path, shift=(0.0, 0.015))
        p = _run("peaks", "shift", "--peaks", str(peaks), "--out", str(tmp_path / "q.txt"),
                 "--auto", "--spectrum", str(spec), "--quick", "--format", "json")
        dy = json.loads(p.stdout)["dy"]
        assert abs(dy % 0.03) < 1e-6 or abs(dy % 0.03 - 0.03) < 1e-6, dy
        assert dy != pytest.approx(0.015, abs=1e-3), "0.015 is not on the --quick grid"

    def test_the_full_grid_resolves_it(self, tmp_path):
        spec, peaks = _dataset(tmp_path, shift=(0.0, 0.015))
        p = _run("peaks", "shift", "--peaks", str(peaks), "--out", str(tmp_path / "f.txt"),
                 "--auto", "--spectrum", str(spec), "--format", "json")
        assert json.loads(p.stdout)["dy"] == pytest.approx(0.015, abs=0.008)

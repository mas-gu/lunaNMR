# ABOUTME: The T1/T2 page's T1rho -> T2 conversion reads whitespace peak lists and reads its
# ABOUTME: fitted times in the page's delay units. Real page, synthetic results files.

import os
import sys
from pathlib import Path

import pandas as pd
import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

_DIR = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_DIR))

# residue, 15N ppm, 1H ppm, T1 (ms), T1rho (ms)
RESIDUES = [("3LysH", 126.463, 8.179, 712.0, 95.0),
            ("8TyrH", 116.317, 6.929, 655.5, 81.25),
            ("10ValH", 120.288, 8.202, 698.0, 102.5)]


def _results_file(path, experiment, column, per_ms):
    """A fit-results file of `experiment` whose times are the ms values times `per_ms`."""
    lines = [f"Residue\tA\t{experiment}\tC\tA_err\t{experiment}_err\tC_err\tSuccess\tWindowRatio"]
    for row in RESIDUES:
        value = row[column] * per_ms
        lines.append(f"{row[0]}\t1.000000e+06\t{value:.6e}\t0.000000e+00\t1.000000e+03\t"
                     f"{value * 0.01:.6e}\t0.000000e+00\tYes\t3.000")
    path.write_text("\n".join(lines) + "\n")
    return str(path)


@pytest.fixture(scope="module")
def app():
    from PySide6.QtWidgets import QApplication
    return QApplication.instance() or QApplication(sys.argv)


@pytest.fixture
def page(app):
    from dynamiXs_gui import DynamiXsMainWindow
    window = DynamiXsMainWindow()
    yield window.t1t2_page
    window.deleteLater()


@pytest.fixture
def sparky_peaks(tmp_path):
    path = tmp_path / "list.peaks"
    rows = ["      Assignment         w1         w2  ", ""]
    rows += [f"{name:>16} {n:10.3f} {h:10.3f}" for name, n, h, _, _ in RESIDUES]
    path.write_text("\n".join(rows) + "\n")
    return str(path)


def test_a_sparky_peak_list_loads(page, sparky_peaks):
    peaks = page._load_peak_list(sparky_peaks)
    assert list(peaks["residue"]) == [r[0] for r in RESIDUES]
    assert list(peaks["N_ppm"]) == [r[1] for r in RESIDUES]


def test_a_whitespace_text_peak_list_loads(page, tmp_path):
    path = tmp_path / "list.txt"
    rows = ["Assignment  15N  1H"] + [f"{name}  {n}  {h}" for name, n, h, _, _ in RESIDUES]
    path.write_text("\n".join(rows) + "\n")
    peaks = page._load_peak_list(str(path))
    assert list(peaks["residue"]) == [r[0] for r in RESIDUES]
    assert list(peaks["N_ppm"]) == [r[1] for r in RESIDUES]


def _converted(page, tmp_path, sparky_peaks, units, per_ms):
    """The T2 table the page derives from fits whose times are written in `units`."""
    out = tmp_path / units
    out.mkdir()
    page.output_dir = str(out)
    page.peak_list_file = sparky_peaks
    page.units_combo.setCurrentText(units)
    page.field1_carrier_spin.setValue(118.0)
    page.session_results["field1"]["t2"] = None
    page.session_results["field1"]["t1"] = {
        "results_file": _results_file(out / "T1_fit_results.txt", "T1", 3, per_ms)}
    page.session_results["field1"]["t1rho"] = {
        "results_file": _results_file(out / "T1rho_fit_results.txt", "T1rho", 4, per_ms)}
    page._calculate_t2_from_t1rho()
    table = out / f"{page._run_prefix()}_field1_T2_from_T1rho.csv"
    assert table.exists(), page.results_text.toPlainText()
    return pd.read_csv(table)


@pytest.mark.parametrize("units, per_ms", [("s", 1e-3), ("us", 1e3)])
def test_times_are_read_in_the_pages_delay_units(page, tmp_path, sparky_peaks, units, per_ms):
    in_ms = _converted(page, tmp_path, sparky_peaks, "ms", 1.0)
    in_units = _converted(page, tmp_path, sparky_peaks, units, per_ms)

    assert len(in_ms) == len(RESIDUES)
    assert list(in_units["residue"]) == list(in_ms["residue"])
    for column in ("R1", "R1rho", "R2", "T2", "T2_err"):
        assert in_units[column].to_numpy() == pytest.approx(in_ms[column].to_numpy(), rel=1e-6)
    # a 700 ms T1 is a rate near 1.4 per second, whatever unit its fit was written in
    assert in_units["R1"].between(1.0, 2.0).all()

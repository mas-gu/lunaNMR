# ABOUTME: The T1/T2 page's T1rho -> T2 conversion reads whitespace peak lists and reads each
# ABOUTME: fit's times in the units that fit records. Real page; synthetic files and real fits.

import json
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

_DIR = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_DIR))

# residue, 15N ppm, 1H ppm, T1 (ms), T1rho (ms)
RESIDUES = [("3LysH", 126.463, 8.179, 712.0, 95.0),
            ("8TyrH", 116.317, 6.929, 655.5, 81.25),
            ("10ValH", 120.288, 8.202, 698.0, 102.5)]

# what a time in ms is multiplied by to be written in each unit
PER_MS = {"ms": 1.0, "s": 1e-3, "us": 1e3}


def _fit_entry(folder, experiment, column, units):
    """A session entry of an `experiment` fit made in `units`: its results file, and the
    fit JSON recording those units."""
    lines = [f"Residue\tA\t{experiment}\tC\tA_err\t{experiment}_err\tC_err\tSuccess\tWindowRatio"]
    for row in RESIDUES:
        value = row[column] * PER_MS[units]
        lines.append(f"{row[0]}\t1.000000e+06\t{value:.6e}\t0.000000e+00\t1.000000e+03\t"
                     f"{value * 0.01:.6e}\t0.000000e+00\tYes\t3.000")
    results = folder / f"{experiment}_fit_results.txt"
    results.write_text("\n".join(lines) + "\n")
    fit_json = folder / f"{experiment}_fit_data.json"
    fit_json.write_text(json.dumps({"metadata": {"time_units": units}}))
    return {"results_file": str(results), "json_file": str(fit_json)}


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


def _matrix(path, column, delays_ms, units):
    """A series intensity matrix decaying with each residue's time constant, its delay
    headers written in `units`."""
    rows = []
    for number, residue in enumerate(RESIDUES, start=1):
        row = {"Peak_Number": number, "Assignment": residue[0],
               "Reference_X": residue[2], "Reference_Y": residue[1]}
        for delay in delays_ms:
            row[f"{delay * PER_MS[units]:g}"] = 1.0e6 * np.exp(-delay / residue[column])
        rows.append(row)
    pd.DataFrame(rows).to_csv(path, index=False)
    return str(path)


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


def _converted(page, tmp_path, sparky_peaks, t1_units, t1rho_units):
    """The T2 table the page derives from a T1 fit made in `t1_units` and a T1rho fit made
    in `t1rho_units`."""
    out = tmp_path / f"{t1_units}_{t1rho_units}"
    out.mkdir()
    page.output_dir = str(out)
    page.peak_list_file = sparky_peaks
    page.field1_carrier_spin.setValue(118.0)
    page.session_results["field1"]["t2"] = None
    page.session_results["field1"]["t1"] = _fit_entry(out, "T1", 3, t1_units)
    page.session_results["field1"]["t1rho"] = _fit_entry(out, "T1rho", 4, t1rho_units)
    page._calculate_t2_from_t1rho()
    table = out / f"{page._run_prefix()}_field1_T2_from_T1rho.csv"
    assert table.exists(), page.results_text.toPlainText()
    return pd.read_csv(table)


@pytest.mark.parametrize("t1_units, t1rho_units", [("s", "s"), ("us", "us"), ("s", "ms")])
def test_times_are_read_in_the_units_each_fit_records(
        page, tmp_path, sparky_peaks, t1_units, t1rho_units):
    in_ms = _converted(page, tmp_path, sparky_peaks, "ms", "ms")
    in_units = _converted(page, tmp_path, sparky_peaks, t1_units, t1rho_units)

    assert len(in_ms) == len(RESIDUES)
    assert list(in_units["residue"]) == list(in_ms["residue"])
    for column in ("R1", "R1rho", "R2", "T2", "T2_err"):
        assert in_units[column].to_numpy() == pytest.approx(in_ms[column].to_numpy(), rel=1e-6)
    # a 700 ms T1 is a rate near 1.4 per second, whatever unit its fit was written in
    assert in_units["R1"].between(1.0, 2.0).all()


def _fitted_then_converted(app, page, tmp_path, sparky_peaks, units, selector_after):
    """The T2 table the page derives from T1 and T1rho series it fitted in `units`, the
    delay-units selector moved to `selector_after` between the fits and the conversion."""
    from test_t1t2_page_all_runs import _run_to_the_end
    out = tmp_path / units
    out.mkdir()
    page.output_dir = str(out)
    page.peak_list_file = sparky_peaks
    page.field1_carrier_spin.setValue(118.0)
    page.session_results["field1"]["t2"] = None
    page.field1_t1_file = _matrix(out / "T1.csv", 3, [0, 50, 100, 200, 400, 800, 1200], units)
    page.field1_t1rho_file = _matrix(out / "T1rho.csv", 4, [0, 10, 20, 40, 80, 160, 240], units)
    page.units_combo.setCurrentText(units)
    for selection in ("Field 1 - T1", "Field 1 - T1ρ"):
        page.analysis_combo.setCurrentText(selection)
        page._start_analysis()
        _run_to_the_end(app, page)
    page.units_combo.setCurrentText(selector_after)
    page._calculate_t2_from_t1rho()
    table = out / f"{page._run_prefix()}_field1_T2_from_T1rho.csv"
    assert table.exists(), page.results_text.toPlainText()
    return pd.read_csv(table)


def test_a_fit_converts_in_its_own_units_wherever_the_selector_stands(
        app, page, tmp_path, sparky_peaks):
    in_ms = _fitted_then_converted(app, page, tmp_path, sparky_peaks, "ms", "s")
    in_s = _fitted_then_converted(app, page, tmp_path, sparky_peaks, "s", "ms")

    assert len(in_ms) == len(RESIDUES)
    # a 700 ms T1 is a rate near 1.4 per second, a 90 ms T1rho one near 11
    for table in (in_ms, in_s):
        assert table["R1"].between(1.0, 2.0).all()
        assert table["R1rho"].between(9.0, 13.0).all()
    for column in ("R1", "R1rho", "R2", "T2"):
        assert in_s[column].to_numpy() == pytest.approx(in_ms[column].to_numpy(), rel=1e-4)

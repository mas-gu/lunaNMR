# ABOUTME: The Kd survey reports as selected the residues its selection file selects,
# ABOUTME: also when the file keeps rulings a person made after an earlier survey.

import re
import sys
from pathlib import Path

import numpy as np
import pytest

_KD_DIR = Path(__file__).resolve().parent.parent / "dynamiXs_Kd"
sys.path.insert(0, str(_KD_DIR))

LIGAND = [0.0, 12.5, 25.0, 50.0, 100.0, 150.0]


def _tidy(path):
    """A six-point titration of five residues moving and decaying with the ligand."""
    lines = ["spectrum_name,assignment,ppm_x,ppm_y,height,volume"]
    for n, name in enumerate(("3LysH", "8TyrH", "10ValH", "11LeuH", "14GlyH"), start=1):
        for point, ligand in enumerate(LIGAND):
            bound = ligand / (ligand + 60.0)
            lines.append(f"{point},{name},{7.5 + 0.1 * n + 0.02 * n * bound:.4f},"
                         f"{110.0 + n + 0.1 * n * bound:.4f},"
                         f"{1000.0 * n * np.exp(-ligand / 80.0):.1f},"
                         f"{10000.0 * n * np.exp(-ligand / 80.0):.1f}")
    path.write_text("\n".join(lines) + "\n")
    return str(path)


def _survey(tmp_path):
    from kd_survey import run_kd_survey_with_params
    return run_kd_survey_with_params({
        "input_csv_file": str(tmp_path / "series_analysis_tidy.csv"),
        "output_dir": str(tmp_path / "survey"), "output_prefix": "kd",
        "concentrations": LIGAND, "protein_conc": 50.0})


def _stated_in_file(text):
    return int(re.search(r"# (\d+) of \d+ residues selected", text).group(1))


@pytest.mark.filterwarnings("ignore:Covariance of the parameters could not be estimated")
def test_the_count_follows_the_selection_file_after_a_person_edited_it(tmp_path):
    from kd_survey import parse_residues_file
    _tidy(tmp_path / "series_analysis_tidy.csv")
    first = _survey(tmp_path)
    selection = Path(first["residues_file"])
    assert first["n_selected"] == _stated_in_file(selection.read_text()) == 5

    # two residues ruled out by hand, as the file's own header tells its reader to
    edited = [("# " + line if line.split()[:1] in (["8TyrH"], ["11LeuH"]) else line)
              for line in selection.read_text().splitlines()]
    selection.write_text("\n".join(edited) + "\n")

    second = _survey(tmp_path)

    text = selection.read_text()
    assert parse_residues_file(text) == ["3LysH", "10ValH", "14GlyH"]
    assert _stated_in_file(text) == 3
    assert second["n_selected"] == 3

# ABOUTME: `export kd` blanks a residue's ref->point intensity by the same reference rule
# ABOUTME: the fit applied, including a --ref-max-ratio other than the default.

import contextlib
import csv
import io
import json
import shutil
from pathlib import Path

import numpy as np
import pytest

TIDY = Path(__file__).resolve().parents[2] / "data_test" / "data_binding" / "tidy"
ABSOLUTE = "0,6.25,12.5,25,50,100,150"

pytestmark = pytest.mark.skipif(not TIDY.exists(), reason="data_test not present")


def _run(*args):
    from lunaNMR.cli import main
    out = io.StringIO()
    with contextlib.redirect_stdout(out):
        code = main([*args, "--format", "json"])
    assert code == 0
    return json.loads(out.getvalue())


def test_the_exported_ratios_keep_the_residues_the_fit_kept(tmp_path):
    tidy = tmp_path / "DNAJA1_HSPA8_tidy.csv"
    shutil.copy(TIDY / "DNAJA1_HSPA8_tidy.csv", tidy)
    fit = _run("kd", "--input", str(tidy), "--out", str(tmp_path / "fit"), "--p0", "50",
               "--conc", ABSOLUTE, "--observable", "intensity", "--ref-max-ratio", "100")
    data = json.loads(Path(fit["json_file"]).read_text())

    # Residues whose reference sits 10-100x below their own maximum: dropped at the
    # default cap of 10, kept by the cap this fit was run with.
    kept_by_this_cap = []
    for f in data["fits"]:
        h = np.array(f["series"]["height"], dtype=float)
        if 10.0 < np.nanmax(h) / h[0] <= 100.0:
            assert f["intensity"]["success"], f["residue"]
            kept_by_this_cap.append(f["residue"])
    assert kept_by_this_cap

    _run("export", "kd", "--json", fit["json_file"], "--out", str(tmp_path / "export"),
         "--kind", "ref-bars", "--observable", "intensity")
    table = tmp_path / "export" / "data" / "intensity_ref_vs_point.csv"
    rows = {r[0]: r[1:] for r in csv.reader(table.open())}
    for name in kept_by_this_cap:
        assert any(cell != "" for cell in rows[name]), name

# ABOUTME: Two t1t2 / methyl-t2 runs sharing an --out must keep a fit JSON each.
# ABOUTME: A JSON named without the prefix let the second run overwrite the first's record.

"""
Each run's `<prefix>_dynamixs.json` names its fit JSON, and `project build` copies that
file into the bundle. When the name came from --field-name alone (default 'field1'),
a 600 MHz and a 700 MHz run into one folder wrote the same file, and the 600 MHz record
then pointed at the 700 MHz fit.
"""

import json
import os
from pathlib import Path

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

DATA = Path(__file__).resolve().parents[2] / "data_test"
SERIES = DATA / "data_JA1_relaxation" / "run_20260825_1405" / "series"
METHYL = (DATA / "methyl" / "sampleNBD" / "series_results_sampleNBD_20260818_124619"
          / "peak_intensity_matrix.csv")

pytestmark = pytest.mark.skipif(not (SERIES.exists() and METHYL.exists()),
                                reason="data_test not present")


def _run(argv, capsys):
    from lunaNMR.cli import main
    assert main(argv + ["--format", "json"]) == 0
    return json.loads(capsys.readouterr().out)


def _assert_each_run_keeps_its_own(summaries):
    files = [s['json_file'] for s in summaries.values()]
    assert len(set(files)) == len(files), files
    for freq, summary in summaries.items():
        meta = json.loads(Path(summary['json_file']).read_text())['metadata']
        assert float(meta['field_freq']) == freq


def test_t1t2_runs_in_one_folder_keep_their_own_json(tmp_path, capsys):
    summaries = {}
    for freq in (600.0, 700.0):
        summaries[freq] = _run(
            ["dynamixs", "t1t2", "--exp", "T1", "--out", str(tmp_path),
             "--input", str(SERIES / f"{int(freq)}_T1" / "peak_intensity_matrix.csv"),
             "--prefix", f"{int(freq)}_T1", "--field-freq", str(freq)], capsys)
    _assert_each_run_keeps_its_own(summaries)


def test_methyl_runs_in_one_folder_keep_their_own_json(tmp_path, capsys):
    summaries = {}
    for freq in (600.0, 800.0):
        summaries[freq] = _run(
            ["dynamixs", "methyl-t2", "--out", str(tmp_path), "--input", str(METHYL),
             "--prefix", f"run{int(freq)}", "--field-freq", str(freq)], capsys)
    _assert_each_run_keeps_its_own(summaries)

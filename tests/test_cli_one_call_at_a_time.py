# ABOUTME: cli.main() called from several threads of one process runs the calls one at a time.
# ABOUTME: Overlapping calls swapped fd 1 under each other and sent a JSON summary to stderr.

import contextlib
import json
import threading
from pathlib import Path

import pytest

from lunaNMR import cli

SERIES = (Path(__file__).resolve().parents[2] / "data_test" / "data_JA1_relaxation"
          / "run_20260825_1405" / "series")

pytestmark = pytest.mark.skipif(not SERIES.exists(), reason="data_test not present")


def _argv(exp, out):
    return ["dynamixs", "t1t2", "--input", str(SERIES / f"600_{exp}" / "peak_intensity_matrix.csv"),
            "--out", str(out), "--exp", exp, "--field-freq", "600", "--prefix", "t",
            "--format", "json"]


def _comparable(summary):
    """A summary without the keys that name where it was written."""
    return {k: v for k, v in summary.items()
            if k not in ('results_file', 'json_file')}


def test_two_threads_get_their_own_summaries_on_stdout(tmp_path, capfd, monkeypatch):
    expected = {}
    for exp in ("T1", "T2"):
        assert cli.main(_argv(exp, tmp_path / f"seq_{exp}")) == 0
        expected[exp] = _comparable(json.loads(capfd.readouterr().out))

    active, peak = [0], [0]
    count = threading.Lock()
    engine_stdout = cli._engine_stdout

    def counted(args):
        with count:
            active[0] += 1
            peak[0] = max(peak[0], active[0])
        try:
            with engine_stdout(args):
                yield
        finally:
            with count:
                active[0] -= 1

    monkeypatch.setattr(cli, "_engine_stdout", contextlib.contextmanager(counted))

    codes = {}
    start = threading.Barrier(2)

    def call(exp):
        start.wait()
        codes[exp] = cli.main(_argv(exp, tmp_path / f"par_{exp}"))

    threads = [threading.Thread(target=call, args=(exp,)) for exp in ("T1", "T2")]
    for t in threads:
        t.start()
    for t in threads:
        t.join()

    captured = capfd.readouterr()
    assert codes == {"T1": 0, "T2": 0}
    assert peak[0] == 1
    summaries = [json.loads(line) for line in captured.out.splitlines() if line.strip()]
    assert sorted(s["experiment"] for s in summaries) == ["T1", "T2"]
    for s in summaries:
        assert _comparable(s) == expected[s["experiment"]]
    assert '"ok"' not in captured.err

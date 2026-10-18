# ABOUTME: An unexpected exception still answers --format json with an ok:false envelope.
# ABOUTME: Marked internal, so a caller can tell a bug from bad input; the traceback goes to stderr.

"""
Expected bad-input failures were reported in the envelope; anything else - a bug - escaped
main() as a traceback with nothing on stdout, so a caller parsing stdout got a JSON error
about the empty string instead of the cause.
"""

import json


def _export(tmp_path, capsys, payload):
    from lunaNMR.cli import main
    fit = tmp_path / "fit.json"
    fit.write_text(payload)
    code = main(["export", "kd", "--json", str(fit), "--out", str(tmp_path / "o"),
                 "--format", "json"])
    return code, capsys.readouterr()


def test_a_bug_is_reported_in_the_envelope(tmp_path, capsys):
    code, out = _export(tmp_path, capsys, "[1, 2, 3]")
    assert code == 1
    envelope = json.loads(out.out)
    assert envelope["ok"] is False
    assert envelope["error"]["internal"] is True
    assert envelope["error"]["type"] == "AttributeError"
    assert "Traceback" in out.err


def test_bad_input_is_not_marked_internal(tmp_path, capsys):
    from lunaNMR.cli import main
    code = main(["kd", "--input", str(tmp_path / "missing.csv"), "--out", str(tmp_path / "o"),
                 "--p0", "50", "--conc-units", "absolute", "--format", "json"])
    envelope = json.loads(capsys.readouterr().out)
    assert code == 1
    assert "internal" not in envelope["error"]

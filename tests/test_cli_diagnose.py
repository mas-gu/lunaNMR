# ABOUTME: CLI surface for the pre-flight checks - `series --dry-run --deep` and `lunaNMR diagnose`.
# ABOUTME: The plain dry-run contract must not change: same keys, no spectrum reads.

import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
ng = pytest.importorskip("nmrglue")

from test_spectra_check import _folder, _peaks   # the synthetic-spectrum helpers

ROOT = Path(__file__).resolve().parent.parent


def _run(*argv):
    p = subprocess.run([sys.executable, '-m', 'lunaNMR', *argv],
                       cwd=str(ROOT), capture_output=True, text=True)
    return p


def _dataset(tmp_path):
    peaks = _peaks()
    _folder(tmp_path, 'T1', {'T1_50ms.ft': {}, 'T1_100ms.ft': {}}, peaks)
    _folder(tmp_path, 'T2', {'T2_8ms.ft': {}, 'T2_17ms.ft': {}}, peaks)
    return tmp_path


class TestPlainDryRunUnchanged:
    """Existing callers parse these keys and rely on it being instant."""

    def test_keys_are_preserved(self, tmp_path):
        d = _dataset(tmp_path)
        p = _run('series', '--spectra', str(d / 'T1'), '--peaks', str(d / 'T1/T1_assignment.txt'),
                 '--out', str(tmp_path / 'out'), '--format', 'json', '--dry-run')
        assert p.returncode == 0, p.stderr
        got = json.loads(p.stdout)
        for key in ('command', 'dry_run', 'spectra_found', 'peaks', 'peaks_exists',
                    'mode', 'peak_source', 'parallel', 'output_dir', 'missing_inputs'):
            assert key in got, f'plain --dry-run lost the {key} key'

    def test_does_not_run_the_checks(self, tmp_path):
        d = _dataset(tmp_path)
        p = _run('series', '--spectra', str(d / 'T1'), '--peaks', str(d / 'T1/T1_assignment.txt'),
                 '--out', str(tmp_path / 'out'), '--format', 'json', '--dry-run')
        assert 'checks' not in json.loads(p.stdout)


class TestDeep:
    def test_adds_checks_without_disturbing_the_plain_keys(self, tmp_path):
        d = _dataset(tmp_path)
        p = _run('series', '--spectra', str(d / 'T1'), '--peaks', str(d / 'T1/T1_assignment.txt'),
                 '--out', str(tmp_path / 'out'), '--format', 'json', '--dry-run', '--deep')
        assert p.returncode == 0, p.stderr
        got = json.loads(p.stdout)
        assert got['spectra_found'] == 2
        checks = got['checks']
        assert len(checks['spectra']) == 2
        assert checks['delays']['n_parsed'] == 2

    def test_requires_dry_run(self, tmp_path):
        d = _dataset(tmp_path)
        p = _run('series', '--spectra', str(d / 'T1'), '--peaks', str(d / 'T1/T1_assignment.txt'),
                 '--out', str(tmp_path / 'out'), '--deep')
        assert p.returncode != 0
        assert '--deep' in p.stderr

    def test_findings_are_serialisable(self, tmp_path):
        """A mis-registered list must survive the trip through JSON."""
        peaks = _peaks()
        d = tmp_path / 'T1'
        from test_spectra_check import _spectrum, _write_list
        d.mkdir()
        _spectrum(d / 'T1_50ms.ft', peaks, shift=(0.06, 0.0))
        _write_list(d / 'T1_assignment.txt', peaks)
        p = _run('series', '--spectra', str(d), '--peaks', str(d / 'T1_assignment.txt'),
                 '--out', str(tmp_path / 'out'), '--format', 'json', '--dry-run', '--deep')
        # 1, not 0: the registration FAIL below is exactly what --deep is for, and a
        # dry-run that reports one while exiting 0 does not gate anything. The JSON is
        # still emitted, which is what this test is about.
        assert p.returncode == 1, p.stderr
        findings = json.loads(p.stdout)['checks']['findings']
        assert any(f['severity'] == 'FAIL' and f['check'] == 'registration' for f in findings)


class TestDiagnoseSubcommand:
    def test_reports_every_experiment(self, tmp_path):
        d = _dataset(tmp_path)
        p = _run('diagnose', str(d), '--quick', '--format', 'json')
        assert p.returncode == 0, p.stderr
        got = json.loads(p.stdout)
        assert sorted(e['experiment'] for e in got['experiments']) == ['T1', 'T2']

    def test_compares_assignment_sets_across_experiments(self, tmp_path):
        """The check a single-experiment --deep structurally cannot do."""
        _folder(tmp_path, 'T1', {'T1_50ms.ft': {}}, _peaks(6))
        _folder(tmp_path, 'T2', {'T2_8ms.ft': {}}, _peaks(6)[:5])
        p = _run('diagnose', str(tmp_path), '--quick', '--format', 'json')
        assert p.returncode == 0, p.stderr
        got = json.loads(p.stdout)
        assert got['assignments']['common'] == 5
        assert got['assignments']['union'] == 6

    def test_text_output_is_human_readable(self, tmp_path):
        d = _dataset(tmp_path)
        p = _run('diagnose', str(d), '--quick')
        assert p.returncode == 0, p.stderr
        assert 'T1' in p.stdout and 'T2' in p.stdout

    def test_missing_root_is_an_error(self, tmp_path):
        p = _run('diagnose', str(tmp_path / 'nope'), '--format', 'json')
        assert p.returncode == 1


class TestTheGateActuallyGates:
    """`diagnose && series` is the documented pre-flight. It only gates if a FAIL
    changes the exit code -- otherwise the shell sails straight past the failure the
    command exists to catch.
    """

    def _misregistered(self, tmp_path):
        """A peak list belonging to a different frame: the failure diagnose is for."""
        from test_spectra_check import _spectrum, _write_list
        peaks = _peaks()
        d = tmp_path / 'T1'
        d.mkdir()
        _spectrum(d / 'T1_50ms.ft', peaks, shift=(0.06, 0.0))
        _write_list(d / 'T1_assignment.txt', peaks)
        return tmp_path

    def test_a_fail_exits_nonzero(self, tmp_path):
        p = _run('diagnose', str(self._misregistered(tmp_path)), '--quick')
        assert p.returncode == 1, p.stdout

    def test_a_clean_dataset_exits_zero(self, tmp_path):
        p = _run('diagnose', str(_dataset(tmp_path)), '--quick')
        assert p.returncode == 0, p.stderr

    def test_warnings_alone_do_not_fail(self, tmp_path):
        """WARN includes 'no peak list in this folder', routine in a mixed dataset
        root. Failing on it makes the gate unusable, so it gets removed instead."""
        _folder(tmp_path, 'T1', {'T1_50ms.ft': {}}, _peaks(6))
        _folder(tmp_path, 'T2', {'T2_8ms.ft': {}}, _peaks(6)[:5])
        p = _run('diagnose', str(tmp_path), '--quick', '--format', 'json')
        assert p.returncode == 0, p.stderr
        assert any(f['severity'] == 'WARN' for f in json.loads(p.stdout)['findings'])

    def test_strict_fails_on_warnings_too(self, tmp_path):
        _folder(tmp_path, 'T1', {'T1_50ms.ft': {}}, _peaks(6))
        _folder(tmp_path, 'T2', {'T2_8ms.ft': {}}, _peaks(6)[:5])
        p = _run('diagnose', str(tmp_path), '--quick', '--strict')
        assert p.returncode == 1, p.stdout

    def test_the_deep_dry_run_gates_the_same_way(self, tmp_path):
        """--deep runs the same checks, so a FAIL there means the same thing."""
        d = self._misregistered(tmp_path)
        p = _run('series', '--spectra', str(d / 'T1'),
                 '--peaks', str(d / 'T1/T1_assignment.txt'),
                 '--out', str(tmp_path / 'out'), '--dry-run', '--deep')
        assert p.returncode == 1, p.stdout


class TestDiagnoseDoesNotAdvertiseADryRun:
    """diagnose reads spectra, fits nothing and writes nothing. A --dry-run of a
    read-only command has nothing to mean, and it was accepted and ignored."""

    def test_dry_run_is_rejected(self, tmp_path):
        p = _run('diagnose', str(_dataset(tmp_path)), '--dry-run')
        assert p.returncode == 2
        assert 'unrecognized arguments' in p.stderr

    def test_format_is_still_accepted(self, tmp_path):
        p = _run('diagnose', str(_dataset(tmp_path)), '--quick', '--format', 'json')
        assert p.returncode == 0, p.stderr

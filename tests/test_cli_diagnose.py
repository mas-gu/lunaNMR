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
        assert p.returncode == 0, p.stderr
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

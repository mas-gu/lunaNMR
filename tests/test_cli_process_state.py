# ABOUTME: cli.main() called repeatedly in one long-lived process leaves no state behind.
# ABOUTME: Logging handlers, signal handlers and similar process-wide settings survive a call unchanged.

import logging
import signal
import threading
from pathlib import Path

import pytest

from lunaNMR import cli

_BATCH_LOGGER = 'lunaNMR.batch_processing.batch_processor'
_BATCH_PACKAGE = 'lunaNMR.batch_processing'


def _process_state():
    root = logging.getLogger()
    return {'root_handlers': list(root.handlers), 'root_level': root.level,
            'batch_handlers': list(logging.getLogger(_BATCH_LOGGER).handlers),
            'batch_propagate': logging.getLogger(_BATCH_LOGGER).propagate,
            'batch_package_handlers': list(logging.getLogger(_BATCH_PACKAGE).handlers),
            'batch_package_propagate': logging.getLogger(_BATCH_PACKAGE).propagate,
            'sigint': signal.getsignal(signal.SIGINT),
            'sigterm': signal.getsignal(signal.SIGTERM)}


class TestBatchLeavesTheProcessAsItFoundIt:
    """`batch` added a stdout handler to the root logger on every call and replaced the
    host's SIGINT/SIGTERM handlers with sys.exit(0). After one call, every later log
    record in the process landed on stdout ahead of a --format json summary, and a
    SIGTERM to the host exited it with status 0."""

    def test_list_presets_twice(self, capsys):
        before = _process_state()
        assert cli.main(['batch', '--list-presets']) == 0
        assert cli.main(['batch', '--list-presets']) == 0
        assert _process_state() == before

    def test_a_dry_run_over_a_folder(self, tmp_path, monkeypatch, capsys):
        """The dry run builds a BatchProcessor, which adds its own console handler and a
        log file in the working directory."""
        monkeypatch.chdir(tmp_path)
        (tmp_path / 'spectra').mkdir()
        # One real spectrum: a dry run that finds none is a failure
        dnaja1 = Path(__file__).resolve().parents[2] / 'data_test' / 'data_binding' / 'DNAJA1_HSPA8'
        if not dnaja1.exists():
            pytest.skip('data_test not present')
        (tmp_path / 'spectra' / '2D_DNAJA1_0o0.ft').symlink_to(dnaja1 / '2D_DNAJA1_0o0.ft')
        before = _process_state()
        assert cli.main(['batch', str(tmp_path / 'spectra'), '--dry-run']) == 0
        assert _process_state() == before

    def test_a_later_log_record_stays_off_stdout(self, capsys):
        cli.main(['batch', '--list-presets'])
        capsys.readouterr()
        logging.getLogger('some.host.module').warning('host warning')
        assert 'host warning' not in capsys.readouterr().out

    def test_it_runs_off_the_main_thread(self, capsys):
        """signal.signal raises ValueError outside the main thread, so a host running
        tool calls in a worker thread could not run `batch` at all."""
        outcome = {}

        def call():
            try:
                outcome['rc'] = cli.main(['batch', '--list-presets'])
            except BaseException as exc:
                outcome['exc'] = exc

        worker = threading.Thread(target=call)
        worker.start()
        worker.join()
        assert outcome == {'rc': 0}


class TestAUsageErrorReturnsFromMain:
    """argparse ends a bad command line with SystemExit, a BaseException: a host that
    calls main() and catches Exception was ended by its caller's typo."""

    def test_a_bad_choice_returns_2(self, tmp_path, capsys):
        assert cli.main(['kd', '--input', 'in.csv', '--out', str(tmp_path / 'o'),
                         '--p0', '50', '--observable', 'cps']) == 2

    def test_a_validator_rejection_returns_2(self, tmp_path, capsys):
        assert cli.main(['series', '--spectra', str(tmp_path), '--peaks', 'p.txt',
                         '--out', str(tmp_path / 'o'), '--deep']) == 2

    def test_help_returns_0(self, capsys):
        assert cli.main(['--help']) == 0
        assert 'usage: lunaNMR' in capsys.readouterr().out


_TIDY = ("spectrum_name,assignment,ppm_x,ppm_y,height,volume\n"
         "0,K14,7.5,110.0,200.0,2000.0\n1,K14,7.51,110.05,150.0,1500.0\n"
         "2,K14,7.53,110.2,90.0,900.0\n0,A17,8.0,120.0,100.0,1000.0\n"
         "1,A17,8.02,120.1,80.0,800.0\n2,A17,8.05,120.3,50.0,500.0\n")


class TestMatplotlibSettingsAreRestored:
    """`kd --survey` and `export kd` set pdf.fonttype for the whole process, so every
    later PDF the host wrote depended on which command had run before it."""

    def _kd_inputs(self, tmp_path):
        csv = tmp_path / 'series_analysis_tidy.csv'
        csv.write_text(_TIDY)
        return ['--input', str(csv), '--p0', '50', '--conc', '0,25,50']

    def test_export_kd(self, tmp_path, capsys):
        import matplotlib
        assert cli.main(['kd', *self._kd_inputs(tmp_path), '--out', str(tmp_path / 'kd')]) == 0
        with matplotlib.rc_context({'pdf.fonttype': 3}):
            assert cli.main(['export', 'kd', '--json',
                             str(tmp_path / 'kd' / 'data' / 'kd_kd_fit_data.json'),
                             '--out', str(tmp_path / 'figs')]) == 0
            assert matplotlib.rcParams['pdf.fonttype'] == 3

    def test_kd_survey(self, tmp_path, capsys):
        import matplotlib
        with matplotlib.rc_context({'pdf.fonttype': 3}):
            assert cli.main(['kd', *self._kd_inputs(tmp_path), '--out', str(tmp_path / 'kd'),
                             '--survey']) == 0
            assert matplotlib.rcParams['pdf.fonttype'] == 3

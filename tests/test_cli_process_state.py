# ABOUTME: cli.main() called repeatedly in one long-lived process leaves no state behind.
# ABOUTME: Logging handlers, signal handlers and similar process-wide settings survive a call unchanged.

import logging
import signal
import threading

from lunaNMR import cli

_BATCH_LOGGER = 'lunaNMR.batch_processing.batch_processor'


def _process_state():
    root = logging.getLogger()
    return {'root_handlers': list(root.handlers), 'root_level': root.level,
            'batch_handlers': list(logging.getLogger(_BATCH_LOGGER).handlers),
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

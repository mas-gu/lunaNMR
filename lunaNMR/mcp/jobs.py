# ABOUTME: Background jobs for the MCP server: a CLI run in a detached process, its state kept on disk.
# ABOUTME: Status is read from files, so a job can be read by a later server; only the starter can cancel it.

import asyncio
import json
import os
import signal
import subprocess
import sys
import threading
import time
from pathlib import Path

from lunaNMR.mcp import summaries
from lunaNMR.mcp.data_root import PathRefused
from lunaNMR.mcp.cli_call import trimmed

JOB_FOLDER = 'mcp_job'
PROGRESS_LINES = 3
LINE_CHARS = 300
STDERR_TAIL_LINES = 40
KILL_GRACE_S = 5

ENRICHMENTS = {'series': ('series_quality', summaries.series_quality),
               'kd': ('kd_quality', summaries.kd_quality),
               't1t2': ('relaxation_quality', summaries.relaxation_quality),
               'hetnoe': ('hetnoe_quality', summaries.hetnoe_quality)}


class NotAJob(ValueError):
    """The folder holds no job record."""


class TooManyJobs(RuntimeError):
    """The server is already running as many jobs as it allows."""


class NotCancellable(RuntimeError):
    """The job was not started by this server, so its process cannot be told apart from a reused id."""


def _write_json(path, content):
    temporary = path.with_name(path.name + '.tmp')
    temporary.write_text(json.dumps(content))
    os.replace(temporary, path)


def _read_json(path):
    try:
        return json.loads(path.read_text())
    except (OSError, ValueError):
        return None


def _is_number(value):
    return isinstance(value, (int, float)) and not isinstance(value, bool)


def _read_record(path):
    """The job record, or None when the file is missing or is not one this server wrote."""
    record = _read_json(path)
    if (isinstance(record, dict) and isinstance(record.get('kind'), str)
            and isinstance(record.get('pid'), int) and _is_number(record.get('started'))):
        return record
    return None


def _read_result(path):
    result = _read_json(path)
    if isinstance(result, dict) and isinstance(result.get('exit_code'), int) \
            and _is_number(result.get('finished')):
        return result
    return None


def _tail(path, lines):
    try:
        text = path.read_text(errors='replace')
    except OSError:
        return []
    return [line[:LINE_CHARS] for line in text.splitlines()[-lines:]]


def _process_exists(pid):
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True
    return True


class JobManager:
    def __init__(self, max_running_jobs=2):
        self.max_running_jobs = max_running_jobs
        self._started = {}

    def _running_here(self):
        return [folder for folder, process in self._started.items() if process.poll() is None]

    def start(self, kind, cli_args, out_folder, files=None):
        """Create the output folder, launch the command detached, and return the job's path.

        `files` maps names to text written into the job folder first, for inputs the command
        must read from a file this server owns.
        """
        if len(self._running_here()) >= self.max_running_jobs:
            raise TooManyJobs(
                f'{len(self._running_here())} jobs already running; wait for one or cancel it')
        job_folder = Path(out_folder) / JOB_FOLDER
        try:
            job_folder.mkdir(parents=True)
        except FileExistsError:
            raise PathRefused(f'already exists, choose a new output folder: {out_folder}')
        for name, text in (files or {}).items():
            (job_folder / name).write_text(text)
        process = subprocess.Popen(
            [sys.executable, '-m', 'lunaNMR.mcp.job_runner', str(job_folder), *cli_args],
            stdin=subprocess.DEVNULL, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
            start_new_session=True)
        threading.Thread(target=process.wait, daemon=True).start()
        _write_json(job_folder / 'job.json', {
            'kind': kind, 'args': cli_args, 'pid': process.pid, 'started': time.time()})
        self._started[Path(out_folder)] = process
        return Path(out_folder)

    def status(self, out_folder):
        job_folder = Path(out_folder) / JOB_FOLDER
        record = _read_record(job_folder / 'job.json')
        if record is None:
            raise NotAJob(f'no job record in {out_folder}')
        status = {'job': str(out_folder), 'kind': record['kind'], 'started': record['started']}
        result = _read_result(job_folder / 'result.json')
        if (job_folder / 'cancelled.json').exists():
            status['status'] = 'cancelled'
        elif result is not None:
            status['status'] = 'done' if result['exit_code'] == 0 else 'failed'
        elif self._is_alive(Path(out_folder), record['pid']):
            status['status'] = 'running'
            status['progress'] = _tail(job_folder / 'stderr.log', PROGRESS_LINES)
        else:
            status['status'] = 'failed'
            status['note'] = 'the job process ended without recording a result'
        finished = result['finished'] if result else time.time()
        status['elapsed_s'] = round(finished - record['started'], 1)
        if status['status'] in ('done', 'failed'):
            self._add_outcome(status, record['kind'], Path(out_folder), job_folder)
        return status

    def _is_alive(self, out_folder, pid):
        process = self._started.get(out_folder)
        if process is not None:
            return process.poll() is None
        return _process_exists(pid)

    @staticmethod
    def _add_outcome(status, kind, out_folder, job_folder):
        envelope = _read_json(job_folder / 'stdout.json')
        if envelope is not None:
            status['result'] = trimmed(envelope)
        if status['status'] == 'failed':
            status['stderr_tail'] = _tail(job_folder / 'stderr.log', STDERR_TAIL_LINES)
        else:
            if kind in ENRICHMENTS:
                key, enrich = ENRICHMENTS[kind]
                status[key] = enrich(out_folder)
            status['files'] = summaries.output_files(out_folder)

    async def cancel(self, out_folder):
        out_folder = Path(out_folder)
        job_folder = out_folder / JOB_FOLDER
        if _read_record(job_folder / 'job.json') is None:
            raise NotAJob(f'no job record in {out_folder}')
        process = self._started.get(out_folder)
        if process is None:
            raise NotCancellable('this server did not start the job')
        if process.poll() is None:
            _write_json(job_folder / 'cancelled.json', {'cancelled': time.time()})
            os.killpg(process.pid, signal.SIGTERM)
            for _ in range(KILL_GRACE_S * 10):
                if process.poll() is not None:
                    break
                await asyncio.sleep(0.1)
            else:
                os.killpg(process.pid, signal.SIGKILL)
        return self.status(out_folder)

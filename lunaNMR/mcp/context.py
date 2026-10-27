# ABOUTME: What every MCP tool shares: the data root, the job manager, and the wait-then-job rule.
# ABOUTME: Tools raise the refusal types below and the server turns them into error results.

import asyncio
import time

from lunaNMR.mcp import cli_call
from lunaNMR.mcp.data_root import PathRefused
from lunaNMR.mcp.jobs import NotAJob, NotCancellable, TooManyJobs
from lunaNMR.mcp.names import InvalidName
from lunaNMR.mcp.results import result

POLL_SECONDS = 0.25

REFUSALS = (PathRefused, TooManyJobs, NotAJob, NotCancellable, InvalidName,
            cli_call.InvalidResidueName)


class ToolContext:
    def __init__(self, root, jobs, wait_seconds):
        self.root = root
        self.jobs = jobs
        self.wait_seconds = wait_seconds

    def job_envelope(self, status):
        status['job'] = str(self.root.relative(status['job']))
        return result({'ok': True, **status})

    async def run(self, args):
        return result(await cli_call.run_cli(args))

    async def launch(self, kind, cli_args, out_folder, dry_run, files=None):
        """A dry run only validates, so it answers at once; a real run is a job."""
        if dry_run:
            return await self.run(cli_args)
        return await self.start_job(kind, cli_args, out_folder, files=files)

    async def start_job(self, kind, cli_args, out_folder, files=None):
        """Start a job and wait up to `wait_seconds` for it: a quick analysis answers in one
        call, a slow one returns while still running and is polled with `job_status`."""
        started = self.jobs.start(kind, cli_args, out_folder, files=files)
        deadline = time.monotonic() + self.wait_seconds
        while True:
            status = self.jobs.status(started)
            if status['status'] != 'running' or time.monotonic() >= deadline:
                return self.job_envelope(status)
            await asyncio.sleep(POLL_SECONDS)

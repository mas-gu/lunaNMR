# ABOUTME: MCP tools to follow and stop the background jobs the other tools start.
# ABOUTME: Status is read from files, so it works for a job an earlier server started.

from mcp.types import CallToolResult

from lunaNMR.mcp.context import REFUSALS
from lunaNMR.mcp.results import CANCELS_A_JOB, READS_ONLY, failure


def register(server, ctx):
    @server.tool(annotations=READS_ONLY)
    async def job_status(job: str) -> CallToolResult:
        """State of a job started by one of the tools that run an analysis: running (with
        the last log lines), done, failed or cancelled. A finished job carries the command's
        summary and what its exit code hides: pinned and undetected points for a series, the
        shared Kd and its pool for a Kd fit. Works for jobs from an earlier server run."""
        try:
            return ctx.job_envelope(ctx.jobs.status(ctx.root.existing(job)))
        except REFUSALS as error:
            return failure(error)

    @server.tool(annotations=CANCELS_A_JOB)
    async def cancel_job(job: str) -> CallToolResult:
        """Stop a running job started by this server. The output folder keeps whatever was
        written so far; start a new run in a new folder."""
        try:
            return ctx.job_envelope(await ctx.jobs.cancel(ctx.root.existing(job)))
        except REFUSALS as error:
            return failure(error)

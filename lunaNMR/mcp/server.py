# ABOUTME: The lunaNMR MCP server: typed tools over the headless CLI, confined to one data folder.
# ABOUTME: Run it as `lunanmr-mcp --root <dataset folder>`; it speaks MCP over stdio.

import sys

from mcp.server.mcpserver import MCPServer

from lunaNMR.mcp import (check, tools_dynamixs, tools_jobs, tools_kd, tools_peaks,
                         tools_project, tools_series)
from lunaNMR.mcp.command_line import build_parser
from lunaNMR.mcp.context import ToolContext
from lunaNMR.mcp.data_root import DataRoot, PathRefused
from lunaNMR.mcp.jobs import JobManager

WAIT_SECONDS = 20


def build_server(root, max_running_jobs=2, wait_seconds=WAIT_SECONDS):
    server = MCPServer(
        'lunaNMR',
        instructions='NMR peak fitting and analysis on the data folder this server was '
                     'started with. Every path is relative to that folder.')
    ctx = ToolContext(root, JobManager(max_running_jobs), wait_seconds)
    for tools in (tools_series, tools_kd, tools_peaks, tools_dynamixs, tools_project, tools_jobs):
        tools.register(server, ctx)
    return server


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.check:
        sys.exit(check.report(check.run_checks()))
    if args.root is None:
        parser.error('--root is required unless --check is given')
    try:
        root = DataRoot(args.root)
    except PathRefused as error:
        parser.error(str(error))
    build_server(root).run()


if __name__ == '__main__':
    main()

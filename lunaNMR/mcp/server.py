# ABOUTME: The lunaNMR MCP server: typed tools over the headless CLI, confined to one data folder.
# ABOUTME: Run it as `lunanmr-mcp --root <dataset folder>`; it speaks MCP over stdio.

import argparse

from mcp.server.mcpserver import MCPServer

from lunaNMR.mcp import cli_call
from lunaNMR.mcp.data_root import DataRoot, PathRefused


def _refused(error):
    return {'ok': False, 'error': {'type': 'PathRefused', 'message': str(error)}}


def build_server(root):
    server = MCPServer(
        'lunaNMR',
        instructions='NMR peak fitting and analysis on the data folder this server was '
                     'started with. Every path is relative to that folder.')

    @server.tool()
    async def diagnose(path: str = '.', quick: bool = False) -> dict:
        """Read-only pre-flight over a dataset folder: per-spectrum registration offset and
        capture rate, delays, peak-list hygiene, repeat-acquisition scale. Run it before
        any fit. Changes nothing."""
        try:
            folder = root.existing(path)
        except PathRefused as error:
            return _refused(error)
        return await cli_call.run_cli(cli_call.diagnose_args(folder, quick=quick))

    return server


def main(argv=None):
    parser = argparse.ArgumentParser(prog='lunanmr-mcp', description=__doc__)
    parser.add_argument('--root', required=True,
                        help='Dataset folder; the only place the server reads or writes')
    args = parser.parse_args(argv)
    try:
        root = DataRoot(args.root)
    except PathRefused as error:
        parser.error(str(error))
    build_server(root).run()


if __name__ == '__main__':
    main()

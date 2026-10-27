# ABOUTME: Builds `python -m lunaNMR` argument vectors and runs them in a subprocess for the MCP server.
# ABOUTME: The CLI's JSON envelope is the contract; a run that prints none becomes an error envelope.

import asyncio
import json
import sys

STDERR_TAIL_LINES = 40


def diagnose_args(path, quick=False):
    args = ['diagnose', str(path), '--format', 'json']
    if quick:
        args.append('--quick')
    return args


async def run_cli(args):
    """Run one CLI command and return its JSON envelope.

    A subprocess keeps the engines' prints, and any spawn worker's, off the server's
    stdout, which carries the MCP protocol.
    """
    process = await asyncio.create_subprocess_exec(
        sys.executable, '-m', 'lunaNMR', *args,
        stdout=asyncio.subprocess.PIPE, stderr=asyncio.subprocess.PIPE)
    stdout, stderr = await process.communicate()
    try:
        return json.loads(stdout.decode())
    except ValueError:
        tail = '\n'.join(stderr.decode().splitlines()[-STDERR_TAIL_LINES:])
        return {'ok': False,
                'error': {'type': 'NoJsonOutput',
                          'message': f'exit {process.returncode}, no JSON on stdout. {tail}'}}

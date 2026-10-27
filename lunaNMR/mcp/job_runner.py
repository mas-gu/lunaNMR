# ABOUTME: Runs one lunaNMR CLI command as a detached job and records how it ended.
# ABOUTME: Stdlib only: it must outlive the MCP server that started it and import nothing heavy.

import json
import os
import subprocess
import sys
import time
from pathlib import Path


def main(job_dir, cli_args):
    job_dir = Path(job_dir)
    with open(job_dir / 'stdout.json', 'wb') as stdout, open(job_dir / 'stderr.log', 'wb') as stderr:
        completed = subprocess.run([sys.executable, '-m', 'lunaNMR', *cli_args],
                                   stdout=stdout, stderr=stderr)
    result = {'exit_code': completed.returncode, 'finished': time.time()}
    temporary = job_dir / 'result.json.tmp'
    temporary.write_text(json.dumps(result))
    os.replace(temporary, job_dir / 'result.json')
    return completed.returncode


if __name__ == '__main__':
    sys.exit(main(sys.argv[1], sys.argv[2:]))

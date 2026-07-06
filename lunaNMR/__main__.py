# ABOUTME: Enables `python -m lunaNMR <subcommand>` by delegating to the CLI dispatcher.
# ABOUTME: Keeps process exit codes propagated from lunaNMR.cli.main.

import sys

from lunaNMR.cli import main

if __name__ == "__main__":
    sys.exit(main())

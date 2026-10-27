<!-- ABOUTME: Install and client setup for the lunaNMR MCP server, with its tools, guardrails and limits. -->
<!-- ABOUTME: Written for a user connecting an AI app (Claude Code, Claude Desktop) to their own lunaNMR install. -->

# lunaNMR MCP server

Lets an AI app run lunaNMR analyses on **your machine** through typed tools. Spectra never
leave the computer; only the short summaries a tool returns go to the model.

The server runs each analysis by calling the lunaNMR command line (`docs/CLI.md`) in a
subprocess, so results are identical to running the commands yourself.

## Install

The server needs Python 3.10 or newer (the MCP SDK's requirement); the rest of lunaNMR runs on 3.9.

```bash
pip install -e ".[mcp]"          # from the lunaNMR_v1o0 folder, inside your lunaNMR environment
lunanmr-mcp --check              # verifies the modules and a real CLI call
```

`--check` prints one line per problem (`FAIL numba: cannot import numba: ...`). Fix those
first; the server needs the same environment as the CLI (numba, nmrglue, scipy).

## Connect a client

The client must start **the Python of the environment lunaNMR is installed in**. Give that
Python's full path; a bare `python` usually finds a different one.

Claude Code:

```bash
claude mcp add lunanmr -- /path/to/env/bin/python -m lunaNMR.mcp.server --root /path/to/dataset
```

Claude Desktop (`claude_desktop_config.json`):

```json
{"mcpServers": {"lunanmr": {
  "command": "/path/to/env/bin/python",
  "args": ["-m", "lunaNMR.mcp.server", "--root", "/path/to/dataset"]}}}
```

`--root` is the one folder the server may read or write. Start it on the folder that holds
the spectra, peak lists and where results should go.

Tested in Claude Code (connect, every tool, refusals, cancel). Not yet tried in Claude Desktop; its
config above follows the documented format.

## Tools

| Tool | Does | Writes |
|---|---|---|
| `diagnose` | Pre-flight over a dataset: registration, capture, delays, peak lists | nothing |
| `plan_series` | Dry run of a series fit: checks inputs, reports the plan | nothing |
| `start_series` | Starts a series fit in the background; returns a job | new folder |
| `kd_survey` | Rates every residue of a titration and proposes a selection | new folder, plus a residue list and settings file **beside the input** |
| `kd_fit` | Starts a Kd fit in the background; returns a job | new folder |
| `job_status` | State of a job; on completion what the exit code hides | nothing |
| `cancel_job` | Stops a job started by this server | marks the job cancelled |
| `export_kd` | Figures and tables from a Kd fit JSON | new folder |

Typical order: `diagnose` -> `plan_series` -> `start_series` -> `job_status` (until done) ->
`kd_survey` -> `kd_fit` -> `job_status` -> `export_kd`.

A finished series job reports how many points were pinned on a fit bound or fell below
3 sigma, and which residues; a finished Kd job reports each shared Kd with its `reliable`
flag and how many residues are in the shared fit. Read those before trusting a number.

## Guardrails

- **One data root.** Every path is resolved (symlinks included) and refused if it leaves
  `--root`.
- **No overwriting.** An output folder that already exists is refused; choose a new name.
- **Links stay inside.** A symbolic link anywhere in a folder a tool will scan, or beside a Kd
  input, that points outside `--root` blocks that folder: scanning follows links on its own, so
  the folder would otherwise read outside the root. Remove the link or copy the data in.
- **Residue names are names.** `kd_fit(residues=[...])` accepts letters, digits and `_ . + -`
  (up to 64 characters, not starting with `-`). The server writes them to a selection file it
  owns; they are never passed as a file path. To use a survey's selection, pass `residues_file`.
- **Bounded output.** Text from the dataset (log lines, residue names, reasons) is cut to a fixed
  length, and a result over 200,000 characters is replaced by an error.
- **Failures are failures.** A refused path or a failed run comes back with the protocol's
  error flag and the reason, not as a success.

## Limits

- A job's state is kept in `<out>/mcp_job/`, so a restarted server can still read it. Only
  the server that started a job can cancel it.
- Cancelling stops the work and keeps whatever was written so far.
- Cancelling uses POSIX process groups (macOS, Linux). Windows is not handled.
- There is no time limit on a job.
- Relaxation (T1/T2/model-free) and 1D integration are not exposed yet; use the command
  line for those.

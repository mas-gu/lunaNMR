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

Every analysis the command line offers is a tool, except the two listed under "Not exposed".

**Check a dataset**

| Tool | Does | Writes |
|---|---|---|
| `diagnose` | Pre-flight over a dataset: registration, capture, delays, peak lists | nothing |
| `plan_series` | Dry run of a series fit: checks inputs, reports the plan | nothing |

**Peak lists and series**

| Tool | Does | Writes |
|---|---|---|
| `start_series` | Fits a series (time, titration, temperature) | new folder |
| `peaks_shift` | Moves a peak list into register, by hand or measured against a spectrum | one new file |
| `integrate_1d` | Measures peaks across a 1D series | new folder |

**Binding (Kd)**

| Tool | Does | Writes |
|---|---|---|
| `kd_survey` | Rates every residue of a titration and proposes a selection | new folder, plus a residue list and settings file **beside the input** |
| `kd_fit` | Fits the Kd (CSP and intensity) | new folder |
| `export_kd` | Figures and tables from a Kd fit | new folder |

**Relaxation and dynamics**

| Tool | Does | Writes |
|---|---|---|
| `dynamixs_t1t2` | Mono-exponential T1 or T2 | new folder |
| `dynamixs_methyl_t2` | Bi-exponential methyl T2 | new folder |
| `dynamixs_hetnoe` | Heteronuclear NOE from the two planes | new folder |
| `dynamixs_t1rho` | T1rho converted to R2 with the tilt-angle correction | new folder |
| `dynamixs_density` | Reduced spectral density mapping | new folder |
| `dynamixs_modelfree` | Model-free (Lipari-Szabo), single or dual field | new folder |

**Project bundles**

| Tool | Does | Writes |
|---|---|---|
| `project_build` | Assembles a `.lunaNMR` bundle from a run | new folder |
| `project_inventory` | Lists what a bundle holds | nothing |
| `project_export` | Copies a bundle to a plain folder | new folder |

**Jobs**

| Tool | Does | Writes |
|---|---|---|
| `job_status` | State of a job; on completion what the exit code hides | nothing |
| `cancel_job` | Stops a job started by this server | marks the job cancelled |

**How a run answers.** A tool that runs an analysis starts it as a job and waits about 20
seconds. A quick analysis (hetNOE, methyl T2, density, a Kd fit) comes back finished in the same
call; a slow one (a series, a dual-field model-free, about 40 s) comes back `running` with a
job to poll with `job_status`. Every one of them takes `dry_run`, which checks the inputs and
reports what would run, without running or writing anything.

Typical order: `diagnose` -> `plan_series` -> `start_series` -> `kd_survey` -> `kd_fit` ->
`export_kd`; or, for relaxation, `start_series` per experiment -> `dynamixs_hetnoe` ->
`dynamixs_modelfree`.

**What a finished job tells you beyond success:**
- series: points pinned on a fit bound or below 3 sigma, and which residues;
- Kd: each shared Kd with its `reliable` flag, how many residues are in the shared fit, and the
  residues that did not fit, with the reason;
- T1/T2: the residues that fail the usual gate (fitted, finite T, error below the value);
- hetNOE: the median and how many residues fall outside the possible range;
- model-free: residues kept at each stage and any fit that hit the iteration cap;
- methyl T2: how many fits were unidentifiable (collapsed to a single exponential);
- every job: the files it wrote.

Read these before trusting a number. In the T1/T2 result, `mean_t2` is the mean time constant
of whichever experiment was fitted, T1 included.

## Not exposed

- `project remove` deletes bundle contents irreversibly. The Project Browser in the GUI does it
  with confirmation.
- `batch` is a pass-through to another parser with no flags of its own.
- `project export --force` would overwrite an existing folder, and `dynamixs t1t2 --no-json`
  would drop the file the server reads results from.

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
- **Names are names.** A prefix or field name becomes part of an output filename, so it must be
  1-64 letters, digits, `_`, `.` or `-` and start with a letter or digit; `../x` is refused.
- **Failures are failures.** A refused path or a failed run comes back with the protocol's
  error flag and the reason, not as a success.

## Limits

- A job's state is kept in `<out>/mcp_job/`, so a restarted server can still read it. Only
  the server that started a job can cancel it.
- Cancelling stops the work and keeps whatever was written so far.
- Cancelling uses POSIX process groups (macOS, Linux). Windows is not handled.
- There is no time limit on a job.
- `project_build` needs `spectrum` and `peaks` named explicitly, inside the data folder, so a
  bundle never reads paths recorded in a run's own metadata.
- No real T1rho series is in the test data; that tool is tested for its inputs and outputs with a
  stand-in matrix, not for its physics.

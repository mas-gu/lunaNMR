# AGENTS.md — driving lunaNMR from an agent

LunaNMR is an NMR peak-analysis and integration suite. Every analysis the GUI performs is
also reachable headlessly, with no display, through one CLI.

```bash
cd lunaNMR_v1o0
python -m lunaNMR --help          # `lunanmr` also works after `pip install -e .`
```

## Read these before running anything

| Document | What it is |
|---|---|
| [`docs/CLI_AGENT.md`](docs/CLI_AGENT.md) | **Start here.** The machine contract: subcommand table, `--format json` output shapes, exit codes, input file shapes, and the gotchas that silently corrupt results rather than raising. |
| [`docs/CLI_AGENTS_DEEP/`](docs/CLI_AGENTS_DEEP/) | Long-form runbooks — the phase structure, physical QC bands, and worked end-to-end flows for relaxation (`RELAXATION_*.md`) and affinity/titration (`AFFINITY_*.md`). |
| [`docs/CLI.md`](docs/CLI.md) | The human reference: every subcommand and flag. `--help` is authoritative when they disagree. |

The gotchas are the point. Most failure modes here are **not** exceptions — a wrong peak
list produces zero heights and R²≈0.2 fits that read as noisy data; a mislabelled
concentration unit produces a confident Kd wrong by a constant factor; an unparseable delay
produces a complete relaxation table built on times never measured. Each is documented, and
each exits 0.

## Before a real run

```bash
python -m lunaNMR diagnose <dataset_root> --format json
```

Read-only, seconds, writes nothing. It checks peak-list registration, capture rate, delay
parsing, hetNOE plane identity and the cross-experiment residue sets — the questions every
later step depends on.

## Conventions

- `--format json` puts one JSON summary object on stdout; all engine output goes to stderr.
- `--dry-run` validates inputs and prints the plan without running.
- Exit codes are listed in `docs/CLI.md` under Global conventions.

## Repository layout

- `lunaNMR/` — the package (`cli.py` is the CLI entry point; `core/` holds the fitting engines).
- `modules/` — optional analysis modules (dynamiXs relaxation, Kd titration, 1D integration).
- `docs/` — everything above.
- `CLAUDE.md` (one level up) — engineering conventions for contributors, not usage.

## Tests

```bash
make test        # all three test roots, headless
```

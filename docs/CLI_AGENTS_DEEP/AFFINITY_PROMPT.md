# Task: full binding titration analysis, spectra → Kd

Run the complete analysis on this dataset end to end, and **record exact wall-clock timing
for every stage**. Work autonomously — do not stop to ask except at the two points below
and on a genuine blocker (see *When to stop*).

---

## Setup

- **Dataset root:** given to you by the user. This prompt lives in `docs/CLI_AGENTS_DEEP/`,
  not in the dataset — do not infer the root from this file's location. Each experiment
  is its own subfolder with its spectra and its own assignment file.
- Read `docs/CLI_AGENT.md` in the lunaNMR package first — it is the flag and output-path
  reference, and reading it once is cheaper than discovering flags by trial.
- **Python:** the environment carrying nmrglue, scipy, PySide6 and numba. Confirm before
  starting — the system interpreter will not have them.
- **lunaNMR CLI:** `python -m lunaNMR ...`, run from the lunaNMR package root.
- **Output:** create a **new** folder per experiment,
  `<YYYYMMDD_HHMM>_kd_fit_<experiment folder name>/`, inside that experiment's folder.
  **Never re-use or overwrite a previous run's folder** — each run is a record of what one
  set of settings produced, and comparing runs is how you tell a real change from a
  reformatted one. Do not modify the raw spectra or the assignment files.

**Read `AFFINITY_PLAYBOOK.md` (beside this file) first** — it has the commands, the checks and the traps.
Establish everything from the spectra; assume no prior result.

---

## Timing — the point of this run

Record epoch seconds at the start and end of every stage. Write `run_<stamp>/TIMING.md`
as you go — append after each stage, do not buffer to the end; a crash must not lose the
timings.

```
| stage | start | end | elapsed_s | compute_s | notes |
```

Split **compute time** (subprocesses running) from **agent time** (your reasoning,
inspection, decisions) where you can — the second is what matters for comparing agents.
Finish with total wall-clock, total compute, number of CLI invocations, number of spectra
integrated, and number of residues surviving to a reported Kd.

---

## Speed — use the whole machine

Measure the core count first (`sysctl -n hw.ncpu`) and size the work to it.

- Pass `--parallel` to every `series` run.
- Separate experiments are independent: integrate them concurrently, don't queue them.
- Overlap compute with reasoning — launch integrations in the background and prepare the
  downstream steps while they run.
- Don't oversubscribe: each `--parallel` job tries to use the whole machine.
- Where stages overlap, say so in TIMING.md — per-stage times are contended and their sum
  will exceed true wall-clock.

---

## Stages

**1. Diagnose** — `python -m lunaNMR diagnose <root> --mode titration`. Read-only. Record
the report verbatim, and add any playbook Phase A check it does not cover — in particular
the per-residue reference-point sanity check, which `diagnose` does not perform.
**Stop and show the user.**

**2. Integrate** — one `series` run per experiment with its own peak list,
`--mode titration --peak-source cascade`. Cascade because peaks move; reference-mode
tracking loses movers at high titrant.

**3. Survey** — `kd --survey`. Produces the vs-sequence figures, the evidence table and a
proposed selection with reasons. **Stop and show the user** the figure and the proposed
exclusions. This is the decision the whole workflow exists for; do not make it yourself.

**4. Fit** — `kd --residues <the user's selection>`, both observables.

**5. Export** — `export kd`, all figure kinds.

**6. QC** — playbook Phase D, in order: `quality_warnings`, then `reliable`, then pool
reconciliation, then the resolvable-window sanity band. **Report failures; do not tune
them away.** A failure here does not stop the run.

**7. Report** — `run_<stamp>/RESULTS.md`: per observable, the shared Kd with its error and
`reliable` flag, pool size, the residues in the pool, every gate's exclusion count with
reasons, and every QC verdict. State plainly whether the dataset supports a quotable Kd.

---

## Pre-authorised decisions

So you don't block. **Log each one you apply, with the measured number that justified it.**

1. **Concentration units** — if titration points parse to values small relative to `--p0`,
   they are equivalents. Pass `--conc-units equivalents` and log the resolved
   concentrations.
2. **Peak source** — `cascade` for titration series.
3. **Observables** — run both CSP and intensity unless the input carries no positions.
4. **Accepting the survey verbatim** — only if the user has explicitly said so. Otherwise
   the selection is theirs.

Anything **not** on this list that would change the data: don't do it. Report it instead.

**Read-only analysis is always in bounds, and is wanted.** Extra cross-checks, fault
isolation, re-deriving a number a different way — none of it changes the data, so none of
it needs permission. The unit-invariance check in the playbook is a good example: it costs
one extra fit and independently tells you which per-residue Kd values were never
determined. The most useful findings usually come from a diagnostic nobody asked for.

---

## When to stop and ask

Two mandatory stops: after **diagnose**, and after **survey**. Both exist because the
decision belongs to the user, not because you are unsure.

Otherwise run to completion. QC failures are **results, not blockers** — record the
measured number, name the threshold it missed, carry it into RESULTS.md, and keep going.

Stop early only if you genuinely cannot produce output:

- A peak list is mis-registered against its own spectra (shift > 50% of the search window).
- Capture at the **reference** point is below 50%.
- A CLI command fails twice for the same reason.

---

## Deliverables

```
run_<stamp>/
├── TIMING.md      per-stage wall-clock + compute, appended live
├── RESULTS.md     Kd per observable, pool composition, QC verdicts
├── LOG.md         what you ran, what you decided, what the user chose
├── scripts/       every step re-runnable
├── series/        integration output
└── kd/            survey, fits, figures  (PDFs at top, CSV/JSON under data/)
```

## Ground rules

- Measure, don't infer from filenames.
- Each experiment gets **its own** peak list, and **its own timestamped output folder**.
  A rerun never overwrites the run before it.
- Never report a Kd without naming its observable and its `reliable` flag.
- A residue absent from a pool has a stated reason — quote it rather than the count.
- Report what actually happened, including failures and anything you skipped.
- Don't claim a stage succeeded without showing the number that proves it.

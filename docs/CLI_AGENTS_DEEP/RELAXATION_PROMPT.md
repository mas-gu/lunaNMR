# Task: full ¹⁵N relaxation analysis, spectra → model-free

Run the complete analysis on this dataset end to end, and **record exact wall-clock timing
for every stage**. Work autonomously — do not stop to ask unless you hit a genuine blocker
(see *When to stop*).

---

## Setup

- **Dataset root:** given to you by the user. This prompt lives in `docs/CLI_AGENTS_DEEP/`,
  not in the dataset — do not infer the root from this file's location.
- Read `docs/CLI_AGENT.md` in the lunaNMR package first — it is the flag and output-path
  reference, and reading it once is cheaper than discovering flags by trial.
- **Python:** the environment carrying nmrglue, scipy, PySide6 and numba. Confirm before
  starting — the system interpreter will not have them.
- **lunaNMR CLI:** `python -m lunaNMR ...`, run from the lunaNMR package root.

- **Output:** create a **new** folder `run_<YYYYMMDD_HHMM>/` in the dataset root and put
  everything there. Do not modify the raw spectra or the peak lists.

**Read `RELAXATION_PLAYBOOK.md` (beside this file) first** — it has the commands, the checks and the
failure modes. There is no pre-existing inventory or prior result for this dataset:
establish everything yourself from the spectra.

The playbook's Phase A ends with "stop and show the user". **For this run that is
overridden** — diagnose, record the report in your log, then keep going under the
pre-authorised decisions below.

---

## Timing — the point of this run

Record epoch seconds at the start and end of every stage. Write
`run_<stamp>/TIMING.md` as you go (append after each stage, don't buffer to the end — a
crash must not lose the timings).

For each stage record: name, start, end, elapsed seconds, and what was produced.
Also split **compute time** (subprocesses running) from **agent time** (your reasoning,
inspection, decisions) where you can — the second is what matters for comparing agents.

```
| stage | start | end | elapsed_s | compute_s | notes |
```

Finish with a total, and a one-line summary: total wall-clock, total compute, number of
CLI invocations, number of spectra fitted.

---

## Speed — use the whole machine

Measure the core count first (`sysctl -n hw.ncpu` on macOS, `nproc` on Linux) and size the
work to it. Wall-clock is being measured, so idle cores are a cost.

- Pass `--parallel` to every `series` run. `modelfree` and `density` parallelise internally.
- Independent work that should run concurrently, not one at a time: the six `series`
  integrations, the T1/T2 fits, and the model-free variants.
- Don't oversubscribe — each `--parallel` job tries to use the whole machine, so a few
  concurrent jobs is the useful range, not all six at once.
- Overlap compute with your own reasoning: launch the long integrations in the background
  and prepare the downstream steps while they run.
- Where stages overlap, say so in TIMING.md — per-stage elapsed times are contended and
  their sum will exceed the true wall-clock.

---

## Stages

**1. Diagnose** — run `python -m lunaNMR diagnose <dataset_root>` (read-only; it fits and
writes nothing). Record its report verbatim in the log, and add any playbook Phase A check it
does not cover.

**2. hetNOE** — integrate both fields, compute I_sat/I_unsat per residue, gate to usable
residues, write the model-free sat/unsat inputs.

**3. T1** — integrate and fit both fields.

**4. T2** — integrate and fit both fields.

**5. Physical QC** — R1 ratio, R2 ratio, τc per field. Report pass/fail against the
playbook thresholds. **Report failures; do not tune them away.** A failure here does not
stop the run.

**6. Model-free** — dual-field, plus each field alone as a cross-check.

**7. Report** — `run_<stamp>/RESULTS.md`: τc, S² median/IQR, τe, flexible residues,
exchange residues, residue counts surviving each stage, and every QC verdict.

---

## Pre-authorised decisions

So you don't block, these are approved in advance. **Log each one you apply, with the
measured number that justified it, and keep the uncorrected version alongside.**

1. **hetNOE plane identity** — if planes are named `001`/`002`, determine saturated by
   intensity ratio over high-S/N peaks (saturated is lower). Symlink to explicit names.
2. **hetNOE usable gate** — both planes non-zero, both fits Excellent/Good, ratio in
   −0.2…1.05.
3. **Split T2 sub-series** — if the shared-delay ratio deviates from 1 by >5%,
   cross-normalise the low sub-series onto the other. Keep the raw fit too.
4. **hetNOE error floor** — 0.044 relative per plane (→ hetNOE_err ≈ 0.05).
5. **Fit filtering** — keep fits with finite `T`, `T > 0`, `T_err < T`.

Anything **not** on this list that would change the data: don't do it. Report it instead.

**Read-only analysis is always in bounds, and is wanted.** Extra measurements, cross-checks,
fault isolation, re-deriving a number a different way — none of that changes the data, so
none of it needs permission. The restriction above is on *modifying* data, not on thinking.
The most useful findings usually come from a diagnostic nobody asked for.

---

## When to stop and ask

Run to completion. QC failures are **results, not blockers** — record the measured number,
name the threshold it missed, carry it into RESULTS.md, and keep going. Deciding whether a
number is trustworthy is the user's call, not yours.

Stop only if you genuinely cannot produce output at all:

- A peak list is mis-registered against its own spectra (shift > 50% of the 0.070 ppm window).
- Any stage yields fewer than 30 usable residues.
- A CLI command fails twice for the same reason.

---

## Deliverables

```
run_<stamp>/
├── TIMING.md      per-stage wall-clock + compute, appended live
├── RESULTS.md     final numbers and QC verdicts
├── LOG.md         what you ran, what you decided, what you corrected
├── scripts/       every step re-runnable
├── series/        integration output
├── fits/          T1/T2 fits
├── modelfree/     dual + per-field
└── results/       per-residue CSVs and QC figures
```

## Ground rules

- Measure, don't infer from filenames.
- Each experiment gets **its own** peak list — hetNOE lists are offset from T1/T2 lists and
  are not interchangeable.
- Report what actually happened, including failures and anything you skipped.
- Don't claim a stage succeeded without showing the number that proves it.

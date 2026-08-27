# Runbook: titration → Kd (lunaNMR)

For an agent running a binding titration on any dataset of this shape.

**Four phases, and they do not mix.**

| | A — DIAGNOSE | B — SURVEY | C — FIT | D — QC |
|---|---|---|---|---|
| changes anything? | **no** | no Kd is produced | yes | **no** |
| output | a report | evidence + a proposed selection | fits and figures | verdicts |
| ends with | **STOP. Show the user.** | **STOP. The user chooses.** | results | pass/fail |

**The rule: fitting every residue to a Kd is meaningless.** Most residues do not bind.
The agent proposes a selection with reasons; the human disposes. A Kd produced without
that step is a number nobody chose.

---

## Phase A — diagnose

```bash
python -m lunaNMR diagnose <dataset_root> --mode titration
```

Read-only. Reports per-spectrum registration offset, capture rate, median S/N, parsed
titration points, and peak-list hygiene.

| check | pass | if it fails |
|---|---|---|
| **peak-list registration** | shift ≈ 0.000, 0.000 | wrong list or mis-referenced. Symptom is Height 0 / R² ≈ 0.2 shoulder fits that read as "noisy data" |
| **capture rate at the reference point** | >85% at zero titrant | peaks outside the search window. Judge **only** at the reference — a low rate further along a titration is exchange broadening, i.e. the experiment working |
| **assignment hygiene** | no whitespace, no `dummy_*`, no duplicates | whitespace breaks exact-string merging silently |
| **titration points parse** | one value per spectrum, ascending | the trailing token needs a **decimal separator** (`_0o5`, `_2.4`), not merely digits — `sample_2.ft` and `ref_100.ft` do **not** parse and are **silently dropped**, not errored. Name points `_0o0`, `_0o5`, `_2o0`. |
| **reference peak sanity** | every residue's I(0) comparable to its own series max | a reference orders of magnitude low rescales that residue entirely, and survives a `quality='Good'` label |

**Capture collapse at high titrant is expected, not a fault.** The most strongly perturbed
residues broaden past detection first. Record the profile; do not treat it as a windowing
problem. But note the consequence for Phase D: the plateau that pins Kd is exactly what
gets lost.

**Then stop.** Present the report. Propose fixes and their magnitude; do not apply them.

---

## Phase B — survey, then let the user choose

```bash
python -m lunaNMR kd --survey --input <tidy.csv> --out <RESULTS> --p0 <conc> \
    [--conc-units equivalents] --format json
```

Produces `<out>/<prefix>_<obs>_vs_sequence.pdf` (prefix first, one per observable),
`<out>/data/<prefix>_survey.csv`, and an editable `<prefix>_residues.txt` **beside the
input CSV** — not in `--out`, so a different output path cannot lose it.

Selection file: one residue per line, `#` comments a line out, deleting the `#` puts it
back. Trailing comments are **evidence, not exclusions**. Re-running `--survey` merges —
human decisions survive, evidence refreshes, disagreements are marked `[kept]`/`[dropped]`
and never re-suggested. So iterating is cheap and safe.

Present the vs-sequence figure and the proposed exclusions **with their reasons**. Do not
apply your own judgement in place of the user's.

---

## Phase C — fit only what was chosen

```bash
python -m lunaNMR kd --input <tidy.csv> --out <RESULTS> --p0 <conc> \
    --residues <selection.txt> --observable csp,intensity [--conc-units equivalents]

python -m lunaNMR export kd --json <RESULTS>/data/<prefix>_kd_fit_data.json --out <RESULTS>
```

`--dry-run` first. Output layout: `*.pdf` and `*_kd_results.txt` at the top of `--out`,
every CSV and JSON under `<out>/data/`.

**Never write two runs into the same folder.** `--out` must be
`<YYYYMMDD_HHMM>_kd_fit_<dataset folder name>`, created fresh inside the dataset folder.
Overwriting in place destroys the only copy of what the previous settings produced, so a
result cannot be compared against the run before it — and a run that fails midway leaves a
folder that is half old and half new, with nothing to say which file is which. Keep every
run; they are small, and the figures are the record of what a given set of thresholds
actually did.

**The two observables are different models and are not interchangeable.**

| | CSP | intensity |
|---|---|---|
| model | full 1:1 quadratic isotherm, uses [P]₀ | phenomenological decay + fitted plateau |
| Kd is | thermodynamic | **apparent** — ignores [P]₀ |
| needs | 4-decimal fitted positions | works on existing data |

A CSP fit on a series whose positions were rounded to 0.1 ppm is meaningless — re-run
`series` first. Intensity is unaffected.

---

## Phase D — QC, before believing anything

**1. Read `metadata.quality_warnings` first.** When the *typical* residue fits a Kd the
titration cannot resolve, the tool says so in plain language and refuses to exclude
outliers. That verdict means the whole dataset is unusable, not that a few residues are
bad. No Kd from such a run should be reported at any confidence.

**2. Check `global.<obs>.reliable`.** False means the shared Kd is pinned at a bound or its
relative error exceeds 30%. A pinned Kd is the optimizer hitting a wall, not a measurement.

**3. Reconcile the pool.** `metadata.csp_pool_excluded` names the gate that removed every
residue — failed fit, no CSP at the last point, below the significance threshold, R² too
low, outside the resolvable window, statistical outlier. Pool + excluded must equal every
fitted residue. If a pool is surprisingly small, the answer is in that map, not in the
figures. **Placeholders are not in it**: `dummy_*` rows are dropped before fitting, so they
are never *fitted* residues — their count is `metadata.n_excluded_dummy`, a separate key.

**4. Sanity band.** The resolvable window is roughly one decade either side of the
titration's own concentrations. A Kd outside it was not measured, however good its R².

---

## The traps, and why each one bites

| trap | why it is invisible |
|---|---|
| **Equivalents read as concentrations** | filenames commonly encode equivalents (`_3o0`). Read as absolute they give a confident Kd wrong by the factor `--p0`, with `success=True` and no warning. Pass `--conc-units equivalents` |
| **A broken reference point** | I(0) far below the residue's own series max rescales *every* point for that residue. It is not a bad point, it is a bad denominator — and it passes a `quality='Good'` label |
| **Failed points dropped residue-wise** | a residue that failed at one titration point still has good data at the others. Dropping the whole residue silently discards it |
| **σ significance is not scale-free** | at an early titration point where nothing has happened, the spread is tiny and most residues clear it. The criterion is only meaningful at the **highest** titrant point |
| **Outlier tests based on spread** | the bad fits inflating the spread are the ones the test should catch, so the gate widens to admit them. Apply the physical window **first**, then statistics on the survivors |
| **R² on a truncated titration** | R² measures fit to the sampled points and cannot constrain a plateau extrapolated beyond the last one. A residue can fit an impossible Δδ_max at R² > 0.95 |
| **`Kd_err` / `r_squared` of `null`** | a degenerate fit has a singular covariance. An unbounded error is legitimate — the fit did not fail. Any consumer must handle `None` |
| **"unmeasured" ≠ "did not shift"** | residues that vanish at high titrant are *enriched* for the binding surface, since the most perturbed broaden first. A significant-residue list is a **floor** on the interface, not a description of it |

**One free diagnostic worth more than any threshold:** refit with all concentrations and
[P]₀ scaled by the same factor. A well-posed Kd scales exactly; an unconstrained one does
not, because the optimizer's path depends on scale. No ground truth, no threshold, nothing
to tune — anything that fails to reproduce was never measured.

---

## Speed

- `--dry-run` validates inputs and prints the plan without running.
- Settings persist per dataset and are read back automatically — concentrations,
  `conc_units`, `--p0`, `--alpha`, `--intensity-from` and every threshold, recorded as you
  typed them. Don't re-enter them, and don't re-derive what a prior run already recorded.
  A survey writes `<input_dir>/<prefix>_kd_params.json`; a fit writes
  `<out>/data/<prefix>_kd_params.json`, which the next run finds when `--out` is the
  conventional `<input_dir>/kd_analysis`. The run summary's `params_source` names the file
  that was read, or `null` if none was — check it rather than assuming.
- Independent datasets are independent processes: run them concurrently, don't queue them.
- The survey is cheap (per-residue fits only, no global). Iterate on the selection there
  rather than re-running full fits.
- Read `docs/CLI_AGENT.md` once at the start instead of discovering flags by trial.

## Rules

- Measure, never infer from filenames.
- Diagnose → report → **stop**. Survey → propose → **stop**. The user decides both times.
- Report QC failures rather than tuning them away. A number that missed a threshold is a
  result; deciding whether to trust it is the user's call.
- Never quote a Kd without stating which observable produced it and whether `reliable` was
  true.

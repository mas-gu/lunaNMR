# Runbook: ¹⁵N relaxation → model-free (lunaNMR)

For an agent running T1/T2/hetNOE on any dataset of this shape.

**Two phases, and they do not mix.**

| | Phase A — DIAGNOSE | Phase B — RUN |
|---|---|---|
| changes anything? | **no** — read-only | yes |
| output | a report | fits and results |
| ends with | **STOP. Show the user.** | results + QC |

**The rule: diagnose, report, stop. Never fix what you find on your own.**
Corrections (rescaling, exclusions, error floors, shifted peak lists) change the science.
Raise them, show the numbers, let the user decide. Then run.

---

## Phase A — diagnose

```bash
python -m lunaNMR diagnose <dataset_root> [--quick]
```

Read-only. Fits nothing, writes nothing. For **every** spectrum in every experiment folder it
reports the registration offset, capture rate and median S/N; per folder it reports the parsed
delays and repeats, peak-list hygiene, hetNOE saturated/unsaturated identity and the scale of
any repeat acquisitions; then it compares the assignment sets across experiments and prints
`FAIL` / `WARN` lines. Delays and peak lists are parsed with the same code the pipeline uses,
so the report cannot certify something the run will not reproduce.

What it checks and why each matters:

| check | pass | if it fails |
|---|---|---|
| **peak-list registration** | shift = 0.000, 0.000 | wrong list, or mis-referenced. hetNOE usually needs its own list — typically ~0.07 ppm ¹H / ~0.7 ppm ¹⁵N off the T1/T2 frame, i.e. 10× the 0.070 ppm search window. Symptom: Height 0 or R²≈0.2 shoulder fits, which read as "noisy data" |
| **capture rate** | >85% at the **reference point** (shortest delay) | peaks outside the search window. Only meaningful before decay — further down the series a low rate is the experiment working. Judge by this, never by "zeros in the intensity matrix" — that metric can move the wrong way |
| **assignment hygiene** | no leading whitespace, no `dummy_*`, no duplicates | whitespace breaks exact-string merging (silently → 0 common residues) |
| **residue sets across fields** | identical | merge is an intersection; the smallest set is the ceiling |
| **hetNOE plane identity** | sat/unsat ≈ 0.75–0.85 | files named `001`/`002` don't say which is saturated. Decided by intensity ratio over high-S/N peaks; saturated is lower, never >1 |
| **split T2 sub-series** | shared-delay ratio ≈ 0.99 | ~0.85–0.90 is instrumental (scans / receiver gain). The low sub-series usually supplies the whole tail, and a free baseline compensates by stretching T2 |

**Then stop.** Present the report. If there are FAILs, propose the fix and its magnitude —
do not apply it.

---

## Phase B — run

Only after the user has seen Phase A.

```bash
# 1. integrate — one run per experiment, each with ITS OWN peak list
python -m lunaNMR series --spectra <dir> --peaks <list> --out <out> \
    --mode time --peak-source reference --parallel --format json

# 2. fit T1 / T2
python -m lunaNMR dynamixs t1t2 --input <out>/peak_intensity_matrix.csv \
    --out <fits> --exp {T1,T2} --field-freq <MHz> --time-units ms --format json

# 3. model-free  (--sat/--unsat are headerless: residue,intensity[,error])
python -m lunaNMR dynamixs modelfree [--dual] \
    --f1-t1 ... --f1-t2 ... --f1-noe-sat ... --f1-noe-unsat ... --field1-freq 600 \
    [--f2-... --field2-freq 700] --out <mf> --prefix <name> --format json
```

`--dry-run` first. A matrix written by `series` is already normalised to ms (`_2.4s` becomes
`2400.0`), and every relaxation subcommand now **reads that from `series_metadata.json`**
rather than asking you: the units flags are for tables `series` did not write, which have no
sidecar. A flag that contradicts the sidecar is refused rather than believed — a T1 in s
beside a T2 in ms puts R1 out by 1000×, and neither value looks wrong afterwards.

**`--time-units` does not mean the same thing everywhere.** On `t1t2`/`methyl-t2` it only
*labels* the output; on `t1rho` and via `--f{1,2}-t{1,2}-units` on `modelfree` it **rescales**.
The defaults differ to match (`s` vs `ms`) and are deliberately not unified, since changing
either moves numbers that have already been published. The full table is in `CLI_AGENT.md`.

Notes: symlink hetNOE planes to explicit `*_saturated.ft` / `*_unsaturated.ft` names before
integrating, and keep NMRPipe intermediates (`test.ft`) out of the folder glob.
Filter fits on finite `T`, `T > 0`, `T_err < T`.

To act on one sub-series alone — cross-normalising a re-acquired tail, say — read the
column → spectrum map from `series_metadata.json`. Matrix columns are delays, not filenames,
and which of two files at the same delay owns the bare label rather than the `_2` one follows
a case-sensitive sort of the basename, so it cannot be inferred from acquisition order.

---

## Phase C — physical QC, before believing anything

R1 = 1000/T1, R2 = 1000/T2 (ms), per residue, common set.

| check | expected | if it fails |
|---|---|---|
| R1(high)/R1(low) | **field-pair dependent — see below** | peak identity / referencing / scaling broken |
| R2(high)/R2(low) | **field-pair dependent — see below** | R2 gains a CSA term ∝ B₀²; failure ⇒ T2 problem or ΔT between datasets |
| τc from R2/R1, each field | agree within ~5% | dual-field assumes ONE global τc |

**The field ratios are not one band.** They follow from the dipolar + CSA rate expressions
and change a lot with the field pair, and only weakly with τc. Computed for a rigid ¹⁵N–H
(rNH 1.02 Å, CSA −160 ppm, S² = 1) over τc = 5–13 ns:

| field pair | R1(high)/R1(low) | R2(high)/R2(low) |
|---|---|---|
| 600 / 700 | **0.81–0.84** | **1.05–1.08** |
| 600 / 800 | **0.68–0.73** | **1.12–1.17** |

Applying a 600/700 band to a 600/800 dataset fails a perfectly good one: 0.71 is correct
there and would look like a gross failure against "0.80–0.90". For any other pair, compute it
rather than interpolating.

S² largely cancels in a field ratio, so these hold for a flexible protein too. Rex does not
cancel — it scales as B₀², so a measured R2 ratio *above* the band suggests exchange rather
than error, while one *below* it points at T2 or at the two datasets not sharing a τc.

τc = (1/4πν_N)·√(6·(R2/R1) − 7). Sanity: hetNOE 0.75–0.85, S² 0.8–0.9 rigid folded.

**R1 agreeing while R2 disagrees isolates the fault to T2** — the single most useful
diagnostic in the workflow. Report such a failure; do not tune it away.

---

## Model-free specifics

- **hetNOE error floor.** Noise-propagated errors understate reality ~5× (a single
  un-repeated sat/unsat pair has no reproducibility measure). ~0.044 relative per plane
  → hetNOE_err ≈ 0.05. Without it hetNOE is over-weighted and S² errors are fake-tight.
  *This is a correction — propose it, don't assume it.*
- **Always run each field alone alongside `--dual`.** Cheap, and it exposes the next point.
- **Rex inflation = forced global τc.** Dual giving Rex>1 on >50% of residues while
  single-field gives 20–25% (χ² ~25× higher) means the fit is absorbing a τc mismatch, not
  reporting chemistry. Trust single-field for per-residue Rex.
- **Decompose S² variance before interpreting per residue:**
  `real_spread = sqrt(observed_sd² − mean_err²)`. If noise is ~half the variance, report
  "uniform S² ≈ X with a few outliers", not a residue-by-residue profile. Cross-field S²
  correlation near r ≈ 0.2 with matching medians is the signature.

---

## Rules

- Measure, never infer from filenames.
- Diagnose → report → **stop**. The user decides on every correction.
- Log every intervention, especially anything that modifies data. Keep the uncorrected
  version alongside.
- Report physical-check failures rather than tuning them away.

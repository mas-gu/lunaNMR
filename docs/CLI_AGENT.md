<!-- ABOUTME: Agent-facing quick reference for driving the headless lunaNMR CLI programmatically. -->
<!-- ABOUTME: Condensed from docs/CLI.md + cli.py for an LLM/coding agent, not a human end-user. -->

# lunaNMR CLI — Agent Reference

Scannable condensation of `docs/CLI.md` + `cli.py` for an agent. Run from `lunaNMR_v1o0/`.

## Invocation
```bash
python -m lunaNMR <subcommand> [flags]   # == lunanmr <subcommand> [flags]
```
No Qt/display needed (`__init__` is lazy, matplotlib forced to `Agg`). Every analysis subcommand accepts `--format {text,json}` and `--dry-run`.

## Subcommands
| Command | Does | Required inputs | Key outputs |
|---|---|---|---|
| `diagnose` | Read-only pre-flight over a dataset: registration, capture, delays, peak lists, cross-experiment residue set | `<root>` | stdout report / JSON |
| `series` | Fit a multi-spectrum series/titration | `--spectra <dir/glob>` `--peaks` `--out` | `peak_intensity_matrix.csv`, `_volume_matrix`, `comprehensive_peak_tracking.csv`, `series_analysis_tidy.csv`, per-spectrum CSVs |
| `dynamixs t1t2` | Mono-exp T1/T2 fit | `--input <matrix.csv>` `--out` `--exp {T1,T2}` | `<prefix>_fit_results.txt` + JSON |
| `dynamixs methyl-t2` | Bi-exp Tugarinov-Kay methyl T2 | `--input` `--out` | `<prefix>_fit_results.txt` + JSON |
| `dynamixs hetnoe` | I_sat/I_unsat per residue | `--sat` `--unsat` `--out` | `<prefix>_hetnoe.csv`, `<prefix>_hetnoe.pdf` |
| `dynamixs t1rho` | Fit a T1rho series and convert it to R2 (tilt-angle corrected, per residue) | `--input` `--t1` `--peaks` `--omega1` `--carrier` `--theta` `--out` | `<prefix>_r2_from_t1rho.csv` |
| `dynamixs density` | Reduced spectral density mapping | `--input <table.csv>` `--out` `--field1-freq` | `<prefix>_results.csv` (or `_basic`/`_detailed` for `--dual`) + `_plots.pdf` |
| `dynamixs modelfree` | T1/T2→hetNOE→density→Lipari-Szabo | `--f1-t1 --f1-t2 --f1-noe-sat --f1-noe-unsat --field1-freq --out` | model-free CSVs, per-field hetNOE.pdf, density plots |
| `kd` | Kd titration (CSP quadratic / intensity decay). `--survey` first, then `--residues` to fit a human-chosen subset | `--input <titration.csv>` `--out` `--p0 <conc>` `[--conc-units {absolute,equivalents}]` `[--survey]` `[--residues FILE|list]` `[--selection PATH]` `[--csp-sigma-multiple 1.0]` `[--kd-outlier-z 3.0]` | `data/<prefix>_kd_fit_data.json`, `<prefix>_kd_results.txt`; with `--survey`: `<prefix>_{csp,intensity}_vs_sequence.pdf`, `data/<prefix>_survey.csv`, and `<prefix>_residues.txt` **beside the input** |
| `export kd` | Figures + summary from a kd fit JSON | `--json` `--out` `[--kind curves,ref-bars,kd-bars,global-fit]` `[--prefix]` | `[prefix_]summary.csv`, `[prefix_]<obs>_fits.pdf` / per-residue PNGs (curves, per fitted obs), `[prefix_]<obs>_ref_vs_point.pdf` (ref→point bars, both obs by default, PDF only), `[prefix_]<obs>_kd_vs_residue.pdf` (per-residue Kd + global line, PDF only), `[prefix_]<obs>_global_fit.pdf` (per-residue data + shared-Kd curve, R²(global) per panel, PDF only) |
| `project inventory` | List a bundle's contents | `<bundle.lunaNMR>` | stdout listing |
| `project remove` | Delete bundle-relative paths | `<bundle> <rel/path>…` | (mutates bundle) |
| `batch` | Folder-wide detect + Voigt/PS2D fit | `<folder>` (flags pass through to `batch_processing`) | per-spectrum fit outputs |

**Doing T1/T2/hetNOE?** Follow the Recipe below in order — it is the canonical workflow and
states every requirement that is not visible from `--help`. The Copy-paste blocks after it are
just the flag surface.

## Recipe: ¹⁵N relaxation → model-free

The dominant workflow. Follow it in order — do not reinvent it. Each step states its own
non-obvious requirement.

**0. Diagnose first.** Read-only, seconds, changes nothing. It answers the questions the
later steps depend on — which hetNOE plane is saturated, whether any repeat acquisitions
disagree, whether every peak list matches its own spectra, and whether the residue sets
agree across fields.
```bash
python -m lunaNMR diagnose <dataset_root> [--quick] [--format json]
```
A wrong peak list does not raise — it produces zero heights and R²≈0.2 shoulder fits that
read as noisy data, after minutes of integration. This is what catches it first. For a single
experiment the same checks are available as `series --dry-run --deep`.

**1. Every experiment gets its OWN peak list.** hetNOE lists are typically offset ~0.07 ppm ¹H / ~0.7 ppm ¹⁵N from the T1/T2 frame — 1× and 10× the 0.070 ppm search window. Using the wrong one gives Height 0 or R²≈0.2 shoulder fits that read as "noisy data". Never widen the window to compensate.

**2. Name the hetNOE planes before integrating.** Two planes are not a delay series, so the matrix columns become the filename stems. Symlink to `*_saturated.ft` / `*_unsaturated.ft`, deciding which is which by the **measured** intensity ratio over high-S/N peaks (saturated is the lower one; 0.75–0.85 typical, never >1). `001`/`002` filenames say nothing.

**3. Integrate — one `series` call per experiment** (6 calls for two fields):
```bash
python -m lunaNMR series --spectra <exp_dir> --peaks <its_own_list.txt> \
  --out RUN/series/<name> --mode time --peak-source reference --parallel --format json
```
Matrix columns are series values (ms), not filenames. **`series_metadata.json`** beside the matrices maps every column → spectrum → value; read it rather than inferring, and check its `n_value_unparsed` is 0.

**4. Build the hetNOE inputs.** Split the hetNOE matrix's two columns into headerless `residue,intensity,error` files. **Write the third column** — without it the error defaults to a flat `noe*0.02` (≈0.016), ~3× too tight, which over-weights hetNOE and fakes tight S² errors; a realistic floor is ~0.044 relative per plane → hetNOE_err ≈ 0.05. Gate residues using `series_analysis_tidy.csv`'s `quality` (both planes Excellent/Good, both intensities > 0, ratio in −0.2…1.05).

**5. Fit T1/T2** — usually unnecessary. `modelfree` refits the matrices internally and already writes `<prefix>_field{1,2}_T{1,2}_fit_results.txt` (per-residue T and T_err), the `_fig{1,2,3}.pdf` decay plots, and R1/R2 in `<prefix>_spectral_density_basic.csv`. Run `t1t2` separately only for a different `--error-method`, or to apply your own fit filter before deciding which residues to keep:
```bash
python -m lunaNMR dynamixs t1t2 --input RUN/series/T1_600/peak_intensity_matrix.csv \
  --out RUN/fits/T1_600 --exp T1 --field-freq 600 --time-units ms --format json
```
Keep fits with finite `T`, `T > 0`, `T_err < T`. Its residue count will not match `modelfree`'s.

**6. Model-free — run `--dual` AND each field alone.** The single-field runs are the control that tells you whether to believe the dual one (step 7).
```bash
python -m lunaNMR dynamixs modelfree --dual \
  --f1-t1 … --f1-t2 … --f1-noe-sat … --f1-noe-unsat … --field1-freq 600 \
  --f2-t1 … --f2-t2 … --f2-noe-sat … --f2-noe-unsat … --field2-freq 700 \
  --out RUN/mf --prefix dual --format json
```
Per-residue results land in `<prefix>_spectral_density_basic.csv`.

**6b. T1ρ instead of T2.** A T1ρ series measures R1ρ, not R2: `R1ρ = R1·cos²θ + R2·sin²θ`.
Recovering R2 needs R1 *and* the spin-lock geometry, and the tilt angle differs per residue
because residues sit at different offsets from the spin-lock carrier — one nominal angle is
not enough (a residue 22 ppm off carrier at ω₁ = 2 kHz tilts to ~56°, not 90°, and its R2 is
~40% higher than the uncorrected value).
```bash
python -m lunaNMR dynamixs t1rho --input <t1rho_matrix.csv> --t1 <t1_matrix.csv> \
  --peaks <list.txt> --omega1 <Hz> --carrier <15N ppm> --theta <deg> \
  --field-freq 600 --out RUN/t1rho --prefix field1 --format json
```
`--omega1` (cnst27), `--carrier` and `--theta` (cnst28) are **acquisition parameters that are
not in the files** — get them from whoever ran the experiment; there is nothing to measure
them from. The output feeds model-free in place of a T2 series:
```bash
python -m lunaNMR dynamixs modelfree --f1-t1 … --f1-r2-table RUN/t1rho/field1_r2_from_t1rho.csv \
  --f1-noe-sat … --f1-noe-unsat … --field1-freq 600 --out RUN/mf
```
`--f{1,2}-t2` and `--f{1,2}-r2-table` are alternatives; giving both is an error. Everything
downstream is identical, so the QC table below still applies — and comparing R2-from-T1ρ
against a directly measured R2 is the sharpest validation available: agreement confirms both,
disagreement localises exchange (R1ρ and R2 weight it differently).

**7. Physical QC — check these before believing anything.** R1 = 1000/T1, R2 = 1000/T2 (ms), per residue, common set.

| check | expect | if it fails |
|---|---|---|
| R1(high)/R1(low) | **field-pair dependent** — table in [`CLI_AGENTS_DEEP/RELAXATION_PLAYBOOK.md`](CLI_AGENTS_DEEP/RELAXATION_PLAYBOOK.md) | peak identity / referencing / scaling broken |
| R2(high)/R2(low) | **field-pair dependent** — same table | a T2 problem, or the two datasets differ in temperature |
| τc = (1/4πν_N)·√(6·R2/R1 − 7), per field | agree within ~5% | there is no single global τc — the dual fit's premise is false |
| hetNOE median | 0.75–0.85 for a rigid fold | plane identity or gating is wrong |
| Rex>1 fraction, dual vs single-field | comparable | dual ≫ single (and χ² ~25× higher) = the dual fit is absorbing a τc mismatch, not reporting chemistry. Take per-residue Rex from the single-field runs. |
| S² spread | `real = √(observed_sd² − mean_err²)` | if noise is ~half the variance (cross-field S² r ≈ 0.2 with matching medians), report a uniform S² with outliers, not a per-residue profile |

**R1 agreeing while R2 disagrees isolates the fault to T2** — the most useful single diagnostic here. Report such a failure; do not tune it away.

## Copy-paste blocks
The flag surface. For T1/T2/hetNOE the Recipe above is the order to run them in.
```bash
# Pre-flight: whole dataset, or one experiment folded into its dry-run
python -m lunaNMR diagnose <root> [--quick] [--sample] [--mode {time,titration}]
python -m lunaNMR series --spectra 600/T1 --peaks 600/T1/peaks.txt --out RUN/T1 \
  --dry-run --deep [--quick] [--sample]

# Series (titration adds --mode titration; --mode time is default)
python -m lunaNMR series --spectra 600/T1 --peaks 600/T1/peaks.txt --out RUN/T1 \
  [--mode {time,titration}] [--peak-source {reference,cascade,detected,independent}] [--parallel]

# T1/T2 fit  (exit 1 if nothing fits)
python -m lunaNMR dynamixs t1t2 --input RUN/T1/peak_intensity_matrix.csv --out OUT --exp T1 \
  [--time-units {s,ms,us}] [--field-freq 600] [--error-method {analytical,bootstrap}] [--bootstrap N] [--no-json]

# hetNOE  (--sat/--unsat headerless: residue,intensity[,error])
python -m lunaNMR dynamixs hetnoe --sat noe_sat.csv --unsat noe_unsat.csv --out OUT [--prefix field1]

# Dual-field model-free (per-series delay units matter — see gotchas)
python -m lunaNMR dynamixs modelfree --dual \
  --f1-t1 600/T1/...matrix.csv --f1-t2 600/T2/...matrix.csv \
  --f1-noe-sat 600/noe_sat.csv --f1-noe-unsat 600/noe_unsat.csv --field1-freq 600 \
  --f2-t1 800/T1/... --f2-t2 800/T2/... --f2-noe-sat ... --f2-noe-unsat ... --field2-freq 800 \
  --out RUN/mf [--f1-t1-units s] [--f1-t2-units ms] [--method {single_087,single_jwh,dual_087,dual_jwh}]

# Kd  (fits csp AND intensity by default)
python -m lunaNMR kd --input series_analysis_tidy.csv --out OUT --p0 200 \
  [--observable csp,intensity] [--alpha 0.14] [--intensity-from {height,volume}] \
  [--conc 0,10,25] [--intensity-scale s0,s1,...] [--bootstrap N]

# Validate before running: prints the plan, runs nothing, exits 1 if inputs missing
python -m lunaNMR <cmd> ... --format json --dry-run
```

## Machine conventions
- **`--format json`** → stdout is one JSON summary object, nothing else; all engine chatter goes to stderr (fd-level, so spawn workers stay clean). Parse stdout for `n_fitted`/`n_successful`, output paths, etc.
- **`--dry-run`** → validates required inputs exist, prints the plan, runs nothing. Exit 1 if any input missing, else 0.
- **`--parallel`** → `series` only (two-pass processor, ~2.7×). **`--no-parallel` → `density` only.** `modelfree` also parallelizes internally but hard-codes `use_multiprocessing=True`, so its CPU use cannot be bounded — `modelfree --no-parallel` exits 2 with `unrecognized arguments`.
- **Exit codes:** `2` = no/invalid subcommand (help to stderr). `1` = bad input (missing file, malformed CSV/JSON, unparseable delay label — reported as `error: …` on stderr, no traceback) OR dry-run with a missing input. `series`/`t1t2`/`methyl-t2` also return `1` when nothing fit. `hetnoe`/`density`/`modelfree`/`kd`/`export` return `0` on any successful run even if some residues failed — check the JSON `n_successful`/`n_fitted`, don't rely on exit code for partial failure.

## Input file shapes
- **Peak list** (`series --peaks`): header `Assignment, Position_X, Position_Y[, Height]` (comma or tab; column aliases standardized).
- **Series spectra** (`series --spectra`): folder or glob; files auto-discovered by extension (`ft ser ft2 ft3 pipe ucsf`), naturally sorted.
- **Relaxation matrix** (`t1t2/methyl-t2 --input`): a series `peak_intensity_matrix.csv` / `comprehensive_peak_tracking.csv`; delay columns parsed by `parse_delay_column` (explicit `ms`/`s`/`us` + optional repeat marker `b` or `_2`, e.g. `..._300msb`). A matrix **produced by `series`** has already had the unit stripped and normalised to ms — its columns are bare numbers (`8`, `2400`, `102_2`); see the units gotcha.
- **hetNOE sat/unsat** (`--sat`/`--unsat`, modelfree `--f{1,2}-noe-*`): headerless `residue,intensity[,error]`. Split from the hetNOE series matrix's two columns — which are the plane **filename stems**, not `…_sat`/`…_unsat`, unless you named the files that way. **Supply the error column:** with only two columns the hetNOE error falls back to a flat `noe * 0.02` (≈0.016 at hetNOE 0.78), ~3× tighter than a realistic ~0.044-per-plane floor gives (≈0.05), which over-weights hetNOE in the model-free fit.
- **Density table** (`density --input`): `Residue,R1,R1err,R2,R2err,hetNOE,hetNOEerr`.
- **Kd input** (`kd --input`): `series_analysis_tidy.csv`, `comprehensive_peak_tracking.csv`, or an intensity matrix (intensity-only, no CSP). Ligand conc from CSV point labels unless `--conc` overrides.
- **Titration spectra** (`series --mode titration`): the point value is a trailing `_<value>` where the value carries a **decimal separator** — `o` or `.` (`_0o0`, `_0o5`, `_2.4`), optionally followed by an extension. **Digits alone do not parse**: `sample_2.ft`, `ref_100.ft` and `experiment_002.ft` all yield no value and are silently dropped from the series. This is stricter than the time-mode delay parser, which keys off the unit (`ms`/`s`/`us`) instead.
- **Kd fit JSON** (`export kd --json`): the self-contained `…_kd_fit_data.json` from `kd` (embeds per-point series + `metadata.protein_conc`; `metadata.name` records the `kd --prefix` used to create it, but `export kd` does **not** read it — pass `export kd --prefix` explicitly to prefix figures/summary.csv, independent of the JSON's own filename).

## Output file shapes (what you parse back)

- **`series` matrices**: `Peak_Number,Assignment,Reference_X,Reference_Y` + one column per spectrum, labelled with the parsed delay (`8`, `2400`, `102_2`) or, when none parses, the filename stem (`hetnoe_600_saturated`). Per spectrum, so one folder can carry both kinds.
- **`series_analysis_tidy.csv`**: `spectrum_name, assignment, peak_number, ppm_x, ppm_y, height, volume, snr, quality, r_squared`. **The only place per-peak fit quality reaches the series output** — gate residues here, the matrices are intensities only.
- **`per_spectrum_results/<original_filename>.csv`**: `…,Detected_Intensity,Height,Volume,LW_X,LW_Y,R_Squared,Quality`. The only output where the filename→data link survives.
- **`series_metadata.json`**: **the column → spectrum → delay map**; nothing else carries it. `value` is `null` when the filename had no parseable delay — that column is named after the file stem instead, and `n_value_unparsed` counts them. `repeat_scale` is measured on the fitted intensities and **reported, never applied**. **It is `null` when the series has no repeat acquisitions** — the common case — and absent entirely when no intensity matrix was produced, so `if 'repeat_scale' in meta:` is True with a `null` payload and crashes on `['ratio']`. Use `meta.get('repeat_scale') or {}`.
  ```json
  {"series_mode": "time", "value_units": "ms", "n_spectra": 14, "n_value_unparsed": 0,
   "columns": [{"column": "8",     "spectrum": "T2_sample_8ms.ft",  "value": 8.0},
               {"column": "8_2",   "spectrum": "bT2_sample_8ms.ft", "value": 8.0},
               {"column": "T2_odd","spectrum": "T2_odd.ft",         "value": null}],
   "repeat_scale": {"ratio": 0.8731, "scale": 1.1454, "deviation_percent": 12.7,
                    "needs_normalisation": true,
                    "per_value": {"8.0": {"first": "8", "second": "8_2",
                                          "ratio": 0.8592, "n_peaks": 29}}}}
  ```
- **`processing_summary.csv`**: processing order — **not** the matrix column order.
- **model-free `<prefix>_spectral_density_basic.csv`**: one row per residue; column names differ single vs dual (see gotchas).
- **`diagnose --format json`** (and `series --dry-run --deep`'s `checks` key): `experiments[]` each with `spectra[]` (`dx`, `dy`, `capture`, `capture_rate`, `median_snr`, `weak`, `decayed`, `is_reference`), `delays` (`parsed`, `unparsed`, `repeats`), `peak_list`, `hetnoe`, `subseries`; plus a flat `findings[]` of `{severity, check, message}` with severity `FAIL`/`WARN`.

## Gotchas (silent-corruption risks)
- **modelfree units rescale, t1t2 units don't.** `modelfree --f{1,2}-t{1,2}-units {ms,s,us}` (default `ms`) DRIVE the T→R conversion — a T1 in seconds with a T2 in ms needs `--f1-t1-units s`, else R1 is off by 1000×. Contrast `t1t2 --time-units` (default `s`) which only *labels* output.
- **CSP needs a fresh series run.** Older series output rounded positions to 0.1 ppm; CSP needs 4-decimal precision. Re-run `series` before `kd` CSP. Intensity fits work on existing data.
- **Kd is always bounded ≥ 0** internally (unbounded lands on a degenerate negative-Kd minimum). Don't post-process assuming otherwise.
- **Non-finite fit values serialize to JSON `null`.** A degenerate Kd fit (e.g. only 3 points) is `success=true` with finite `Kd` but `Kd_err`/`r_squared` = `null`. Treat `Kd`/`Kd_err`/`r_squared` as possibly `None` — legitimate, not a failure.
- **Search-window default is 0.070 ppm** (1H and 15N). Don't widen 15N for low-S/N spectra (e.g. saturated hetNOE) — wide windows fit peaks to noise → Height 0 → NOE≈0 → residues silently dropped.
- **`--method` field prefix follows `--dual`.** In `modelfree` field count derives from `--dual`; you only pick the `087` vs `jwh` variant (default `087`). A mismatched `single_*`/`dual_*` prefix is ignored.
- **`--intensity-scale` scales height/volume only** (per-point scan-count correction), never positions/CSP; a uniform scale cancels in the I/I₀ ratio.
- **Which duplicate delay gets `_2` is a case-sensitive sort, not semantics.** `sort_files_with_sequence` orders by `(delay, basename)`; the first gets the bare label. `..._600ms` beats `..._b_600ms`, and `T2_8ms` beats `bT2_8ms` only because `T` < `b` in ASCII — a `BT2_` sub-series would steal the bare label. It is not the sort `--spectra` uses to discover files, so `processing_summary.csv` order is not column order either. **Read `series_metadata.json` instead of reasoning about this.**
- **The series delay parser needs an explicit unit at the end of the stem.** `ms`, `s` and `us` all parse and all normalise to ms (`_2.4s` → `2400.0`, `_2400us` → `2.4`), optionally followed by a repeat marker — a single letter (`_300msb`) or a numeric suffix (`_8ms_2`) — which names the same delay as the original, so the two collide into `300`/`300_2`. A **bare number** (`_2400`) does **not** parse: ms or s is unknowable and guessing puts R1 out by 1000×. Unparsed spectra fall back to a stem-named column, mixing label kinds in one matrix; `series_metadata.json`'s `n_value_unparsed` counts them after the fact and `diagnose` reports them before you run.
- **`modelfree --f{1,2}-t{1,2}-units` applies to matrices `series` did not write.** A `series`-produced matrix is already normalised to ms (`value_units` in the sidecar), so the default `ms` is right. Hand-built or DynamiXs-format tables still need the flag — a T1 in s with a T2 in ms puts R1 out by 1000×.
- **model-free CSV columns differ single vs dual.** Dual: `R1_f1`/`R1_f2`, `Rex_f1_err`, `Rex_f1_95CI`, … Single: the bare names `R1`, `Rex_err`, `Rex_95CI` **plus `_f1` aliases of each**, so a reader written against the dual spelling works on both. Shared, never suffixed: `Residue, S2, tc, te, S2_err, te_err, chi2, fit_success`. CSVs written before the aliases existed carry only the bare names, and a `_f1` selection on those returns an **empty selection** rather than raising — which reads as "no exchange anywhere".
- **`export kd` curves vs ref-bars differ on observable.** `curves` (`<obs>_fits.pdf`) render only fitted observables (follow the JSON / `--observable`); `ref-bars` (`<obs>_ref_vs_point.pdf`) are model-free (raw series) and default to BOTH csp+intensity even for a single-observable fit — `--observable` still restricts. Ref-bars need embedded per-point `series` + ≥2 points.

## Kd workflow: survey → choose → fit

Fitting every residue to a Kd is meaningless — most residues do not bind. The agent
proposes, a human disposes, and only then is a Kd reported.

```
1. python -m lunaNMR kd --survey --input tidy.csv --out RESULTS --p0 50 --conc-units equivalents
2. read RESULTS/<prefix>_csp_vs_sequence.pdf; edit <input_dir>/<prefix>_residues.txt
3. python -m lunaNMR kd --input tidy.csv --out RESULTS --p0 50 --conc-units equivalents \
       --residues <input_dir>/<prefix>_residues.txt
4. python -m lunaNMR export kd --json RESULTS/data/<prefix>_kd_fit_data.json --out RESULTS
```

**Selection file.** One residue per line; `#` comments a line out, deleting the `#` puts it
back. Trailing comments are evidence, never exclusions. It lives **beside the input CSV**,
not in `--out`, so a different `--out` cannot lose it; `--selection` overrides. Re-running
`--survey` MERGES: human decisions survive, evidence refreshes, disagreements are marked
`[kept]`/`[dropped]` and never re-suggested.

**Output layout.** `--out` holds `*.pdf` and `*_kd_results.txt`; every CSV and JSON is in
`<out>/data/`. `export kd --json` therefore points into `data/`.

**Why a residue is missing from a shared CSP Kd.** Five gates, each naming its exclusions
in `metadata.csp_pool_excluded` — `dummy_*`, `quality='Failed'`, CSP below the significance
threshold, R² < 0.8, Kd outside the resolvable window or a robust-z outlier.

**Read `metadata.quality_warnings` before quoting any Kd.** When the *typical* residue fits
a Kd outside what the titration can resolve, the tool says so in plain language and refuses
to exclude outliers — there is no trustworthy centre to measure deviation from, and the
whole dataset is unusable rather than a few residues. `--kd-outlier-z 0` disables the
statistical gate but never this verdict.

**Survey thresholds — flags, with the measurement behind each default.** The reason a
default is what it is, is what you need before changing it.

| flag | default | measured basis |
|---|---|---|
| `--csp-sigma-multiple` | `1.0` | multiples of the trimmed CSP spread at the last titration point a residue must exceed to enter the shared fit. `2.0` is a common stricter variant |
| `--kd-outlier-z` | `3.0` | robust z (median/MAD on log10 Kd). `0` disables the statistical test but never the credibility verdict or the resolvable-window gate |
| `--ref-max-ratio` | `10.0` | a reference intensity this far below the residue's own series max is a broken denominator. Across 129 residues legitimate ratios topped out at 1.30 and the next value was 49.9 — 10 sits in the empty band |
| `--dd-runaway-ratio` | `10.0` | fitted plateau over the residue's own largest CSP. Cited as evidence, never used to exclude: the flagged fraction ran 17% on one dataset and 43% on another, so it describes the experiment, not the residue |
| `--noise-quantile` | `0.25` | quantile of the max-CSP distribution below which a residue counts as a non-mover. Derived per dataset, because the noise floor is a property of the spectra |

All five persist per dataset and follow the same precedence as everything else.

**Settings persist per dataset.** `conc_units`, `csp_sigma_multiple`, `kd_outlier_z`,
`noise_quantile`, `dd_runaway_ratio` and `ref_max_ratio` are written to
`<input_dir>/<prefix>_kd_params.json` and read back automatically. Precedence: an explicit
flag beats the persisted value beats the schema default.

**⚠ Equivalents are not concentrations.** Titration filenames commonly encode equivalents
(`..._3o0.ft` = 3.0 eq). Read as absolute they give a confident Kd wrong by the factor
`--p0`, with `success=True` and no error. Pass `--conc-units equivalents`, or `--conc` with
explicit values. The default is `absolute`, so this trap is avoidable, **not closed**.

## Deep guides

`docs/CLI_AGENTS_DEEP/` holds the long-form runbooks — the phase structure, the physical QC
thresholds, and the traps that are invisible from the CLI surface. Read the playbook for
your analysis type before starting; it is the difference between running the commands and
knowing whether the answer is real.

| file | for |
|---|---|
| `AFFINITY_PLAYBOOK.md` | binding titrations → Kd: diagnose → survey → fit → QC, and the traps (equivalents read as concentrations, broken reference points, spread-based gates that widen to admit the garbage inflating them) |
| `AFFINITY_PROMPT.md` | a ready task prompt for an end-to-end titration run, with timing and the two mandatory stop points |
| `RELAXATION_PLAYBOOK.md` | ¹⁵N T1/T2/hetNOE → model-free: the diagnose/run split, field-pair ratio bands, τc cross-checks |
| `RELAXATION_PROMPT.md` | a ready task prompt for an end-to-end relaxation run |

Both playbooks share one rule worth stating here: **diagnose, report, stop.** Corrections
change the science, so raise them with the measured number and let the user decide.

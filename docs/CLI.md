<!-- ABOUTME: Command reference for the headless lunaNMR CLI (python -m lunaNMR). -->
<!-- ABOUTME: Subcommand flags plus the end-to-end dual-field model-free workflow. -->

# lunaNMR CLI Reference

Headless access to the analysis pipeline. Every reasonable GUI analysis is scriptable — no
display required.

```bash
python -m lunaNMR <subcommand> [flags]      # console script `lunanmr` needs `pip install -e .`
python -m lunaNMR --version
python -m lunaNMR --help
```

The CLI imports no Qt: `lunaNMR/__init__.py` is lazy, and matplotlib is forced to `Agg`.

## Global conventions

- `--format {text,json}` — machine-readable run summary on stdout (`json` stays clean even with
  multiprocessing-spawn workers).
- `--dry-run` — validate inputs and print the plan without running; exits non-zero if an input
  is missing.

## Subcommands

| Command | Purpose |
|---|---|
| `diagnose` | Read-only pre-flight over a dataset before running anything |
| `series` | Process a multi-spectrum series/titration into intensity/position matrices |
| `dynamixs t1t2` | Mono-exponential T1 or T2 relaxation fit |
| `dynamixs methyl-t2` | Bi-exponential (Tugarinov-Kay) methyl T2 fit |
| `dynamixs hetnoe` | Heteronuclear NOE = I_sat / I_unsat per residue (+ QC plot) |
| `dynamixs t1rho` | Fit a T1ρ series and convert it to R2 |
| `dynamixs density` | Reduced spectral density mapping (single or dual field) |
| `dynamixs modelfree` | Integrated pipeline: T1/T2 fit → hetNOE → density → Lipari-Szabo |
| `kd` | Kd titration (CSP quadratic / intensity decay) |
| `export kd` | CSP/intensity figures + summary from a saved Kd fit JSON |
| `project inventory` / `remove` | Inspect / prune a `.lunaNMR` bundle |
| `batch` | Folder-wide peak detect + Voigt/PS2D fit |

### diagnose
```bash
python -m lunaNMR diagnose <dataset_root> [--quick] [--sample] [--mode {time,titration}]
```
Read-only: reads spectra, fits nothing, writes nothing. Per experiment folder it reports the
parsed delays and repeats, peak-list hygiene, the registration offset and capture rate of
every spectrum, which hetNOE plane is saturated, and the scale of any repeat acquisitions —
then compares the assignment sets across experiments and prints `FAIL`/`WARN` lines.

Run it before the pipeline. A peak list belonging to a different experiment does not raise;
it produces zero heights and R²≈0.2 shoulder fits that read as noisy data, and you only find
out after several minutes of integration.

Two things it does that nothing else can. **Registration** slides the whole peak list over
the spectrum and reports the offset that best aligns it — the pipeline assumes the list is
right. **Cross-experiment residue sets** are compared, because merging is an intersection, so
the smallest set is the ceiling on every downstream result.

Everything about labels and identity — delay values, peak-list contents — is parsed with the
same code `series` uses, so the report cannot certify something the run will not reproduce.
Intensities are deliberate approximations (the maximum in a small box, not a fit); a
sub-series ratio from here is an estimate to be re-measured on fitted intensities before it
is applied.

`--quick` uses a coarser registration grid. `--sample` assesses only the first and last
spectrum of each experiment instead of all of them.

**Capture rate is only judged at the reference point** — the shortest delay, where nothing
has decayed yet. A shortfall there means the list or the data is wrong. Further down a
relaxation series the same number just means the experiment worked, so it is reported without
a warning.

For a single experiment, the same checks fold into its dry-run:
```bash
python -m lunaNMR series --spectra <dir> --peaks <list> --out <out> --dry-run --deep
```
`--deep` adds a `checks` object to the dry-run summary and leaves every existing key
untouched. It reads the spectra, so it takes seconds rather than being instant — plain
`--dry-run` is unchanged. The cross-experiment comparison is not available this way; that
needs `diagnose`.

### series
```bash
python -m lunaNMR series --spectra <folder> --peaks <list.txt> --out <dir> \
    [--mode {time,titration}] [--peak-source {reference,cascade,detected,independent}] [--parallel]
```
- Peak list: `Assignment, Position_X, Position_Y[, Height]` header.
- `--mode time` (relaxation, default) or `titration` (forces cascade tracking + wider margins).
- Filename values: `time` mode reads a trailing delay with an explicit unit (`_300ms`, `_2.4s`,
  `_2400us`). `titration` mode reads a trailing `_<value>` whose value carries a **decimal
  separator** — `o` or `.` (`_0o0`, `_0o5`, `_2.4`). Digits alone (`sample_2.ft`, `ref_100.ft`)
  parse in neither mode; those spectra get a stem-named column instead of a value.
- `--peak-source reference` (default) holds list positions fixed across the series — correct for
  T1/T2/hetNOE.
- Writes `peak_intensity_matrix.csv`, `_volume_matrix`, `comprehensive_peak_tracking.csv`,
  `series_analysis_tidy.csv`, `series_metadata.json`, per-spectrum CSVs.

**Repeat acquisitions are measured for you.** When two spectra share a series value,
`series_metadata.json`'s `repeat_scale` gives the median ratio of the second to
the first over the strong peaks, per shared value and overall, plus the `scale` that would
cross-normalise them and a `needs_normalisation` flag (>5% off). It is measured on the
intensities the run just fitted, which is the number worth acting on — a pre-flight check can
only estimate it from box maxima. It is **reported and never applied**: rescaling one
sub-series changes the science and stays your decision. The key is always written but is
`null` when there are no repeats, so read it as `meta.get('repeat_scale') or {}`.

**Column labels, and how to get back to the filename.** The CLI always extracts delays, so
each matrix column is labelled with the delay parsed from the filename (`8`, `102`, `2400`)
when one parses, and with the filename stem (`hetnoe_600_saturated`) when none does. Repeated
delays get `_2`, `_3`. **`series_metadata.json`**, written beside the matrices, holds the map: per column, the source
spectrum, the parsed value, and whether it parsed at all. Use it — `processing_summary.csv`
lists spectra in processing order, which is **not** the column order (under `--parallel` the two
differ). This matters as soon as you need to treat one sub-series differently from another
(cross-normalising a re-acquired tail, say).

Values are normalised to **ms**: `_2.4s` becomes `2400.0` and `_2400us` becomes `2.4`. The
unit must be explicit (`ms`, `s` or `us`) at the end of the stem, optionally followed by a
repeat marker — a single letter (`_300msb`) or a numeric suffix (`_8ms_2`), both of which
name the same delay as the original and so collide into `300`/`300_2`. A **bare number**
(`_2400`) does not parse: ms or s is unknowable, and guessing would put R1 out by 1000×.
Those spectra get a column named after the file stem instead; the sidecar's
`n_value_unparsed` counts them and should be 0.

Within one delay group the file that gets the **bare** label is chosen by a plain
case-sensitive sort of the basename (`sort_files_with_sequence`), so `..._600ms` precedes
`..._b_600ms`, and `T2_8ms` precedes `bT2_8ms` only because `T` < `b` in ASCII. A sub-series
prefixed `BT2_` would sort first and take the bare label instead. Don't assume the "original"
owns it.

**Per-peak fit quality is in `series_analysis_tidy.csv`**, not in the matrices — columns
`quality` (`Excellent`/`Good`/`Fair`/`Poor`/`Failed`) and `r_squared`, one row per
(spectrum × peak). Gate residues from there.

**hetNOE folders:** two planes are not a delay series, so the columns come out as the
filename stems. Name (or symlink) the planes explicitly — `*_saturated.ft` /
`*_unsaturated.ft` — before integrating, and decide which is which by measuring the intensity
ratio over high-S/N peaks (saturated is lower, typically 0.75–0.85), never from a `001`/`002`
filename.

### dynamixs t1t2 / methyl-t2
```bash
python -m lunaNMR dynamixs t1t2 --input <matrix.csv> --out <dir> --exp {T1,T2} \
    [--prefix field1] [--field-name field1] [--field-freq 600] [--time-units {s,ms,us}] \
    [--error-method {analytical,bootstrap}] [--bootstrap 1000] [--no-json]
python -m lunaNMR dynamixs methyl-t2 --input <matrix.csv> --out <dir> \
    [--prefix field1] [--field-name field1] [--field-freq 600] [--time-units {s,ms,us}] \
    [--error-method {analytical,bootstrap}] [--bootstrap 1000]
```
Input = a series `peak_intensity_matrix.csv` / `comprehensive_peak_tracking.csv`. Delay column
headers are parsed by `parse_delay_column` (explicit `ms`/`s`/`us` unit + optional repeat
marker `b` or `_2`). Exit code 1 if nothing fits.

### dynamixs hetnoe
```bash
python -m lunaNMR dynamixs hetnoe --sat <sat.csv> --unsat <unsat.csv> --out <dir> [--prefix field1]
```
`--sat`/`--unsat` are headerless `residue,intensity[,error]` CSVs. Writes `<prefix>_hetnoe.csv`
and a QC plot `<prefix>_hetnoe.pdf` (I_sat/I_unsat vs residue number).

### dynamixs t1rho
```bash
python -m lunaNMR dynamixs t1rho --input <t1rho_matrix.csv> --t1 <t1_matrix.csv> \
    --peaks <list.txt> --omega1 <Hz> --carrier <ppm> --theta <deg> \
    [--field-freq 600] [--time-units ms] [--error-method {analytical,bootstrap}] \
    --out <dir> [--prefix field1]
```
A T1ρ experiment measures R1ρ, which mixes longitudinal and transverse relaxation:
`R1ρ = R1·cos²θ + R2·sin²θ`. Recovering R2 therefore needs the T1 series as well as the
spin-lock geometry, which is why `--t1` is required.

**The tilt angle is computed per residue.** Residues sit at different offsets from the
spin-lock carrier, so they are tilted differently and a single nominal θ is wrong for most of
them — hence `--peaks`, which supplies each residue's ¹⁵N shift. At ω₁ = 2 kHz a residue
22 ppm off carrier tilts to ~56° rather than 90°, and its R2 comes out ~40% higher than the
uncorrected value.

`--omega1` (cnst27), `--carrier` and `--theta` (cnst28) are acquisition parameters. They are
**not recoverable from the spectra or the filenames** — they have to come from whoever ran the
experiment.

Writes `<prefix>_T1_fit_results.txt`, `<prefix>_T1rho_fit_results.txt` and
`<prefix>_r2_from_t1rho.csv` (`residue, R1, R1_err, R1rho, R1rho_err, theta, R2, R2_err,
T2, T2_err`). Feed that table to `modelfree` with `--f{1,2}-r2-table` in place of
`--f{1,2}-t2`; the two are alternatives and giving both is an error.

When a dataset has both a T2 and a T1ρ series, run both and compare R2. Agreement validates
each independently; disagreement localises exchange, since R1ρ and R2 weight it differently.

### dynamixs density
```bash
python -m lunaNMR dynamixs density --input <table.csv> --out <dir> [--field1-freq 600] \
    [--dual --input2 <table2.csv> --field2-freq 700] [--no-087] [--monte-carlo --n-samples 1000] \
    [--rnh 1.015] [--csa -172] [--no-parallel] [--no-plot] [--prefix spectral_density]
```
Input table columns: `Residue, R1, R1err, R2, R2err, hetNOE, hetNOEerr`. Outputs
J(0)/J(ωN)/J(ωH), S², τc, τe, Rex.

### dynamixs modelfree
Runs the whole chain per field from raw series matrices.
```bash
python -m lunaNMR dynamixs modelfree [--dual] \
    --f1-t1 <T1_matrix.csv> --f1-t2 <T2_matrix.csv> \
    --f1-noe-sat <sat.csv> --f1-noe-unsat <unsat.csv> --field1-freq 600 \
    [--f2-t1 ... --f2-t2 ... --f2-noe-sat ... --f2-noe-unsat ... --field2-freq 700] \
    [--f1-t1-units {ms,s,us}] [--f1-t2-units ...] [--f2-t1-units ...] [--f2-t2-units ...] \
    [--method {single_087,single_jwh,dual_087,dual_jwh}] [--prefix modelfree] --out <dir>
```
- **Field count follows `--dual`** — the `--method` prefix is derived; you only pick the
  `087`/`jwh` variant. Omit `--method` for the `087` default.
- **Per-series time units:** the T→R conversion needs each series' delay units. A T1 series in
  seconds with a T2 in ms → `--f1-t1-units s`. Cleanly `ms`-labelled headers (`..._300msb`)
  need no flag, and neither does a matrix produced by `series` — those are already normalised
  to ms (`series_metadata.json` records `value_units`). The flag is for tables `series` did not
  write: a hand-built or DynamiXs-format T1 in seconds beside a T2 in ms puts R1 out by 1000×.
- **Outputs:** `<prefix>_field{1,2}_hetnoe.pdf`, `<prefix>_spectral_density_basic.csv` (one row
  per residue — this is the per-residue result), `_detailed.csv`, and
  `<prefix>_spectral_density_plots.pdf`.
- **Column names differ between single- and dual-field runs.** Dual writes `R1_f1`/`R1_f2`,
  `Rex_f1_err`, `Rex_f1_95CI`. Single writes the bare names `R1`, `Rex_err`, `Rex_95CI` **and
  `_f1` aliases of each**, so one reader handles both shapes. `Residue, S2, tc, te, S2_err,
  te_err, chi2, fit_success` are never suffixed. CSVs written before the aliases existed have
  only the bare names, and a `_f1` selection on those returns an empty selection instead of
  raising — which reads as a real result ("no residue has Rex > 1"), not an error.
- **It re-fits T1/T2 internally** from the matrices, so its residue count is its own and will
  not match a `t1t2` run you filtered yourself. Run `t1t2` separately for per-residue T values.
- **Run each field alone alongside `--dual`.** Dual assumes one global τc; when the two
  datasets disagree the fit absorbs the difference as exchange, inflating Rex on most residues
  while median χ² rises by an order of magnitude. Comparing the Rex>1 fraction and χ² between
  `--dual` and the two single-field runs is what tells you whether the dual fit is reporting
  chemistry or a τc mismatch.
- Also accepts the density physical flags (`--rnh 1.015`, `--csa -172`) and fit-init overrides
  `--init-amp 5.0`, `--init-t1 800`, `--init-t2 100`, `--n-bootstrap 1000`,
  `--error-method {analytical,bootstrap}`, `--n-monte-carlo 50`.

### kd / export kd
```bash
python -m lunaNMR kd --input <titration.csv> --out <dir> --p0 <conc> [--prefix kd] \
    [--observable csp,intensity] [--alpha 0.14] [--intensity-from {height,volume}] \
    [--conc c0,c1,...] [--bootstrap 0] [--intensity-scale s0,s1,...]
python -m lunaNMR export kd --json <kd_fit.json> --out <dir> \
    [--fig-format pdf,png] [--per-page 20] [--kind curves,ref-bars,kd-bars,global-fit] \
    [--observable csp,intensity] [--summary-only] [--prefix DNAJA1_HSPA8]
```
CSP uses the full 1:1 quadratic isotherm; intensity uses exponential decay + plateau. All fits
bound Kd ≥ 0. `export kd` writes `summary.csv` plus, per `--kind`, all prefixed with `--prefix`
if given (default: none, i.e. unprefixed `summary.csv`/`<obs>_fits.pdf`/...):
- `curves` — `<obs>_fits.pdf` (multi-page grid of per-residue binding fits; `--fig-format png`
  writes one file per residue under `<obs>/` instead).
- `ref-bars` — `<obs>_ref_vs_point.pdf`, the observable per residue between point 0 and each
  later titration point (one page per point; **PDF only**). Matches the GUI viewer's export.
- `kd-bars` — `<obs>_kd_vs_residue.pdf`, per-residue Kd bars + the shared global-Kd line
  (**PDF only**). The intensity Kd is an apparent decay constant (see the global-fit caveat),
  not a thermodynamic dissociation constant.
- `global-fit` — `<obs>_global_fit.pdf`, per-residue observed data + the **single shared-Kd**
  global-model curve (one Kd for all residues, per-residue amplitudes), 20 panels/page, PDF
  only. Each panel's title carries R²(global) = the data vs the shared-Kd curve, so you can
  see which residues the one Kd fits poorly. Needs the run's `global[<obs>]` (≥2 residues).

`--kind` default `curves,ref-bars,kd-bars,global-fit`; `--fig-format` default `pdf`. Ref-bars
need the embedded per-point `series` + ≥2 points (both present in any `kd`-produced JSON).

### project
```bash
python -m lunaNMR project inventory <bundle.lunaNMR>
python -m lunaNMR project remove <bundle.lunaNMR> <rel/path> [<rel/path> ...]
```

### batch
```bash
python -m lunaNMR batch <folder> [--preset 1H] [--config ...] [--optimize] ...
```
Flags pass straight through to `batch_processing`.

## End-to-end example: dual-field model-free (600 + 800 MHz)

Data layout — one series folder per experiment, each with its own peak list:
```
600_WT/{600_T1,600_T2,600_hetNOE}/…ms.ft   +  <folder>/600_DNAJA1.txt
800_WT/{800_T1,800_T2,800_hetNOE}/…ms.ft   +  <folder>/800_DNAJA1.txt
```

**Stage 1 — process all six series** (repeat with each folder's own peak list):
```bash
python -m lunaNMR series --spectra 600_WT/600_T1 --peaks 600_WT/600_T1/600_DNAJA1.txt --out RUN/600/T1
# … 600 T2, 600 hetNOE, 800 T1, 800 T2, 800 hetNOE
```

**Stage 2 — hetNOE sat/unsat CSVs.** The hetNOE matrix columns are the two plane **filename
stems** (there is no delay to parse), so name or symlink the planes explicitly first — e.g.
`hetnoe_600_saturated.ft` / `hetnoe_600_unsaturated.ft`. Which plane is saturated is decided
by measuring the intensity ratio over high-S/N peaks (saturated is lower, typically 0.75–0.85),
not by a `001`/`002` filename. Then split the two columns into headerless
`residue,intensity[,error]` files (`noe_sat.csv`, `noe_unsat.csv`) per field. **Supply the
third column.** With no error column `calculate_hetnoe_from_intensities` falls back to a flat
`noe_err = noe * 0.02` — about 0.016 at a hetNOE of 0.78, roughly 3× tighter than a realistic
per-plane error (~0.044 relative per plane → hetNOE_err ≈ 0.05). Under-stated hetNOE errors
over-weight hetNOE in the model-free fit and make the S² errors look artificially tight.

**Stage 3 — model-free** (delays already `ms`, so no units flags):
```bash
python -m lunaNMR dynamixs modelfree --dual \
  --f1-t1 RUN/600/T1/peak_intensity_matrix.csv --f1-t2 RUN/600/T2/peak_intensity_matrix.csv \
  --f1-noe-sat RUN/600/noe/noe_sat.csv --f1-noe-unsat RUN/600/noe/noe_unsat.csv \
  --f2-t1 RUN/800/T1/peak_intensity_matrix.csv --f2-t2 RUN/800/T2/peak_intensity_matrix.csv \
  --f2-noe-sat RUN/800/noe/noe_sat.csv --f2-noe-unsat RUN/800/noe/noe_unsat.csv \
  --field1-freq 600 --field2-freq 800 --out RUN/modelfree_out
```

**Gotcha — the peak-detection window.** The default 1H/15N search window is **0.070 ppm** in
both dimensions. A wide 15N window ruins the low-S/N saturated hetNOE spectrum (peaks fit to
noise → Height 0 → NOE≈0 → residues dropped in model-free). Keep it tight.

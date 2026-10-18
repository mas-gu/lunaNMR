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
  is missing. Accepted by every subcommand that writes something. `diagnose` takes `--format`
  only: it is read-only, so a dry-run has nothing to mean. `project inventory` likewise takes
  `--format` only; `project export`, `remove` and `build` take both `--format` and `--dry-run`.
  `batch` takes `--format {text,json}` (under json its own text goes to stderr) and its own
  `--dry-run`, which lists the files and exits 1 when it finds none.

### Exit codes

| code | meaning |
|---|---|
| `0` | ran to completion |
| `1` | bad input (missing or unwritable path, malformed CSV/JSON, no parseable delay) — reported as `error: …` on stderr, no traceback; **or** a dry-run whose inputs are missing |
| `2` | no subcommand, or a command line the parser rejects: an unrecognised flag, a bad value, or a flag combination a subcommand refuses (argparse's usage report on stderr) |

Under `--format json` both `1` and `2` also print an `{"ok": false, …}` object on stdout; a
rejected command line reports `error.type` `UsageError`. `ok` is `true` exactly when the exit
code is `0`. A run that returns `1` after producing its summary — the findings below, or a
dry-run with a missing input — keeps the summary's keys and adds `ok: false` with an `error`
of type `MissingInput`, `NothingFitted`, `PreflightFailed` or (`batch`) `BatchFailed`. A
`series` that fits nothing reports `NothingFitted` with no summary keys; any other failure
names the exception's class (`ValueError`, `FileNotFoundError`, …), with an `error.details`
object where the failure has more to say (see `kd`), and one no input should
cause is marked `error.internal: true`, with its traceback on stderr. Called in-process,
`lunaNMR.cli.main()` returns these codes — including `0` for `--help`/`--version` — rather
than raising `SystemExit`.

### One call at a time per process

A run under `--format json` points fd 1 and `sys.stdout` at stderr for the whole process
while its engine runs, so two runs in one process cannot overlap safely. `cli.main()` therefore
holds a process-wide lock: a second thread that calls it waits until the first call returns.
A host that needs calls to run concurrently — an MCP server, for one — should run each call in
its own process (`python -m lunaNMR …` as a subprocess), not on threads of its own process.

Beyond that, some subcommands report *findings* through the exit code:

- **`diagnose`** returns `1` when any finding is `FAIL`, so `diagnose && series` stops there.
  `WARN` alone returns `0` — it covers routine cases like a folder with no peak list — unless
  you pass `--strict`. `series --dry-run --deep` runs the same checks and gates the same way.
- **`series`**, **`t1t2`**, **`methyl-t2`** return `1` when nothing fitted, **`t1rho`** when the
  R2 table comes out empty.
- **`hetnoe`**, **`density`**, **`modelfree`**, **`kd`**, **`export kd`** return `0` on any run
  that completed, even with residues dropped. Check `n_fitted`/`n_successful` against the
  matching total in the JSON summary; the exit code is not a partial-failure signal.

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
| `integrate-1d` | Measure peaks across a 1D spectrum series |
| `peaks shift` | Apply a rigid ppm offset to a peak list, measured or given |
| `project build` / `inventory` / `export` / `remove` | Build / inspect / extract / prune a `.lunaNMR` bundle |
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

**Which peak list a spectrum is checked against.** A folder with one `*.txt` uses it for
every spectrum. A folder with several pairs each spectrum with the list of the same stem
(`state.ft` ↔ `state.txt`), each pair reported as its own experiment named `folder:stem`,
since one folder often holds several states each with its own list. A spectrum with no list
of its name in such a folder is not checked, and a `WARN` names it and the lists that were
there — judging it against some other list reads a wrong list as decay.

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

### project export
```bash
python -m lunaNMR project export <bundle.lunaNMR> --out <folder> [--force] [--dry-run]
```
Copies a bundle's contents into a plain folder, so an analysis can be read without the
app. A bundle is already a directory of CSV/JSON, so nothing is converted — the point is
getting it out headlessly and completely.

**Everything is copied**, and anything `project inventory` does not classify is *named* in
the summary rather than absorbed: a file in no category is also a file `project remove`
cannot target, so it is worth seeing. Refuses a non-empty `--out` unless `--force`.

### project build
```bash
python -m lunaNMR project build <run_dir> --out <bundle.lunaNMR> [--spectrum <file>] [--peaks <file>]
```
Turns an automated run's output folders into a project the GUI opens. Under `<run_dir>` it
takes, wherever they sit - symlinked folders included, each folder read once however many
links reach it, so a link back into the run does not loop:

- every series (`batch_results.json`, written by `series`);
- every Kd fit (`<prefix>_kd_fit_data.json`, written by `kd`). A Kd fit whose prefix
  matches a series' folder name is linked to that series' tidy CSV;
- every DynamiXs run, from the `<prefix>_dynamixs.json` record each `dynamixs` command
  writes beside its outputs (command, absolute inputs, parameters, summary). All
  `dynamixs t1t2` fits become ONE T1/T2 analysis - the GUI page holds both fields, field 1
  the lower frequency. A `dynamixs t1rho` run joins that analysis: it fills its field's
  T1rho slot with the field's tilt, spin-lock field and carrier, and its T1 fit and
  derived T2 fill their slots only where no `t1t2` fit did. Each `modelfree`,
  `methyl-t2`, `hetnoe` and `density` run becomes its own analysis, on the page of that
  name. Runs made before the record existed carry none and are not found.

It saves through `ProjectManager.save_project` itself, handed a stand-in main window, so
the bundle is exactly what the GUI would have saved. The GUI reopens the first series'
reference spectrum and the peak list it was fitted against unless `--spectrum`/`--peaks`
say otherwise. Both are **referenced by absolute path, not copied**: the bundle needs the
data where it was. With no `--spectrum`, a reference spectrum that no longer exists where
the series read it is refused (exit 1) with a message to pass `--spectrum`; so is a given
`--spectrum` that does not exist. DynamiXs results are copied into the bundle, so they reopen after the
run folder is gone. Refuses an existing `--out`, and leaves nothing behind when a build
fails. `--dry-run` performs the whole build into a scratch folder and reports what the
bundle would hold.

A T1/T2 experiment fitted twice at one field is refused rather than one silently replacing
the other: the page holds one. Two model-free runs of one `--prefix` are both bundled,
the second named after where its record sits. Every DynamiXs analysis keeps the CLI runs
it came from (`analysis_metadata.cli_runs`), so inputs the page has no field for - an
`--f*-r2-table` - stay on record.

Runs that share an `--out` folder are bundled apart: each analysis gets only the files its
own run wrote - every `dynamixs` output starts with `<prefix>_`, and a file belongs to the
longest prefix of that folder's runs it starts with, so `a` does not claim `a_b_*`. A
folder shared by two methyl-t2 runs otherwise gave each analysis both fit JSONs, and its
viewer, which opens every fit in its folder, showed the other run's fits as its own.

A series is stored slim, by `series` and by the GUI alike: each fit keeps its parameters
(`all_peaks`) and its region's ppm axes, and the viewers rebuild the region and surfaces
from the spectrum when a fit is opened - bit-identical to storing them. ~3 MB for a
7-point titration, against ~80 MB whole.

### integrate-1d
```bash
python -m lunaNMR integrate-1d --spectra <folder> --out <dir> \
    (--peaks 5.00,5.30 | --detect [--min-snr 10]) [--value area] [--prefix integration]
```
Measures each peak in every spectrum of a 1D series. Peaks are re-matched at each point
within a drift window bounded by the peak's own linewidth, and matching is competitive —
two peaks can never collapse onto one maximum, the failure that silently turns an
intensity ratio into 1.0.

Writes a wide `<prefix>_series.csv` (one row per spectrum, one column per peak) and a long
`<prefix>_table.csv` carrying every observable per measurement.

**Quantify with `area`, the region sum.** On the reference 53-point series a 1:1 conversion
conserves to **6.4% for the region sum, 11.5% for height and 38% for the fitted area**.
Height fails because a broadening peak loses height at constant area. The fitted area fails
because a close neighbour caps the window at ~1.4 linewidths, leaving the σ/γ split
unconstrained while the analytic area integrates tails the data never sampled — and
**R² stays above 0.97 throughout without validating that extrapolation**, because it only
measures the fit inside the window. Use `fit_area`, `r_squared` and `fwhm` diagnostically.

The summary reports `quantitative_observable` and `diagnostic_observables` so a caller
does not have to know this, and counts measurements where no maximum was found (those are
read at the expected position rather than allowed to climb a neighbour's flank).

### peaks shift
```bash
python -m lunaNMR peaks shift --peaks <list.txt> --out <shifted.txt> \
    [--dx 0.0] [--dy 0.0] [--auto --spectrum <ref.ft>] [--quick]
```
`diagnose` measures a peak list's registration offset and stops there. This applies it.
The error is usually rigid — one list referenced against a slightly different frame, not a
wrong list — so a single offset corrects the whole file, and without this the only headless
remedy was to re-pick every peak.

`--auto` measures the offset against `--spectrum` with the same routine `diagnose` uses,
then applies it; `--dx`/`--dy` give it explicitly. `--out` must differ from `--peaks`: the
input is often the only copy of an assignment nobody wants to redo.

**A measured offset is checked before it is written.** The search grid spans ±0.09 ppm in
¹H, and it always has a maximum: a list further out than that gets a spurious one inside
the grid, and applying it leaves the list off in both dimensions. So `--auto` re-reads the
spectrum at the shifted positions and refuses (exit 1, `ok: false`, nothing written) when
the median S/N there is below 10 or fewer than 85% of the peaks are captured — the
thresholds `diagnose` uses. On the real reference spectrum with the list 0.15 ppm out in
¹H, the "correction" left a median S/N of 3 and 6/68 peaks captured; a correct one gives
~124 and 64/68. A written `--auto` summary carries `captured` and `median_snr`. Given
`--dx`/`--dy` are applied unchecked.

**`--quick` is a resolution limit, not just a speed knob.** It quantises the ¹⁵N offset to
0.03 ppm instead of 0.015, so a smaller real offset is misreported — as zero, or as a full
step, i.e. twice the truth — and either reads as a confident measurement. Measure on the
full grid before correcting anything.

### series
```bash
python -m lunaNMR series --spectra <folder> --peaks <list.txt> --out <dir> \
    [--mode {time,titration}] [--peak-source {reference,cascade,detected,independent}] [--parallel]
```
- Peak list: `Assignment, Position_X, Position_Y[, Height]` header. Every assignment must be
  unique: the series matches peaks by assignment, so a duplicated label is refused, naming it.
- `--mode time` (relaxation, default) or `titration` (forces cascade tracking with no
  absolute drift limit, and gives each peak its own position window — the shared one, shrunk
  to 0.49× the gap to its nearest neighbour — so no peak can reach another's). The summary reports the
  `peak_source` that ran beside the `requested_peak_source`.
  `--peak-source independent` is refused with `titration`: it detects every point against the
  reference list, which anchors each point to the free state instead of tracking it.
- Filename values: `time` mode reads a trailing delay with an explicit unit (`_300ms`, `_2.4s`,
  `_2400us`). `titration` mode reads a trailing `_<value>` whose value carries a **decimal
  separator** — `o` or `.` (`_0o0`, `_0o5`, `_2.4`). Digits alone (`sample_2.ft`, `ref_100.ft`)
  parse in neither mode; those spectra get a stem-named column instead of a value (and `kd`
  reads that stem's trailing digits as a titration point, so rename them first).
- Spectra run in order of their parsed value — the order a cascade seeds each point from
  the one before — and spectra with no value run last, by name.
- `--peak-source reference` (default) holds list positions fixed across the series — correct for
  T1/T2/hetNOE.
- `--params <file.json>` overrides the detection/fitting knobs this page tells you to tune,
  deep-merged over the defaults so naming one key leaves the rest alone:
  ```json
  {"detection_params": {"search_window_x": 0.05},
   "gui_params":       {"fix_linewidths": true}}
  ```
  **An unknown key is refused, with a suggestion** — a knob that is accepted, does nothing and
  still reports success is the failure mode this whole surface exists to prevent, and it is
  checked during `--dry-run` too. Command-line flags (`--parallel`, `--peak-source`) win over
  the file, so a stale tuning file cannot silently override what you just typed. The run
  summary - and the `--dry-run` summary - records the effective `params`. A key the defaults carry but nothing on the series
  path reads is refused as having no effect: `gui_params.use_centroid_refinement`,
  `centroid_noise_multiplier`, `use_ps2d_multi_peak`, `height_threshold`, `distance_factor`,
  `prominence_threshold`, `smoothing_sigma`, `max_peaks_fit`, `max_optimization_iterations`,
  `fitting_params.min_r_squared`, `processing_options.use_global_optimization`,
  `detection_params.noise_threshold` (series detection runs at its own S/N threshold), and
  `use_parallel_processing` in either section (`--parallel` sets it). With `--mode titration`,
  `processing_options.use_original_reference_for_detection: true` is refused too.
- `"gui_params": {"use_adaptive_optimization": true}` turns on the fitting-window grid
  search (1.0 / 1.5 / 2.0 x the median linewidth). It is **off by default**: on the DNAJA1
  and ERDJ6 titrations it chose 1.0x both times, a window that cuts the peak at half height,
  and ERDJ6's CSP Kd went from 282 µM to pinned at its 1500 µM bound.
- `"gui_params": {"auto_add_dummy_peaks": false}` stops the first spectrum's detection from
  adding the unmatched peaks it finds near listed ones as `dummy_NNN` peaks. It is **on by
  default**; those dummies are fitted at every point and appear in the matrices and the tidy
  CSV. The series dialog takes it from the Expert "Auto-add nearby unmatched peaks" box.
- Writes `peak_intensity_matrix.csv`, `_volume_matrix`, `comprehensive_peak_tracking.csv`,
  `series_analysis_tidy.csv`, `series_metadata.json`, `batch_results.json` (the series as
  the GUI holds it, read by `project build`), per-spectrum CSVs. A spectrum that produced no
  fits is named in the summary's `failed_spectra` (its matrix columns are empty), and one
  some of whose peaks a fallback algorithm fitted in `fitting_fallbacks`.
  `peak_detection_matrix.csv` (and each `{pt}_Detected` tracking column) is 1 where the fit
  stands at least 3 sigma above the noise (`detection_sigma`, the Kd loader's floor), else 0.

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
lists spectrum files, not columns. This matters as soon as you need to treat one sub-series
differently from another (cross-normalising a re-acquired tail, say). The same file records
how the series fitted, so one point can be refitted the same way: `ps2d_config`, and per
column, where that point had them, the `territories` and `linewidth_caps` it was fitted with.

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

The same file carries `detection_sigma` (height ÷ noise), `position_pinned` and
`intensity_pinned` (the fit stopped on a bound). Neither `quality` nor `r_squared` can see
those states. A decayed peak (`detection_sigma` < 3) keeps its intensity — zero is the
measurement in a decay — but has no position for a CSP. A pinned position has no position
for a CSP either, but keeps its intensity. A fit on its intensity floor is no measurement at
all, so `series` writes it as an empty cell in `peak_intensity_matrix.csv` and
`peak_volume_matrix.csv`, and the relaxation fitters skip empty cells. `kd` applies the full
rule when its input is the tidy CSV; on `comprehensive_peak_tracking.csv` it applies both
bound flags but cannot see undetected points. The table and the reasons are in
[CLI_AGENT.md](CLI_AGENT.md) under Gotchas.

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
marker `b` or `_2`), which returns every delay in ms. Exit code 1 if nothing fits.
Writes `<prefix>_fit_results.txt` (`t1t2` adds its `_fig{1,2,3}.pdf` plots), the fit JSON
`<prefix>_<field-name>_<exp>_fit_data.json` (`methyl-t2`: `<prefix>_<field-name>_methylT2_fit_data.json`),
and the `<prefix>_dynamixs.json` record every `dynamixs` subcommand leaves for `project build`.
A series sidecar with `n_value_unparsed > 0` is refused, naming the spectra, unless
`--allow-unparsed-delays` (also on `t1rho` and `modelfree`).

`--time-units` labels the output (summary, fit JSON, plot axes); it does not rescale. Left
out, it is taken from `series_metadata.json` beside the input, else `ms` — the unit the
parser returns and `series` writes. On headers that state a unit, any value but `ms` is
refused: those delays are already in ms.

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
J(0)/J(ωN)/J(ωH), S², τc, τe, Rex in `<prefix>_results.csv` (`--dual`: `_basic.csv` and
`_detailed.csv`) plus `<prefix>_plots.pdf`. S², τe and Rex are fitted by bounded least
squares from a grid of starts, which carries no uncertainty of its own: their errors are
NaN unless `--monte-carlo`.

`--dual` needs `--input2` and `--field2-freq`. `--field2-freq` has no default - there is no
conventional second field - and its absence is a usage error (exit 2), as on `modelfree`.

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
- **Per-series time units:** the T→R conversion needs each series' delay units. Headers that
  state a unit (`..._300msb`, `..._2o4s`) need no flag: the parser reads them as ms whatever
  unit they state, and any flag but `ms` is refused for them — `s` on `..._2o4s` headers would
  convert them a second time. Neither does a matrix produced by `series`, which is already
  normalised to ms (`series_metadata.json` records `value_units`; a flag contradicting it is
  refused too). The flag is for bare-number headers `series` did not write: a hand-built or
  DynamiXs-format T1 in seconds (`0.1, 0.4, …`) beside a T2 in ms → `--f1-t1-units s`, without
  which R1 is out by 1000×.
- **Outputs:** `<prefix>_field{1,2}_T{1,2}_fit_results.txt` (+ `_fig{1,2,3}.pdf`),
  `<prefix>_field{1,2}_input.csv` (the R1/R2/hetNOE table the density step read),
  `<prefix>_field{1,2}_hetnoe.pdf`, `<prefix>_spectral_density_basic.csv` (one row
  per residue — this is the per-residue result), `_detailed.csv`, and
  `<prefix>_spectral_density_plots.pdf`.
- **Column names differ between single- and dual-field runs.** Dual writes `R1_f1`/`R1_f2`,
  `Rex_f1_err`, `Rex_f1_95CI`. Single writes the bare names `R1`, `Rex_err`, `Rex_95CI` **and
  `_f1` aliases of each**, so one reader handles both shapes. `Residue, S2, tc, te, S2_err,
  te_err, chi2, fit_success` are never suffixed. CSVs written before the aliases existed have
  only the bare names, and a `_f1` selection on those returns an empty selection instead of
  raising — which reads as a real result ("no residue has Rex > 1"), not an error.
- **The model-free fit is bounded least squares from a grid of starts** (S², τe and Rex,
  including Rex > 0 starts, at fixed τc), deterministic for the same input. The model's
  hetNOE is predicted from the model's own T1, not the measured one. The fit itself carries
  no uncertainty; the errors come from the `--n-monte-carlo` refits.
- **`n_successful` counts a fit that ran out of evaluations.** A fit whose best start stopped
  on its evaluation cap rather than a convergence test still reports `fit_success`, so the
  summary carries `n_at_iteration_cap` and `at_iteration_cap_residues` beside it, and the
  per-residue CSV an `at_iteration_cap` column. On a real 600/700 MHz dataset none did,
  single- or dual-field.
- **It re-fits T1/T2 internally** from the matrices, so its residue count is its own and will
  not match a `t1t2` run you filtered yourself. Run `t1t2` separately for per-residue T values.
- **Run each field alone alongside `--dual`.** Dual assumes one global τc; when the two
  datasets disagree the fit absorbs the difference as exchange, inflating Rex on most residues
  while median χ² rises by an order of magnitude. Comparing the Rex>1 fraction and χ² between
  `--dual` and the two single-field runs is what tells you whether the dual fit is reporting
  chemistry or a τc mismatch.
- Also accepts the density physical flags (`--rnh 1.015`, `--csa -172`) and fit-init overrides
  `--init-amp` (default: derived from the data), `--init-t1 800`, `--init-t2 100`, `--n-bootstrap 1000`,
  `--error-method {analytical,bootstrap}`, `--n-monte-carlo 50`.

### kd / export kd
```bash
python -m lunaNMR kd --input <titration.csv> --out <dir> --p0 <conc> [--prefix kd] \
    [--conc c0,c1,...] [--conc-units {absolute,equivalents}] \
    [--observable csp,intensity] [--alpha 0.14] [--intensity-from {height,volume}] \
    [--bootstrap 0] [--intensity-scale s0,s1,...] [--survey | --residues <file|list>]
python -m lunaNMR export kd --json <out>/data/<prefix>_kd_fit_data.json --out <dir> \
    [--fig-format pdf,png] [--per-page 20] [--kind curves,ref-bars,kd-bars,global-fit] \
    [--observable csp,intensity] [--summary-only] [--prefix DNAJA1_HSPA8]
```
CSP uses the full 1:1 quadratic isotherm; intensity uses exponential decay + plateau. All fits
bound Kd ≥ 0. `kd` writes `<prefix>_kd_results.txt` in `--out` and the fit JSON and
`<prefix>_kd_params.json` under `<out>/data/`. `--survey` fits nothing global: it writes the
vs-sequence figures, `<out>/data/<prefix>_survey.csv`, and an editable
`<prefix>_residues.txt` beside the input for `--residues`. The thresholds, the survey →
fit workflow and the pool gates are in [CLI_AGENT.md](CLI_AGENT.md#kd-workflow-survey--choose--fit).

**Concentrations are required, with their units.** `--conc` is read in `--conc-units`
(default `absolute`; `equivalents` are multiplied by `--p0`). Without `--conc` they come from
a params file stored for this same input and titration points (the survey's, beside the
input, or `kd_analysis/` beside it); a `--conc-units` contradicting it is refused. With
neither, `--conc-units` alone says the spectrum names are the concentrations. With none of
these the run — survey included — is refused rather than fitting point labels as µM. The
summary's `params_source` names the file read and `params_skipped` every file passed over,
with the reason; a run refused for want of concentrations names them in its message and,
under `--format json`, in `error.details.params_skipped`. Only concentrations, intensity scales and the survey thresholds are restored;
`--p0`, `--alpha`, `--intensity-from` and `--observable` are taken from the command line.

`export kd` writes PDFs in `--out` and its CSV/JSON under `<out>/data/`: `summary.csv` plus,
per `--kind`, all prefixed with `--prefix` if given (default: none, i.e. unprefixed
`summary.csv`/`<obs>_fits.pdf`/...):
- `curves` — `<obs>_fits.pdf` (multi-page grid of per-residue binding fits; `--fig-format png`
  writes one file per residue under `<obs>/` instead).
- `ref-bars` — `<obs>_ref_vs_point.pdf`, the observable per residue between point 0 and each
  later titration point (one page per point). Also writes `<obs>_ref_vs_point.csv` (wide,
  Excel-ready) and `.json`. Matches the GUI viewer's export.
- `kd-bars` — `<obs>_kd_vs_residue.pdf`, per-residue Kd bars for the residues in the shared
  fit's pool, titled "N of M residues", + the shared global-Kd line (**PDF only**). Matches
  the GUI viewer's export. The intensity Kd is an apparent decay constant (see the global-fit caveat),
  not a thermodynamic dissociation constant.
- `global-fit` — `<obs>_global_fit.pdf`, per-residue observed data + the **single shared-Kd**
  global-model curve (one Kd for all residues, per-residue amplitudes), 20 panels/page. Also
  writes `<obs>_global_fit.csv` and `.json` with the per-residue parameters behind the global
  fit. Each panel's title carries R²(global) = the data vs the shared-Kd curve, so you can
  see which residues the one Kd fits poorly. Needs the run's `global[<obs>]` (≥2 residues).

`--kind` default `curves,ref-bars,kd-bars,global-fit`; `--fig-format` default `pdf`. Ref-bars
need the embedded per-point `series` + ≥2 points (both present in any `kd`-produced JSON).

### project
```bash
python -m lunaNMR project build <run_dir> --out <bundle.lunaNMR>
python -m lunaNMR project inventory <bundle.lunaNMR>
python -m lunaNMR project remove <bundle.lunaNMR> <rel/path> [<rel/path> ...] [--dry-run]
```
`remove` refuses a folder with no `project.json`, and refuses `project.json` itself (shown
`(protected)` by `inventory`); a refused path deletes nothing in that call. The summary
separates `removed` from `missing` (named paths that were not there).

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

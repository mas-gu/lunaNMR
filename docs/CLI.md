<!-- ABOUTME: Command reference for the headless lunaNMR CLI (python -m lunaNMR). -->
<!-- ABOUTME: Subcommand flags plus the end-to-end dual-field model-free workflow. -->

# lunaNMR CLI Reference

Headless access to the analysis pipeline. Every reasonable GUI analysis is scriptable — no
display required.

```bash
python -m lunaNMR <subcommand> [flags]      # or the console script:  lunanmr <subcommand>
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
| `series` | Process a multi-spectrum series/titration into intensity/position matrices |
| `dynamixs t1t2` | Mono-exponential T1 or T2 relaxation fit |
| `dynamixs methyl-t2` | Bi-exponential (Tugarinov-Kay) methyl T2 fit |
| `dynamixs hetnoe` | Heteronuclear NOE = I_sat / I_unsat per residue (+ QC plot) |
| `dynamixs density` | Reduced spectral density mapping (single or dual field) |
| `dynamixs modelfree` | Integrated pipeline: T1/T2 fit → hetNOE → density → Lipari-Szabo |
| `kd` | Kd titration (CSP quadratic / intensity decay) |
| `export kd` | CSP/intensity figures + summary from a saved Kd fit JSON |
| `project inventory` / `remove` | Inspect / prune a `.lunaNMR` bundle |
| `batch` | Folder-wide peak detect + Voigt/PS2D fit |

### series
```bash
python -m lunaNMR series --spectra <folder> --peaks <list.txt> --out <dir> \
    [--mode {time,titration}] [--peak-source {reference,cascade,detected,independent}] [--parallel]
```
- Peak list: `Assignment, Position_X, Position_Y[, Height]` header.
- `--mode time` (relaxation, default) or `titration` (forces cascade tracking + wider margins).
- `--peak-source reference` (default) holds list positions fixed across the series — correct for
  T1/T2/hetNOE.
- Writes `peak_intensity_matrix.csv`, `_volume_matrix`, `comprehensive_peak_tracking.csv`,
  `series_analysis_tidy.csv`, per-spectrum CSVs.

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
  need no flag. Outputs `<prefix>_field{1,2}_hetnoe.pdf`, the density basic/detailed CSVs, and
  `<prefix>_spectral_density_plots.pdf`.
- Also accepts the density physical flags (`--rnh 1.015`, `--csa -172`) and fit-init overrides
  `--init-amp 5.0`, `--init-t1 800`, `--init-t2 100`, `--n-bootstrap 1000`,
  `--error-method {analytical,bootstrap}`, `--n-monte-carlo 50`.

### kd / export kd
```bash
python -m lunaNMR kd --input <titration.csv> --out <dir> --p0 <conc> [--prefix kd] \
    [--observable csp,intensity] [--alpha 0.14] [--intensity-from {height,volume}] \
    [--conc c0,c1,...] [--bootstrap 0] [--intensity-scale s0,s1,...]
python -m lunaNMR export kd --json <kd_fit.json> --out <dir> \
    [--fig-format pdf,png] [--per-page 20] [--observable csp,intensity] [--summary-only]
```
CSP uses the full 1:1 quadratic isotherm; intensity uses exponential decay + plateau. All fits
bound Kd ≥ 0. `export kd` writes `summary.csv`; `--fig-format` default `pdf` (multi-page grid
per observable), `png` = one file per residue.

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

**Stage 2 — hetNOE sat/unsat CSVs.** The hetNOE `peak_intensity_matrix.csv` has `…_sat` and
`…_unsat` columns; split each into two headerless `residue,intensity` files (`noe_sat.csv`,
`noe_unsat.csv`) per field.

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

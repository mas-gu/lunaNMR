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
| `series` | Fit a multi-spectrum series/titration | `--spectra <dir/glob>` `--peaks` `--out` | `peak_intensity_matrix.csv`, `_volume_matrix`, `comprehensive_peak_tracking.csv`, `series_analysis_tidy.csv`, per-spectrum CSVs |
| `dynamixs t1t2` | Mono-exp T1/T2 fit | `--input <matrix.csv>` `--out` `--exp {T1,T2}` | `<prefix>_fit_results.txt` + JSON |
| `dynamixs methyl-t2` | Bi-exp Tugarinov-Kay methyl T2 | `--input` `--out` | `<prefix>_fit_results.txt` + JSON |
| `dynamixs hetnoe` | I_sat/I_unsat per residue | `--sat` `--unsat` `--out` | `<prefix>_hetnoe.csv`, `<prefix>_hetnoe.pdf` |
| `dynamixs density` | Reduced spectral density mapping | `--input <table.csv>` `--out` `--field1-freq` | `<prefix>_results.csv` (or `_basic`/`_detailed` for `--dual`) + `_plots.pdf` |
| `dynamixs modelfree` | T1/T2→hetNOE→density→Lipari-Szabo | `--f1-t1 --f1-t2 --f1-noe-sat --f1-noe-unsat --field1-freq --out` | model-free CSVs, per-field hetNOE.pdf, density plots |
| `kd` | Kd titration (CSP quadratic / intensity decay) | `--input <titration.csv>` `--out` `--p0 <conc>` | `<prefix>_kd_fit_data.json`, results.txt |
| `export kd` | Figures + summary from a kd fit JSON | `--json` `--out` | `summary.csv`, `<obs>_fits.pdf` / per-residue PNGs |
| `project inventory` | List a bundle's contents | `<bundle.lunaNMR>` | stdout listing |
| `project remove` | Delete bundle-relative paths | `<bundle> <rel/path>…` | (mutates bundle) |
| `batch` | Folder-wide detect + Voigt/PS2D fit | `<folder>` (flags pass through to `batch_processing`) | per-spectrum fit outputs |

## Copy-paste blocks
```bash
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
- **`--parallel`** → `series` only (two-pass processor, ~2.7×). `density`/`modelfree` parallelize internally and accept `--no-parallel`.
- **Exit codes:** `2` = no/invalid subcommand (help to stderr). `1` = bad input (missing file, malformed CSV/JSON, unparseable delay label — reported as `error: …` on stderr, no traceback) OR dry-run with a missing input. `series`/`t1t2`/`methyl-t2` also return `1` when nothing fit. `hetnoe`/`density`/`modelfree`/`kd`/`export` return `0` on any successful run even if some residues failed — check the JSON `n_successful`/`n_fitted`, don't rely on exit code for partial failure.

## Input file shapes
- **Peak list** (`series --peaks`): header `Assignment, Position_X, Position_Y[, Height]` (comma or tab; column aliases standardized).
- **Series spectra** (`series --spectra`): folder or glob; files auto-discovered by extension (`ft ser ft2 ft3 pipe ucsf`), naturally sorted.
- **Relaxation matrix** (`t1t2/methyl-t2 --input`): a series `peak_intensity_matrix.csv` / `comprehensive_peak_tracking.csv`; delay columns parsed by `parse_delay_column` (explicit `ms`/`s`/`us` + optional repeat marker `b` or `_2`, e.g. `..._300msb`).
- **hetNOE sat/unsat** (`--sat`/`--unsat`, modelfree `--f{1,2}-noe-*`): headerless `residue,intensity[,error]`. Split from a hetNOE series matrix's `…_sat`/`…_unsat` columns.
- **Density table** (`density --input`): `Residue,R1,R1err,R2,R2err,hetNOE,hetNOEerr`.
- **Kd input** (`kd --input`): `series_analysis_tidy.csv`, `comprehensive_peak_tracking.csv`, or an intensity matrix (intensity-only, no CSP). Ligand conc from CSV point labels unless `--conc` overrides.
- **Kd fit JSON** (`export kd --json`): the self-contained `…_kd_fit_data.json` from `kd` (embeds per-point series + `metadata.protein_conc`).

## Gotchas (silent-corruption risks)
- **modelfree units rescale, t1t2 units don't.** `modelfree --f{1,2}-t{1,2}-units {ms,s,us}` (default `ms`) DRIVE the T→R conversion — a T1 in seconds with a T2 in ms needs `--f1-t1-units s`, else R1 is off by 1000×. Contrast `t1t2 --time-units` (default `s`) which only *labels* output.
- **CSP needs a fresh series run.** Older series output rounded positions to 0.1 ppm; CSP needs 4-decimal precision. Re-run `series` before `kd` CSP. Intensity fits work on existing data.
- **Kd is always bounded ≥ 0** internally (unbounded lands on a degenerate negative-Kd minimum). Don't post-process assuming otherwise.
- **Non-finite fit values serialize to JSON `null`.** A degenerate Kd fit (e.g. only 3 points) is `success=true` with finite `Kd` but `Kd_err`/`r_squared` = `null`. Treat `Kd`/`Kd_err`/`r_squared` as possibly `None` — legitimate, not a failure.
- **Search-window default is 0.070 ppm** (1H and 15N). Don't widen 15N for low-S/N spectra (e.g. saturated hetNOE) — wide windows fit peaks to noise → Height 0 → NOE≈0 → residues silently dropped.
- **`--method` field prefix follows `--dual`.** In `modelfree` field count derives from `--dual`; you only pick the `087` vs `jwh` variant (default `087`). A mismatched `single_*`/`dual_*` prefix is ignored.
- **`--intensity-scale` scales height/volume only** (per-point scan-count correction), never positions/CSP; a uniform scale cancels in the I/I₀ ratio.
</content>
</invoke>

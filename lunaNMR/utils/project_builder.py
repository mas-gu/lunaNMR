# ABOUTME: Builds a .lunaNMR bundle from an automated run's output folders, without a GUI.
# ABOUTME: Drives ProjectManager.save_project with a stand-in main window holding the run's results.

"""
An automated run writes plain folders; the GUI opens .lunaNMR bundles. Rather than
write a bundle's files a second way, this hands ProjectManager.save_project - the code
the GUI saves with - a stand-in main window carrying what the run produced. The bundle is
then exactly what the GUI would have saved had it held those results.

What a run contributes, found anywhere under its folder:
  series      every `batch_results.json` the `series` command writes
  Kd          every `<prefix>_kd_fit_data.json` the `kd` command writes
  DynamiXs    every `<prefix>_dynamixs.json` a `dynamixs` command writes: the `t1t2`
              fits become one T1/T2 analysis, each `modelfree` and `methyl-t2` run
              its own

The spectrum and peak list the GUI reopens are the first series' reference spectrum and
the peak list it was fitted against, unless given. Both are referenced by absolute path,
never copied: a bundle moved to another machine needs the data at the same place.
"""

import fnmatch
import json
import os
import shutil
import tempfile
import types
from pathlib import Path

import pandas as pd

from lunaNMR.utils.project_manager import ProjectManager, batch_results_from_dict

KD_FIT_SUFFIX = '_kd_fit_data.json'


def _files_under(run_dir, pattern):
    """Every file under `run_dir` whose name matches `pattern`, sorted.

    A run folder is often assembled by linking result folders in, so symlinked folders
    are followed. Each folder is entered once, by its resolved path: a link back into
    the run would otherwise be walked forever, and a folder linked in twice would give
    every result in it twice.
    """
    found, visited = [], set()
    for folder, subfolders, files in os.walk(run_dir, followlinks=True):
        real = os.path.realpath(folder)
        if real in visited:
            subfolders[:] = []
            continue
        visited.add(real)
        found.extend(Path(folder) / f for f in fnmatch.filter(files, pattern))
    return sorted(found)


def find_series(run_dir):
    """{name: batch_results.json path}, named after the folder each series wrote to."""
    found = {}
    for path in _files_under(run_dir, 'batch_results.json'):
        name = path.parent.name
        if name in found:   # two series folders of one name: tell them apart by place
            name = '_'.join(path.parent.relative_to(run_dir).parts)
        # A folder directly under run_dir has no place beyond its name, and a place
        # can spell another folder's name: number what is still taken.
        base, n = name, 2
        while name in found:
            name, n = f'{base}_{n}', n + 1
        found[name] = path
    return found


def find_kd_fits(run_dir):
    """{name: fit JSON path}, named after the `--prefix` the fit was written with."""
    found = {}
    for path in _files_under(run_dir, f'*{KD_FIT_SUFFIX}'):
        name = path.name[:-len(KD_FIT_SUFFIX)] or path.parent.name
        if name in found:
            name = '_'.join(path.parent.relative_to(run_dir).parts) + '_' + name
        found[name] = path
    return found


def kd_page_state(fit_data, fit_path, input_file=None):
    """The Kd page's saved state (kd_titration_page.get_session_state) for a CLI fit.

    The fit JSON's metadata carries every parameter the page restores. The input CSV is
    the one thing it does not record; without it the page rebuilds the input from the
    per-point series the JSON embeds.

    Needs dynamiXs_Kd on sys.path, as every Kd consumer does (cli._add_modules_path).
    """
    from kd_params import observables_to_combo_index
    meta = fit_data.get('metadata', {})
    points = list(meta.get('points') or [])
    scales = meta.get('intensity_scales') or [1.0] * len(points)
    return {
        'input_file': str(input_file) if input_file else None,
        'input_basename': Path(input_file).name if input_file else None,
        'output_dir': str(Path(fit_path).parent.parent),
        'detected_points': points,
        'concentrations': list(meta.get('concentrations') or []),
        'conc_units': meta.get('conc_units'),
        'intensity_scales': list(scales),
        'p0': meta.get('protein_conc'),
        'alpha': meta.get('alpha'),
        'observable': observables_to_combo_index(meta.get('observables') or []),
        'intensity_value': meta.get('intensity_value', 'height'),
        'bootstrap_iter': meta.get('n_bootstrap', 0),
        'last_json_file': str(fit_path),
        'last_json_folder': str(Path(fit_path).parent),
    }


DYNAMIXS_RECORD_SUFFIX = '_dynamixs.json'

# DynamiXs dialog stack indices (dynamixs_dialog: menu, t1t2, spectral, plotting,
# integrated, methyl_t2).
T1T2_PAGE, SPECTRAL_PAGE, INTEGRATED_PAGE, METHYL_T2_PAGE, HETNOE_PAGE = 1, 2, 4, 5, 6


def find_dynamixs_records(run_dir):
    """Every `<prefix>_dynamixs.json` a dynamixs command left under the run, each with the
    file it was read from as `record_file`."""
    records = []
    for path in _files_under(run_dir, f'*{DYNAMIXS_RECORD_SUFFIX}'):
        record = json.loads(path.read_text())
        record['record_file'] = str(path)
        records.append(record)
    return records


def _analysis_state(name, page_index, analysis_type, page_key, page_state, records):
    """A DynamiXs analysis state, with the CLI runs it came from in its metadata - every
    input on record, including those the page has no field for (an R2 table)."""
    runs = [{key: r[key] for key in ('command', 'args', 'record_file')} for r in records]
    return {'active_page': page_index, 'analysis_name': name,
            'analysis_metadata': {'analysis_name': name, 'analysis_type': analysis_type,
                                  'cli_runs': runs},
            page_key: page_state}


def t1t2_analysis(records, name, staging, t1rho_records=()):
    """ONE T1/T2 page analysis from a run's `dynamixs t1t2` and `dynamixs t1rho` records.

    The GUI page holds both fields and every experiment together; the CLI fits each one
    separately, into its own folder and often on the same default --prefix. Field 1 is
    the lower frequency. Each experiment's fit JSON goes into `json_<field>_<exp>/` under
    `staging` and its results file beside them, both under names starting with that
    key - the experiment's file prefix - so the project save copies one folder per page
    as it does for the GUI, and no experiment can land on another's file.

    A `t1rho` run fills its field's T1rho slot and brings the field's acquisition
    parameters (tilt, spin-lock field, carrier) and peak list. It also fitted a T1 and
    derived a T2; each fills its slot only where no `t1t2` fit did - a fit outranks a
    refit, and a measured T2 a derived one.

    :raises ValueError: for more than two fields - the page holds two - or one
        experiment fitted twice, of which the page can show only one
    """
    from lunaNMR.utils.t1rho_calculator import write_t2_fit_results
    all_records = list(records) + list(t1rho_records)
    freqs = sorted({float(r['args']['field_freq']) for r in all_records})
    if len(freqs) > 2:
        raise ValueError(f"{len(freqs)} fields of T1/T2 fits ({freqs}); the T1/T2 page "
                         f"holds two")
    field_of = {freq: f'field{i + 1}' for i, freq in enumerate(freqs)}
    staging = Path(staging)
    staging.mkdir(parents=True, exist_ok=True)

    session = {field_of[f]: {'frequency_mhz': f, 't1': None, 't2': None, 't1rho': None}
               for f in freqs}
    session.setdefault('field2', None)
    page = {'session_results': session, 'fitted_experiments': [],
            'field2_enabled': len(freqs) == 2,
            'field1_freq': freqs[0], 'output_dir': str(staging), 'prefix': name}
    if len(freqs) == 2:
        page['field2_freq'] = freqs[1]
    refs, source = {}, {}

    def fill(record, field, exp, results_file, json_file, input_file):
        key = f'{field}_{exp}'
        if key in source:
            raise ValueError(f"{exp.upper()} for {field} came from two runs - "
                             f"{source[key]} and {record['record_file']}; the T1/T2 page "
                             f"holds one. Build from a folder holding the one to keep.")
        source[key] = record['record_file']
        json_dir = staging / f'json_{key}'
        json_dir.mkdir()
        staged_json = json_dir / f"{key}_{Path(json_file).name}" if json_file else None
        results = staging / f"{key}_{Path(results_file).name}"
        shutil.copy2(results_file, results)
        if json_file:
            # The JSON records its results file relative to itself, for a refit to
            # update: here, the staged copy beside it.
            fit = json.loads(Path(json_file).read_text())
            fit['metadata']['results_file'] = os.path.relpath(results, json_dir)
            staged_json.write_text(json.dumps(fit, indent=2))
        session[field][exp] = {'output_dir': str(staging), 'json_folder': str(json_dir),
                               'json_file': str(staged_json) if staged_json else None,
                               'file_prefix': key,
                               'results_file': str(results),
                               'n_fitted': record['summary'].get('n_fitted'),
                               'frequency_mhz': float(record['args']['field_freq'])}
        page['fitted_experiments'].append(key)
        page[f'{key}_file'] = refs[f'{key}_file'] = input_file

    for record in records:
        args, summary = record['args'], record['summary']
        fill(record, field_of[float(args['field_freq'])], args['exp'].lower(),
             summary['results_file'], summary.get('json_file'), args['input'])

    for record in t1rho_records:
        args, summary = record['args'], record['summary']
        field = field_of[float(args['field_freq'])]
        fill(record, field, 't1rho', summary['t1rho_results'], summary.get('t1rho_json'),
             args['input'])
        if session[field]['t1'] is None:
            fill(record, field, 't1', summary['t1_results'], summary.get('t1_json'),
                 args['t1'])
        if session[field]['t2'] is None:
            table = staging / f'{field}_T2_from_T1rho.csv'
            shutil.copy2(summary['r2_table'], table)
            results = staging / f'{field}_T2_from_T1rho_fit_results.txt'
            write_t2_fit_results(pd.read_csv(table, dtype={'residue': str}), results)
            session[field]['t2'] = {'output_dir': str(staging), 'json_folder': None,
                                    'file_prefix': f'{field}_T2_from_T1rho',
                                    'results_file': str(results),
                                    'n_fitted': summary.get('n_residues'),
                                    'frequency_mhz': float(args['field_freq']),
                                    'derived_from': 't1rho'}
            page['fitted_experiments'].append(f'{field}_t2')
        prefix = '' if field == 'field1' else 'field2_'
        page[f'{prefix}theta_angle'] = args['theta']
        page[f'{prefix}omega1'] = args['omega1']
        page[f'{field}_carrier'] = args['carrier']
        page['peak_list_file'] = refs['peak_list_file'] = args['peaks']
    page['file_prefixes'] = sorted({e['file_prefix'] for f in session.values() if f
                                    for e in f.values() if isinstance(e, dict)})
    return {'state': _analysis_state(name, T1T2_PAGE, 't1t2', 't1t2', page, all_records),
            'file_refs': {'t1t2': refs}}


def modelfree_analysis(record, name, out):
    """An integrated (model-free) page analysis from one `dynamixs modelfree` record.

    `out` - the run's own outputs (stage_own_outputs) - is used as the page's folder: the
    results viewer opens `<prefix>_spectral_density_detailed.csv` there, which the CLI
    writes under that name. An R2 table given for a field (--f*-r2-table) has no field
    on the page; it is kept with the run in the analysis metadata.
    """
    args, out = record['args'], str(out)
    variant = 'jwh' if (args.get('method') or '').endswith('jwh') else '087'
    page = {
        'session_results': {'analysis_complete': True, 'results_prefix': args['prefix'],
                            'output_dir': out, 'json_folder': None},
        'field1_t1_file': args.get('f1_t1'), 'field1_t2_file': args.get('f1_t2'),
        'field1_noe_sat_file': args.get('f1_noe_sat'),
        'field1_noe_unsat_file': args.get('f1_noe_unsat'),
        'field2_t1_file': args.get('f2_t1'), 'field2_t2_file': args.get('f2_t2'),
        'field2_noe_sat_file': args.get('f2_noe_sat'),
        'field2_noe_unsat_file': args.get('f2_noe_unsat'),
        'field1_freq': args.get('field1_freq'), 'field2_freq': args.get('field2_freq'),
        'enable_dual_field': bool(args.get('dual')),
        'analysis_method': ('dual_' if args.get('dual') else 'single_') + variant,
        'rnh_angstrom': args.get('rnh'), 'csa_ppm': args.get('csa'),
        'initial_amplitude': args.get('init_amp'),
        't1_initial_time': args.get('init_t1'), 't2_initial_time': args.get('init_t2'),
        'bootstrap_iterations': args.get('n_bootstrap'),
        'monte_carlo_iterations': args.get('n_monte_carlo'),
        'output_dir': out, 'prefix': args['prefix'], 'file_prefixes': [args['prefix']],
        'json_folder_name': 'json',
    }
    page = {k: v for k, v in page.items() if v is not None}
    refs = {key: page[key] for key in page if key.endswith('_file')}
    return {'state': _analysis_state(name, INTEGRATED_PAGE, 'model_free', 'integrated',
                                     page, [record]),
            'file_refs': {'integrated': refs}}


def methyl_t2_analysis(record, name, out):
    """A methyl T2 page analysis from one `dynamixs methyl-t2` record.

    `out` - the run's own outputs (stage_own_outputs) - is used as the page's folder: the
    CLI writes the fit JSON (`*_methylT2_fit_data.json`, which the page's viewers look
    for) and the results file there, so both are the page's `last_json_folder` /
    `last_results_file`.
    """
    args, summary, out = record['args'], record['summary'], Path(out)
    page = {'input_file': args['input'], 'output_dir': str(out),
            'last_json_folder': str(out),
            'last_results_file': str(out / Path(summary['results_file']).name),
            'file_prefixes': [args['prefix']],
            'field_freq': args.get('field_freq'), 'bootstrap_iter': args.get('bootstrap')}
    page = {k: v for k, v in page.items() if v is not None}
    return {'state': _analysis_state(name, METHYL_T2_PAGE, 'methyl_t2', 'methyl_t2',
                                     page, [record]),
            'file_refs': {'methyl_t2': {'input_file': args['input']}}}


def hetnoe_analysis(record, name, out):
    """A hetNOE page analysis from one `dynamixs hetnoe` record: its two intensity files,
    and the results file it wrote, from `out` - the run's own outputs."""
    args, summary, out = record['args'], record['summary'], Path(out)
    page = {'mode': 1, 'sat_file': args['sat'], 'unsat_file': args['unsat'],
            'output_dir': str(out), 'prefix': args['prefix'],
            'file_prefixes': [args['prefix']],
            'last_results_file': str(out / Path(summary['output']).name)}
    return {'state': _analysis_state(name, HETNOE_PAGE, 'hetnoe', 'hetnoe', page, [record]),
            'file_refs': {'hetnoe': {'sat_file': args['sat'], 'unsat_file': args['unsat']}}}


def density_analysis(record, name, out):
    """A spectral-density page analysis from one `dynamixs density` record.

    The page has no viewer of its own: it holds the inputs, the parameters and the
    detailed results file the run wrote, from `out` - the run's own outputs.
    """
    args, summary, out = record['args'], record['summary'], Path(out)
    # A dual-field run writes _basic and _detailed tables, a single-field one _results.
    tables = [str(out / Path(o).name) for o in summary.get('outputs', [])
              if o.endswith('.csv')]
    detailed = next((o for o in tables if '_detailed' in Path(o).name),
                    tables[0] if tables else None)
    page = {'session_results': {'analysis_complete': True, 'results_prefix': args['prefix'],
                                'output_dir': str(out), 'results_file': detailed},
            'input_file1': args['input'], 'input_file2': args.get('input2'),
            'output_dir': str(out),
            'mode': 'dual_field' if args.get('dual') else 'single_field',
            'use_087': bool(args.get('use_087')),
            'field1_freq': args.get('field1_freq'), 'field2_freq': args.get('field2_freq'),
            'rnh_angstrom': args.get('rnh'), 'csa_ppm': args.get('csa'),
            'monte_carlo': bool(args.get('monte_carlo')), 'n_samples': args.get('n_samples'),
            'prefix': args['prefix'], 'file_prefixes': [args['prefix']]}
    page = {k: v for k, v in page.items() if v is not None}
    refs = {k: page[k] for k in ('input_file1', 'input_file2') if k in page}
    return {'state': _analysis_state(name, SPECTRAL_PAGE, 'spectral_density', 'spectral',
                                     page, [record]),
            'file_refs': {'spectral': refs}}


def stage_own_outputs(record, records, stage):
    """Copy into `stage` the files of the record's --out folder that its run wrote.

    The project save copies a page's results from its folder by name pattern, so runs
    sharing an --out folder would each be bundled with all of them. Every output a
    dynamixs command writes starts with `<prefix>_`; a file belongs to the longest prefix
    among the folder's records that it starts with, so `a` does not claim `a_b`'s files.
    """
    out = Path(record['args']['out']).resolve()
    prefix = record['args']['prefix']
    longer = [r['args']['prefix'] for r in records
              if Path(r['args']['out']).resolve() == out
              and r['args']['prefix'].startswith(f'{prefix}_')]
    stage = Path(stage)
    stage.mkdir(parents=True)
    for path in sorted(out.iterdir()):
        if (path.is_file() and path.name.startswith(f'{prefix}_')
                and not any(path.name.startswith(f'{p}_') for p in longer)):
            shutil.copy2(path, stage / path.name)
    return stage


def dynamixs_analyses(run_dir, staging):
    """{name: {state, file_refs}} for every DynamiXs analysis the run's records hold.

    A model-free or methyl T2 run is named after its prefix; two runs of one prefix are
    told apart by where their records sit, as two series of one folder name are.

    :return: (analyses, skipped) - `skipped` names the records no GUI page can hold
    """
    records = find_dynamixs_records(run_dir)
    analyses, skipped = {}, []
    t1t2 = [r for r in records if r['command'] == 'dynamixs t1t2']
    t1rho = [r for r in records if r['command'] == 'dynamixs t1rho']
    if t1t2 or t1rho:
        name = f'{Path(run_dir).name}_T1T2'
        analyses[name] = t1t2_analysis(t1t2, name, Path(staging) / name, t1rho)
    per_run = {'dynamixs modelfree': ('modelfree', modelfree_analysis),
               'dynamixs methyl-t2': ('methylT2', methyl_t2_analysis),
               'dynamixs hetnoe': ('hetNOE', hetnoe_analysis),
               'dynamixs density': ('spectralDensity', density_analysis)}
    for record in records:
        if record['command'] in per_run:
            suffix, build = per_run[record['command']]
            name = f"{record['args']['prefix']}_{suffix}"
            if name in analyses:
                place = Path(record['record_file']).parent.relative_to(run_dir)
                name = '_'.join(place.parts) + f'_{suffix}'
            analyses[name] = build(record, name, stage_own_outputs(
                record, records, Path(staging) / 'runs' / name))
        elif record['command'] not in ('dynamixs t1t2', 'dynamixs t1rho'):
            skipped.append(f"{record['command']} ({record['args'].get('prefix')})")
    return analyses, skipped


def _reference_spectrum(batch):
    """The series' first spectrum - its reference point - as an absolute path."""
    folder = batch.metadata.get('data_folder')
    for result in batch.results.values():
        name = result.get('spectrum_file')
        if name:
            path = Path(name) if Path(name).is_absolute() else Path(folder or '.') / name
            return str(path.resolve())
    return None


def _series_at(path):
    """The series in `path`, its output folder taken to be the one holding it.

    `series` writes batch_results.json beside its CSVs, so that is where they are. The
    folder the metadata records is where they were first written: gone when the run was
    moved, and when it was copied, a folder that may since hold a different run.
    """
    batch = batch_results_from_dict(json.loads(path.read_text()))
    folder = path.parent.resolve()
    batch.metadata['output_folder'] = str(folder)
    batch.metadata['csv_path'] = str(folder / Path(
        batch.metadata.get('csv_path') or 'series_analysis_tidy.csv').name)
    return batch


def build_project(run_dir, bundle, spectrum=None, peaks=None):
    """Write `bundle` from what the run under `run_dir` produced.

    :param spectrum: the spectrum the GUI reopens; defaults to the first series'
        reference spectrum
    :param peaks: the peak list it shows; defaults to the one that series was fitted with
    :return: {'series', 'kd', 'dynamixs': [...names], 'dynamixs_skipped': [...records
        no GUI page holds], 'spectrum', 'peaks'}
    :raises ValueError: when the run holds no series, Kd fit or DynamiXs record, or when
        no `spectrum` is given and the series' reference spectrum does not exist
    """
    series_files = find_series(run_dir)
    kd_files = find_kd_fits(run_dir)
    if not (series_files or kd_files or find_dynamixs_records(run_dir)):
        raise ValueError(f"No series (batch_results.json), Kd fit (*{KD_FIT_SUFFIX}) or "
                         f"DynamiXs record (*{DYNAMIXS_RECORD_SUFFIX}) found under {run_dir}")
    with tempfile.TemporaryDirectory(prefix='lunanmr_build_') as staging:
        return _build(run_dir, bundle, spectrum, peaks, series_files, kd_files, staging)


def _build(run_dir, bundle, spectrum, peaks, series_files, kd_files, staging):

    saved_series = {name: _series_at(path) for name, path in series_files.items()}
    first = next(iter(saved_series.values()), None)
    if spectrum is None and first is not None:
        spectrum = _reference_spectrum(first)
        # The bundle reopens it by path, and rebuilds the slim series' fits from it.
        if spectrum and not Path(spectrum).exists():
            raise ValueError(f"The series' reference spectrum {spectrum} does not exist - "
                             f"its spectra have moved since the series ran. Pass "
                             f"--spectrum <file> to name the spectrum the bundle reopens")
    if peaks is None and first is not None:
        peaks = first.metadata.get('peak_list_file')

    kd_analyses = {}
    for name, path in kd_files.items():
        fit_data = json.loads(path.read_text())
        source = saved_series.get(name)
        input_file = source.metadata.get('csv_path') if source is not None else None
        kd_analyses[name] = {'fit_data': fit_data,
                             'meta': kd_page_state(fit_data, path, input_file)}

    dynamixs, skipped = dynamixs_analyses(run_dir, staging)

    from lunaNMR.utils.file_manager import NMRFileManager
    peak_list = (NMRFileManager().load_peak_list(str(peaks)) if peaks
                 else pd.DataFrame())
    window = types.SimpleNamespace(
        current_nmr_file=str(Path(spectrum).resolve()) if spectrum else None,
        current_peak_file=str(Path(peaks).resolve()) if peaks else None,
        integrator=types.SimpleNamespace(peak_list=peak_list),
        saved_series=saved_series,
        kd_analyses=kd_analyses,
        dynamixs_analyses=dynamixs,
    )
    saved, _ = ProjectManager(window).save_project(Path(bundle))
    if not saved:
        raise RuntimeError(f"ProjectManager could not save {bundle}")
    return {'series': list(saved_series), 'kd': list(kd_analyses),
            'dynamixs': list(dynamixs), 'dynamixs_skipped': skipped,
            'spectrum': window.current_nmr_file, 'peaks': window.current_peak_file}

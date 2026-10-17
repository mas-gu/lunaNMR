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

import json
import shutil
import tempfile
import types
from pathlib import Path

import pandas as pd

from lunaNMR.utils.project_manager import ProjectManager, batch_results_from_dict

KD_FIT_SUFFIX = '_kd_fit_data.json'


def find_series(run_dir):
    """{name: batch_results.json path}, named after the folder each series wrote to."""
    found = {}
    for path in sorted(Path(run_dir).rglob('batch_results.json')):
        name = path.parent.name
        if name in found:   # two series folders of one name: tell them apart by place
            name = '_'.join(path.parent.relative_to(run_dir).parts)
        found[name] = path
    return found


def find_kd_fits(run_dir):
    """{name: fit JSON path}, named after the `--prefix` the fit was written with."""
    found = {}
    for path in sorted(Path(run_dir).rglob(f'*{KD_FIT_SUFFIX}')):
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
T1T2_PAGE, INTEGRATED_PAGE, METHYL_T2_PAGE = 1, 4, 5


def find_dynamixs_records(run_dir):
    """Every `<prefix>_dynamixs.json` a dynamixs command left under the run, each with the
    file it was read from as `record_file`."""
    records = []
    for path in sorted(Path(run_dir).rglob(f'*{DYNAMIXS_RECORD_SUFFIX}')):
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


def t1t2_analysis(records, name, staging):
    """ONE T1/T2 page analysis from a run's `dynamixs t1t2` records.

    The GUI page holds both fields and every experiment together; the CLI fits each one
    separately, into its own folder and often on the same default --prefix. Field 1 is
    the lower frequency. Each experiment's fit JSON goes into `json_<field>_<exp>/` under
    `staging` and its results file beside them under a name of the same key, so the
    project save copies one folder per page as it does for the GUI, and no experiment
    can land on another's file.

    :raises ValueError: for more than two fields - the page holds two - or one
        experiment fitted twice, of which the page can show only one
    """
    freqs = sorted({float(r['args']['field_freq']) for r in records})
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
    for record in records:
        args, summary = record['args'], record['summary']
        freq, exp = float(args['field_freq']), args['exp'].lower()
        field = field_of[freq]
        key = f'{field}_{exp}'
        if key in source:
            raise ValueError(f"{args['exp']} at {freq:g} MHz was fitted twice - "
                             f"{source[key]} and {record['record_file']}; the T1/T2 page "
                             f"holds one. Build from a folder holding the one to keep.")
        source[key] = record['record_file']
        json_dir = staging / f'json_{key}'
        json_dir.mkdir()
        if summary.get('json_file'):
            shutil.copy2(summary['json_file'], json_dir)
        results = staging / f"{key}_{Path(summary['results_file']).name}"
        shutil.copy2(summary['results_file'], results)
        session[field][exp] = {'output_dir': str(staging), 'json_folder': str(json_dir),
                               'results_file': str(results),
                               'n_fitted': summary.get('n_fitted'), 'frequency_mhz': freq}
        page['fitted_experiments'].append(key)
        page[f'{key}_file'] = refs[f'{key}_file'] = args['input']
    return {'state': _analysis_state(name, T1T2_PAGE, 't1t2', 't1t2', page, records),
            'file_refs': {'t1t2': refs}}


def modelfree_analysis(record, name):
    """An integrated (model-free) page analysis from one `dynamixs modelfree` record.

    Its output folder is used as the page's own: the results viewer opens
    `<prefix>_spectral_density_detailed.csv` there, which the CLI writes under that name.
    An R2 table given for a field (--f*-r2-table) has no field on the page; it is kept
    with the run in the analysis metadata.
    """
    args = record['args']
    variant = 'jwh' if (args.get('method') or '').endswith('jwh') else '087'
    page = {
        'session_results': {'analysis_complete': True, 'output_dir': args['out'],
                            'json_folder': None},
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
        'output_dir': args['out'], 'prefix': args['prefix'], 'json_folder_name': 'json',
    }
    page = {k: v for k, v in page.items() if v is not None}
    refs = {key: page[key] for key in page if key.endswith('_file')}
    return {'state': _analysis_state(name, INTEGRATED_PAGE, 'model_free', 'integrated',
                                     page, [record]),
            'file_refs': {'integrated': refs}}


def methyl_t2_analysis(record, name):
    """A methyl T2 page analysis from one `dynamixs methyl-t2` record.

    Its output folder is used as the page's own: the CLI writes the fit JSON
    (`*_methylT2_fit_data.json`, which the page's viewers look for) and the results file
    there, so both are the page's `last_json_folder` / `last_results_file`.
    """
    args, summary = record['args'], record['summary']
    page = {'input_file': args['input'], 'output_dir': args['out'],
            'last_json_folder': str(Path(summary['json_file']).parent),
            'last_results_file': summary['results_file'],
            'field_freq': args.get('field_freq'), 'bootstrap_iter': args.get('bootstrap')}
    page = {k: v for k, v in page.items() if v is not None}
    return {'state': _analysis_state(name, METHYL_T2_PAGE, 'methyl_t2', 'methyl_t2',
                                     page, [record]),
            'file_refs': {'methyl_t2': {'input_file': args['input']}}}


def dynamixs_analyses(run_dir, staging):
    """{name: {state, file_refs}} for every DynamiXs analysis the run's records hold.

    A model-free or methyl T2 run is named after its prefix; two runs of one prefix are
    told apart by where their records sit, as two series of one folder name are.

    :return: (analyses, skipped) - `skipped` names the records no GUI page can hold
    """
    records = find_dynamixs_records(run_dir)
    analyses, skipped = {}, []
    t1t2 = [r for r in records if r['command'] == 'dynamixs t1t2']
    if t1t2:
        name = f'{Path(run_dir).name}_T1T2'
        analyses[name] = t1t2_analysis(t1t2, name, Path(staging) / name)
    per_run = {'dynamixs modelfree': ('modelfree', modelfree_analysis),
               'dynamixs methyl-t2': ('methylT2', methyl_t2_analysis)}
    for record in records:
        if record['command'] in per_run:
            suffix, build = per_run[record['command']]
            name = f"{record['args']['prefix']}_{suffix}"
            if name in analyses:
                place = Path(record['record_file']).parent.relative_to(run_dir)
                name = '_'.join(place.parts) + f'_{suffix}'
            analyses[name] = build(record, name)
        elif record['command'] != 'dynamixs t1t2':
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


def build_project(run_dir, bundle, spectrum=None, peaks=None):
    """Write `bundle` from what the run under `run_dir` produced.

    :param spectrum: the spectrum the GUI reopens; defaults to the first series'
        reference spectrum
    :param peaks: the peak list it shows; defaults to the one that series was fitted with
    :return: {'series', 'kd', 'dynamixs': [...names], 'dynamixs_skipped': [...records
        no GUI page holds], 'spectrum', 'peaks'}
    :raises ValueError: when the run holds no series, Kd fit or DynamiXs record
    """
    series_files = find_series(run_dir)
    kd_files = find_kd_fits(run_dir)
    if not (series_files or kd_files or find_dynamixs_records(run_dir)):
        raise ValueError(f"No series (batch_results.json), Kd fit (*{KD_FIT_SUFFIX}) or "
                         f"DynamiXs record (*{DYNAMIXS_RECORD_SUFFIX}) found under {run_dir}")
    with tempfile.TemporaryDirectory(prefix='lunanmr_build_') as staging:
        return _build(run_dir, bundle, spectrum, peaks, series_files, kd_files, staging)


def _build(run_dir, bundle, spectrum, peaks, series_files, kd_files, staging):

    saved_series = {name: batch_results_from_dict(json.loads(path.read_text()))
                    for name, path in series_files.items()}
    first = next(iter(saved_series.values()), None)
    if spectrum is None and first is not None:
        spectrum = _reference_spectrum(first)
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

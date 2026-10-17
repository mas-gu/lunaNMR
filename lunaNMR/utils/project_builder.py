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

The spectrum and peak list the GUI reopens are the first series' reference spectrum and
the peak list it was fitted against, unless given. Both are referenced by absolute path,
never copied: a bundle moved to another machine needs the data at the same place.
"""

import json
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
    :return: {'series': [...names], 'kd': [...names], 'spectrum', 'peaks'}
    :raises ValueError: when the run holds no series and no Kd fit
    """
    series_files = find_series(run_dir)
    kd_files = find_kd_fits(run_dir)
    if not series_files and not kd_files:
        raise ValueError(f"No series (batch_results.json) or Kd fit "
                         f"(*{KD_FIT_SUFFIX}) found under {run_dir}")

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

    from lunaNMR.utils.file_manager import NMRFileManager
    peak_list = (NMRFileManager().load_peak_list(str(peaks)) if peaks
                 else pd.DataFrame())
    window = types.SimpleNamespace(
        current_nmr_file=str(Path(spectrum).resolve()) if spectrum else None,
        current_peak_file=str(Path(peaks).resolve()) if peaks else None,
        integrator=types.SimpleNamespace(peak_list=peak_list),
        saved_series=saved_series,
        kd_analyses=kd_analyses,
    )
    saved, _ = ProjectManager(window).save_project(Path(bundle))
    if not saved:
        raise RuntimeError(f"ProjectManager could not save {bundle}")
    return {'series': list(saved_series), 'kd': list(kd_analyses),
            'spectrum': window.current_nmr_file, 'peaks': window.current_peak_file}

# ABOUTME: Wide CSV + structured JSON writers for the Kd vs-sequence (ref->point) data.
# ABOUTME: Excel-ready table (blank=absent) + JSON mirror (null=absent), built from a fit JSON.

import csv
import json
import math

from kd_models import ref_vs_point_table, global_fit_table

_UNIT = {'csp': 'ppm', 'intensity': 'ratio'}


def _blank_if_nonfinite(v):
    return '' if isinstance(v, float) and not math.isfinite(v) else v


def _none_if_nonfinite(v):
    return None if isinstance(v, float) and not math.isfinite(v) else v


def export_ref_vs_point(base_path, fits, labels, obs, ref=0, alpha=0.14, value='height'):
    """Write '<base_path>.csv' (wide, Excel-ready) and '<base_path>.json' (structured
    mirror) for one observable's per-residue values across titration points, relative
    to the reference point. Returns the list of written file paths.

    CSV: rows = residue (sequence-ordered), columns = titration points, blank where a
    peak is absent. JSON: {observable, unit, reference_point, labels, residues:[{residue,
    values}]} with null where absent.
    """
    header, rows = ref_vs_point_table(fits, labels, obs, ref=ref, alpha=alpha, value=value)

    csv_path = base_path + ".csv"
    with open(csv_path, 'w', newline='') as fh:
        writer = csv.writer(fh)
        writer.writerow(header)
        writer.writerows(rows)

    json_path = base_path + ".json"
    residues = [{'residue': r[0],
                 'values': [None if c == '' else c for c in r[1:]]}
                for r in rows]
    payload = {
        'observable': obs,
        'unit': _UNIT.get(obs, ''),
        'reference_point': ref,
        'labels': list(labels),
        'residues': residues,
    }
    with open(json_path, 'w') as fh:
        json.dump(payload, fh, indent=2)

    return [csv_path, json_path]


def export_global_fit(base_path, fits, global_fit, obs, P0):
    """Write '<base_path>.csv' and '.json' with the per-residue params behind the global
    shared-Kd fit figure: each residue's amplitude(s), the single shared Kd, and the R^2
    of that one-Kd model against the residue's observed points. Returns the written file
    paths, or [] if there is no successful global fit for `obs` (no file is written).

    obs='intensity' columns: Residue, I0, I_inf, global_Kd, R2.
    obs='csp'       columns: Residue, dd_max, global_Kd, R2.
    """
    header, rows = global_fit_table(fits, global_fit, obs, P0)
    if not rows:
        return []

    csv_path = base_path + ".csv"
    with open(csv_path, 'w', newline='') as fh:
        writer = csv.writer(fh)
        writer.writerow(header)
        for r in rows:
            writer.writerow([r[0]] + [_blank_if_nonfinite(c) for c in r[1:]])

    keys = header[1:]
    residues = [{'residue': r[0],
                 **{k: _none_if_nonfinite(c) for k, c in zip(keys, r[1:])}}
                for r in rows]
    g = global_fit.get(obs) or {}
    json_path = base_path + ".json"
    payload = {
        'observable': obs,
        'unit': _UNIT.get(obs, ''),
        'global_Kd': _none_if_nonfinite(g.get('Kd')),
        'global_Kd_err': _none_if_nonfinite(g.get('Kd_err')),
        'residues': residues,
    }
    with open(json_path, 'w') as fh:
        json.dump(payload, fh, indent=2)

    return [csv_path, json_path]

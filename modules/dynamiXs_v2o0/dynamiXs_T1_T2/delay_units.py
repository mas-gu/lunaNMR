# ABOUTME: The delay units of a relaxation series: read from the series_metadata.json beside it,
# ABOUTME: or stated by the caller and refused when they contradict that record or the headers.

import csv
import json
import os

try:
    from delay_parser import parse_delay_column, stated_unit
except ImportError:
    from dynamiXs_T1_T2.delay_parser import parse_delay_column, stated_unit

# The delay units every relaxation fit falls back to: ms, because that is what the data is
# in. The delay parser returns every header in ms, `series` writes its matrices in ms, and
# the fitters themselves default to ms.
TIME_UNITS_DEFAULT = 'ms'

SIDECAR_NAME = 'series_metadata.json'


def series_sidecar(input_path):
    """The `series` run-metadata written beside a matrix, or None.

    `series` records the delay units it normalised to and which spectra had no
    parseable delay. Both were previously re-asserted by hand while the answer sat next
    to the input file.
    """
    path = os.path.join(os.path.dirname(os.path.abspath(input_path)), SIDECAR_NAME)
    if not os.path.isfile(path):
        return None
    try:
        with open(path) as fh:
            meta = json.load(fh)
    except (OSError, ValueError):
        return None            # a damaged sidecar must not block a runnable input
    if not isinstance(meta, dict):
        return None
    meta['path'] = path
    return meta


def headers_stating_units(input_path):
    """The delay headers of a matrix that write out their unit ('..._2o4s', '..._300ms').

    The delay parser has converted those to ms before anything is fitted, whatever unit
    they state.
    """
    if not os.path.isfile(input_path):
        return []          # the fitter reports the missing file in its own words
    with open(input_path, newline='') as fh:
        header = next(csv.reader(fh), [])
    return [h for h in header if stated_unit(h)]


def resolve_time_units(stated, sidecar, stated_by='--time-units', input_path=None):
    """Delay units: the ones `stated`, else the sidecar's, else the documented default.

    Units stated against the sidecar are refused rather than preferred. One of the two
    is wrong, and picking either way is how a T1 in seconds meets a T2 in milliseconds
    and puts R1 out by 1000x.

    So is any unit but ms on a matrix whose delay headers state their unit: the parser
    has already converted those to ms, and a seconds unit converts them a second time.
    `stated_by` names where the units were stated, for the message.
    """
    recorded = (sidecar or {}).get('value_units')
    if stated is None:
        return recorded or TIME_UNITS_DEFAULT
    if recorded and recorded != stated:
        raise ValueError(
            f"{stated_by} says {stated!r} but {sidecar['path']} records the delays as "
            f"{recorded!r}. One of them is wrong; guessing puts R1 out by 1000x. "
            f"Use the recorded units, or re-run `series` if the sidecar is stale."
        )
    stated_headers = (headers_stating_units(input_path)
                      if input_path and stated != 'ms' else [])
    if stated_headers:
        raise ValueError(
            f"{stated_by} says {stated!r} but the delay headers of {input_path} state "
            f"their own units, and are read as ms: {stated_headers[0]!r} is "
            f"{parse_delay_column(stated_headers[0]):g} ms. Reading those ms as "
            f"{stated!r} is a 1000x error. Use ms."
        )
    return stated

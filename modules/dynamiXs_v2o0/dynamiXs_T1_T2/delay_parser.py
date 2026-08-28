# ABOUTME: Shared delay-column header parser for the dynamiXs relaxation fitters.
# ABOUTME: Turns a spectrum-name / matrix column header into a delay value in ms.

import re


def parse_delay_column(col_name):
    """Parse a delay-column header to a delay in ms.

    Handles:
        bare numbers                  "300", "50.5"           → 300, 50.5
        duplicate-measurement suffix  "300_2", "50.5_2"       → 300, 50.5
        descriptive spectrum names    "003_T2_ADDA_3ms"       → 3
                                      "003_T2_ADDA_300ms"     → 300
        filesystem-safe decimal       "600_T1_0o3"            → 0.3
        unit conversion               "..._5s" / "..._500us"  → 5000 / 0.5
        unit + repeat marker          "..._300msb"            → 300
    Returns None when no delay token can be extracted.
    """
    col_str = str(col_name).strip()
    if not col_str:
        return None

    bare = re.fullmatch(r"(\d+(?:\.\d+)?)(?:_\d+)?", col_str)
    if bare:
        return float(bare.group(1))

    # Descriptive name: trailing numeric token, an optional unit (ms/s/us), and an
    # optional repeat-acquisition marker written as a single trailing letter (e.g.
    # the 'b' in '300msb' or '0o0b') or an '_2'-style suffix. 'o' stands in for '.'
    # to keep column headers filesystem-safe.
    tail = re.search(
        r"(?:^|_)(\d+(?:[o.]\d+)?)(ms|s|us)?([a-z]|_\d+)?$",
        col_str,
        flags=re.IGNORECASE,
    )
    if tail:
        num_part = tail.group(1).replace("o", ".")
        try:
            value = float(num_part)
        except ValueError:
            return None
        unit = (tail.group(2) or "").lower()
        # A trailing integer with no unit, no decimal separator and no letter marker
        # is an acquisition index, not a time: '03_2D_sample_ref_001' is a
        # spectrum name, and reading it as 1 ms builds a relaxation table out of
        # delays nobody measured. Any one of the three marks a real delay --
        # '..._102ms', '..._0o3', '..._51b' -- and a matrix written by `series` is
        # normalised to bare numbers, which never reach this branch.
        marker = tail.group(3) or ""
        has_separator = "o" in tail.group(1).lower() or "." in tail.group(1)
        has_letter_marker = bool(marker) and not marker.startswith("_")
        if not (unit or has_separator or has_letter_marker):
            return None
        if unit == "s":
            value *= 1000.0
        elif unit == "us":
            value /= 1000.0
        return value

    try:
        return float(col_str)
    except ValueError:
        return None


def parse_delay_columns(headers):
    """Parse delay-column headers, excluding those that carry no delay.

    Returns (values, keep, dropped): the parsed delays in header order, a boolean mask
    over `headers` selecting them, and the headers that carried none.

    Callers must slice their intensity columns with `keep`. Leaving a None in the array
    instead makes it dtype=object, and the first reduction numpy performs over it raises
    `TypeError: '<=' not supported between instances of 'float' and 'NoneType'` from
    inside numpy — a message that names neither the column nor the file. This is
    reachable on ordinary output: `series` mixes label kinds in one matrix, giving a
    stem-named column to any spectrum whose filename carried no parseable delay.
    """
    import numpy as np
    parsed = [parse_delay_column(h) for h in headers]
    keep = np.array([v is not None for v in parsed], dtype=bool)
    values = np.array([v for v in parsed if v is not None], dtype=float)
    dropped = [str(h) for h, v in zip(headers, parsed) if v is None]
    return values, keep, dropped


def require_delay_start(header_row, reserved, input_file):
    """Index of the first delay column, or raise saying why there is none.

    The failure mode this replaces: with no delay column found, a caller that defaults
    the start index to 0 reads the assignment column as intensities and dies with
    "could not convert string to float: '<a residue label>'" — a residue label, naming neither
    the file nor the cause. A whole experiment can legitimately have no delays (two
    hetNOE planes, or spectra named by acquisition index rather than delay), so this is
    a routine input error and deserves a routine message.
    """
    for i, col in enumerate(header_row):
        if col in reserved:
            continue
        if parse_delay_column(col) is not None:
            return i
    raise ValueError(
        f"No delay column found in {input_file}: no header carries a delay, so there "
        f"is nothing to fit against. Headers seen: {list(header_row)[:12]}. "
        f"Delay columns come from the "
        f"spectrum filenames: they need an explicit unit (_300ms, _2.4s, _2400us), a "
        f"decimal separator (_0o3), or a repeat marker (_300msb). A trailing plain "
        f"integer (_001) is read as an acquisition index, not a time."
    )

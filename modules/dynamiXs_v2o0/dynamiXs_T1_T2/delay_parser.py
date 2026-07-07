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
    Returns None when no delay token can be extracted.
    """
    col_str = str(col_name).strip()
    if not col_str:
        return None

    bare = re.fullmatch(r"(\d+(?:\.\d+)?)(?:_\d+)?", col_str)
    if bare:
        return float(bare.group(1))

    # Descriptive name: take the trailing numeric token (with optional unit and
    # duplicate-measurement suffix, either '_2'-style or a trailing letter like
    # '0o0b'). 'o' stands in for '.' to keep column headers filesystem-safe.
    tail = re.search(
        r"(?:^|_)(\d+(?:[o.]\d+)?)(ms|s|us)?[a-z]?(?:_\d+)?$",
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
        if unit == "s":
            value *= 1000.0
        elif unit == "us":
            value /= 1000.0
        return value

    try:
        return float(col_str)
    except ValueError:
        return None

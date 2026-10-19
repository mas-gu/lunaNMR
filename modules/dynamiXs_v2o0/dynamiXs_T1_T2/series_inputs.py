# ABOUTME: Tells a series' tidy table from the intensity matrix the relaxation fitters read,
# ABOUTME: and words the refusal so it names the matrix written beside the table.

import csv
import os

# Columns only the tidy table carries: one row per spectrum and peak, where a matrix has
# one row per residue and one column per delay.
TIDY_TABLE_COLUMNS = ('spectrum_name', 'assignment')
INTENSITY_MATRIX_NAME = 'peak_intensity_matrix.csv'


def tidy_table_refusal(path, matrix_name=INTENSITY_MATRIX_NAME):
    """Why `path` cannot be fitted as an intensity matrix, or None when it is not a tidy
    table (or cannot be read: the caller reports a missing file its own way)."""
    try:
        with open(path, newline='') as handle:
            header = next(csv.reader(handle), [])
    except (OSError, UnicodeDecodeError):
        return None
    if not all(column in header for column in TIDY_TABLE_COLUMNS):
        return None
    beside = os.path.join(os.path.dirname(os.path.abspath(path)), matrix_name)
    instead = beside if os.path.isfile(beside) else f"the series' {matrix_name}"
    return (f"{os.path.basename(path)} is a series' tidy table (one row per spectrum and "
            f"peak), not an intensity matrix: give {instead} instead")

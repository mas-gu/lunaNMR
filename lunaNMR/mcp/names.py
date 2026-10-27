# ABOUTME: Validates text that becomes part of an output filename, so it cannot name a path or a flag.
# ABOUTME: Prefixes and field names are joined into filenames by the CLI, which would follow a "../".

import re

NAME = re.compile(r'[A-Za-z0-9][A-Za-z0-9_.\-]{0,63}')


class InvalidName(ValueError):
    """A name that could leave the output folder or be read as something else."""


def safe_name(value, what):
    if not NAME.fullmatch(value):
        raise InvalidName(
            f'{what} must be 1-64 letters, digits, "_", "." or "-", starting with a letter or '
            f'digit: {value[:40]!r}')
    return value

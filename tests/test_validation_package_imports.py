# ABOUTME: The validation package exports only names that exist: a star import of it works.
# ABOUTME: Its tools are the modules verify_installation and spectra_check.

import importlib


def test_every_exported_name_exists():
    package = importlib.import_module('lunaNMR.validation')
    for name in getattr(package, '__all__', []):
        assert hasattr(package, name), name


def test_a_star_import_does_not_raise():
    exec('from lunaNMR.validation import *', {})

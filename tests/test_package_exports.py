# ABOUTME: Every name lunaNMR, lunaNMR.processors and lunaNMR.utils export resolves, so a star import works.
# ABOUTME: ParallelVoigtProcessor is reached at lunaNMR.core.parallel_voigt_processor, not the package.

import importlib

import pytest


def test_every_name_the_package_exports_resolves():
    lunanmr = importlib.import_module('lunaNMR')
    assert [name for name in lunanmr.__all__ if not hasattr(lunanmr, name)] == []


def test_a_star_import_of_the_package_succeeds():
    namespace = {}
    exec('from lunaNMR import *', namespace)
    assert 'EnhancedVoigtIntegrator' in namespace


def test_the_package_does_not_name_a_parallel_processor():
    """The name pointed at a module that never defined it; the live class is
    lunaNMR.core.parallel_voigt_processor.ParallelVoigtProcessor."""
    lunanmr = importlib.import_module('lunaNMR')
    assert 'ParallelVoigtProcessor' not in lunanmr.__all__
    assert 'ParallelVoigtProcessor' not in dir(lunanmr)
    with pytest.raises(AttributeError):
        lunanmr.ParallelVoigtProcessor


def test_every_name_the_processors_export_resolves():
    processors = importlib.import_module('lunaNMR.processors')
    assert [name for name in processors.__all__ if not hasattr(processors, name)] == []


def test_every_name_the_utils_export_resolves():
    utils = importlib.import_module('lunaNMR.utils')
    assert [name for name in utils.__all__ if not hasattr(utils, name)] == []

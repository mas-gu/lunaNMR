# ABOUTME: The DynamiXs background workers are the ones its pages start; the generic file-load
# ABOUTME: and batch workers, which no page started, are not offered.

import os
import sys
from pathlib import Path

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

_DIR = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_DIR))


def test_only_the_workers_a_page_starts_are_offered():
    import workers

    assert hasattr(workers, 'AnalysisWorker')
    assert hasattr(workers, 'IntegratedAnalysisParams')
    assert not hasattr(workers, 'FileLoadWorker')
    assert not hasattr(workers, 'BatchProcessingWorker')

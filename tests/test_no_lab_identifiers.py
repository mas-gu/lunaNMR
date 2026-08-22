# ABOUTME: Fails if a tracked file reintroduces a real construct name, mutant code
# ABOUTME: or calibrated spectrometer frequency as example or default data.
"""Demo defaults are the easiest place for real sample data to reach a public repo.

Example blocks in the dynamiXs scripts shipped a mutant-comparison study's own
identifiers -- variant codes, per-construct filenames, and the measured Larmor
frequencies of two specific magnets. Those say what is being studied before it
is published, so they are replaced with neutral placeholders and pinned here.

Round spectrometer classes (600, 700) are fine; it is the calibrated decimals
that identify an instrument.
"""

import re
import subprocess
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]

# Assembled from fragments so this file does not match its own patterns.
_MUTANTS = ["T" + "5D", "T" + "6D"]
_SAMPLES = ["A1_" + "WT", "data_in_" + "WT", "ZZ_" + "WT", "T1_" + "WT",
            "T2_" + "T6D", "WT_" + "density", "rex_analysis_" + "WT"]
_FREQUENCIES = ["700." + "093", "600." + "133"]
_ACQUISITIONS = ["NR_ATP_" + "ref_noCa", "594_" + "ce_D"]

FORBIDDEN = {
    "calibrated spectrometer frequency": _FREQUENCIES,
    "protein variant code": _MUTANTS,
    "per-construct example filename": _SAMPLES,
    "real acquisition or residue label": _ACQUISITIONS,
}

# Binary payloads and this test cannot meaningfully be scanned.
SKIP_SUFFIXES = {".joblib", ".pt", ".ft", ".ft1", ".ucsf", ".png", ".pdf"}
SKIP_PATHS = {
    Path(__file__).relative_to(REPO_ROOT).as_posix(),
    # Minified bundle: base64 payloads coincidentally contain these tokens.
    "docs/index.html",
}


def _tracked_text_files():
    out = subprocess.run(["git", "ls-files", "-z"], cwd=REPO_ROOT,
                         capture_output=True, text=True, check=True).stdout
    for rel in filter(None, out.split("\0")):
        if rel in SKIP_PATHS or Path(rel).suffix in SKIP_SUFFIXES:
            continue
        yield rel


@pytest.mark.parametrize("label,needles", sorted(FORBIDDEN.items()))
def test_no_tracked_file_contains(label, needles):
    pattern = re.compile("|".join(re.escape(n) for n in needles))
    offenders = []
    for rel in _tracked_text_files():
        try:
            text = (REPO_ROOT / rel).read_text(encoding="utf-8")
        except (UnicodeDecodeError, OSError):
            continue
        for lineno, line in enumerate(text.splitlines(), 1):
            match = pattern.search(line)
            if match:
                offenders.append(f"{rel}:{lineno} ({match.group(0)})")

    assert not offenders, (
        f"{len(offenders)} tracked line(s) contain a {label}:\n  "
        + "\n  ".join(offenders[:40])
        + ("\n  ..." if len(offenders) > 40 else "")
    )

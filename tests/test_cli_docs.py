# ABOUTME: Guards the agent-facing CLI documentation against claims the code does not honour.
# ABOUTME: Agents run these documents verbatim, so a stale flag or a duplicated number is a bug.

import re
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
DOCS = REPO_ROOT / "docs"
CLI_AGENT_DOC = DOCS / "CLI_AGENT.md"
RELAXATION_PLAYBOOK = DOCS / "CLI_AGENTS_DEEP" / "RELAXATION_PLAYBOOK.md"


def _agent_docs():
    """Every document an agent is pointed at, including the long-form runbooks."""
    return sorted(list(DOCS.glob("*.md")) + list(DOCS.glob("CLI_AGENTS_DEEP/*.md")))


def _table_rows(text):
    """Yield (first_cell, remaining_cells) for each markdown table row."""
    for line in text.splitlines():
        stripped = line.strip()
        if not stripped.startswith("|"):
            continue
        cells = [c.strip() for c in stripped.strip("|").split("|")]
        if cells:
            yield cells[0], cells[1:]


# A numeric band such as "0.80–0.90" or "1.10-1.25" (en dash or hyphen).
_NUMERIC_BAND = re.compile(r"\d+\.\d+\s*[–-]\s*\d+\.\d+")

# The QC checks whose value depends on the field pair, so no single band is correct.
_FIELD_RATIO_CHECKS = ("R1(high)/R1(low)", "R2(high)/R2(low)")


class TestTheFieldRatioBandsAreStatedOnce:
    """The R1/R2 field ratios follow from the dipolar + CSA rate expressions and change
    a lot with the field pair: 600/700 gives R1 0.81-0.84, 600/800 gives 0.68-0.73.
    A single universal band therefore fails a perfectly good 600/800 dataset. The
    numbers belong in RELAXATION_PLAYBOOK.md's per-field-pair table and nowhere else.
    """

    @pytest.mark.parametrize("doc", _agent_docs(), ids=lambda p: p.name)
    def test_no_qc_row_asserts_a_universal_band(self, doc):
        offenders = [
            (first, rest)
            for first, rest in _table_rows(doc.read_text())
            if first in _FIELD_RATIO_CHECKS and rest and _NUMERIC_BAND.search(rest[0])
        ]
        assert offenders == [], (
            f"{doc.relative_to(REPO_ROOT)} states a universal field-ratio band; "
            f"it depends on the field pair — defer to the per-field-pair table in "
            f"RELAXATION_PLAYBOOK.md: {offenders}"
        )

    def test_the_per_field_pair_table_is_still_there(self):
        """The deferral above is only safe while the table it defers to exists."""
        text = RELAXATION_PLAYBOOK.read_text()
        assert "| 600 / 700 |" in text
        assert "| 600 / 800 |" in text

    def test_the_agent_contract_points_at_the_table(self):
        """An agent reading only CLI_AGENT.md must still be able to reach the numbers."""
        assert "RELAXATION_PLAYBOOK.md" in CLI_AGENT_DOC.read_text()

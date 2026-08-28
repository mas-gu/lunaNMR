# ABOUTME: Guards the agent-facing CLI documentation against claims the code does not honour.
# ABOUTME: Agents run these documents verbatim, so a stale flag or a duplicated number is a bug.

import re
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
DOCS = REPO_ROOT / "docs"
CLI_AGENT_DOC = DOCS / "CLI_AGENT.md"
AGENTS_DOC = REPO_ROOT / "AGENTS.md"
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


class TestTheAgentDocsAreReachable:
    """The agent documentation is only worth writing if an agent can find it. The
    entry points are `--help`, AGENTS.md, and the two README indexes.
    """

    @pytest.fixture(scope="class")
    def top_level_help(self):
        import subprocess
        out = subprocess.run([sys.executable, "-m", "lunaNMR", "--help"],
                             cwd=str(REPO_ROOT), capture_output=True, text=True)
        assert out.returncode == 0, out.stderr
        return out.stdout

    @pytest.mark.parametrize("path", ["docs/CLI_AGENT.md", "docs/CLI_AGENTS_DEEP"])
    def test_help_names_the_agent_docs(self, top_level_help, path):
        """`python -m lunaNMR --help` is the one entry point every agent already runs."""
        assert path in top_level_help

    def test_agents_md_exists_and_names_both(self):
        text = AGENTS_DOC.read_text()
        assert "docs/CLI_AGENT.md" in text
        assert "docs/CLI_AGENTS_DEEP" in text

    def test_the_docs_index_lists_the_agent_reference(self):
        assert "CLI_AGENT.md" in (DOCS / "README.md").read_text()

    def test_the_project_readme_links_the_agent_reference(self):
        assert "CLI_AGENT.md" in (REPO_ROOT / "README.md").read_text()


# --------------------------------------------------------------- the generalised guard

_CMD_LINE = re.compile(r"python -m lunaNMR\s+(.+)")
# A flag token, excluding the shorthands the prose uses for families of flags:
# `--f{1,2}-t{1,2}-units`, `--f2-*` and `--f2-...` each stand for several real flags, and
# matching the fragment before the brace/star/ellipsis reports a flag nobody wrote. Flags
# are written in backticks throughout these documents, so a trailing '.' is a placeholder
# rather than the end of a sentence.
_FLAG = re.compile(r"(?<![\w-])(--[a-z][a-z0-9-]*)(?![-{*.\w])")


def _subcommand_tree(parser=None, prefix=()):
    """Every subcommand path -> its parser, walked from build_parser()."""
    import argparse
    from lunaNMR.cli import build_parser
    parser = parser or build_parser()
    out = {' '.join(prefix): parser} if prefix else {}
    for action in parser._actions:
        if isinstance(action, argparse._SubParsersAction):
            for name, sub in action.choices.items():
                out.update(_subcommand_tree(sub, prefix + (name,)))
    return out


def _accepted_flags(parser):
    return {opt for action in parser._actions
            for opt in action.option_strings if opt.startswith('--')}


def _batch_flags():
    """`batch` is dispatched before argparse sees it, so build_parser() knows none of
    its flags. They come from the parser it delegates to, which is the truth for it."""
    from lunaNMR.batch_processing.cli_interface import CLIInterface
    return _accepted_flags(CLIInterface().create_parser())


def _command_lines():
    """Every runnable `python -m lunaNMR ...` line in the agent-facing documents.

    Continuation backslashes are joined, so a wrapped multi-line invocation is checked
    as the single command an agent would paste.
    """
    lines = []
    for doc in _agent_docs():
        text = doc.read_text().replace('\\\n', ' ')
        for raw in text.splitlines():
            match = _CMD_LINE.search(raw)
            if match:
                lines.append((doc, match.group(1)))
    return lines


def _resolve(argv_text, tree):
    """The longest subcommand path this command line starts with, or None."""
    tokens = [t for t in argv_text.split() if not t.startswith('-')]
    for depth in (3, 2, 1):
        name = ' '.join(tokens[:depth])
        if name in tree:
            return name
    return None


class TestEveryDocumentedCommandIsReal:
    """These documents are read by agents that run what they say verbatim, so a flag
    that does not exist is a bug in the deliverable, not a typo.

    What this catches: a flag on a runnable command line that its subcommand does not
    accept; a flag named anywhere that no subcommand accepts; a subcommand nobody
    documented. What it does NOT catch: a flag attributed to the wrong subcommand in
    prose, when some other subcommand named in the same sentence does accept it -- that
    was the `density`/`modelfree` `--no-parallel` claim, and it needs a reader.
    """

    def test_there_are_command_lines_to_check(self):
        assert len(_command_lines()) > 10, "the extractor stopped matching"

    def test_every_flag_on_a_command_line_is_accepted_by_that_subcommand(self):
        tree = _subcommand_tree()
        universal = {'--help'}
        bad = []
        for doc, line in _command_lines():
            name = _resolve(line, tree)
            if name is None:
                continue                       # `--help`/`--version` only, nothing to bind
            accepted = _accepted_flags(tree[name]) | universal
            if name == 'batch':
                accepted |= _batch_flags()
            for flag in _FLAG.findall(line):
                if flag not in accepted:
                    bad.append(f"{doc.name}: `{name}` does not accept {flag}")
        assert bad == [], "\n".join(bad)

    def test_every_flag_mentioned_anywhere_exists_somewhere(self):
        """Catches a flag that was renamed or removed, wherever it is written."""
        tree = _subcommand_tree()
        known = set().union(*(_accepted_flags(p) for p in tree.values()))
        known |= _batch_flags() | {'--help', '--version'}
        # Installation instructions quote pip's flags, not this CLI's.
        pip_flags = {'--upgrade', '--user', '--editable'}
        bad = sorted({f"{doc.name}: {flag}" for doc in _agent_docs()
                      for flag in _FLAG.findall(doc.read_text())
                      if flag not in known and flag not in pip_flags})
        assert bad == [], "\n".join(bad)

    def test_every_subcommand_is_documented(self):
        """The direction the Kd-scoped version could not check: a subcommand nobody
        wrote down may as well not exist."""
        text = (DOCS / "CLI.md").read_text()
        # Word-boundary, not substring: `dynamixs t1rho` is a substring of a typo'd
        # `dynamixs t1rhoX`, so a plain `in` check certifies the misspelling.
        missing = sorted(name for name in _subcommand_tree()
                         if not re.search(r"\b" + re.escape(name) + r"\b", text))
        assert missing == [], missing

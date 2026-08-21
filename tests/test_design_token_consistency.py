# ABOUTME: Guards the GUI against design-system drift - unit mismatches and off-palette colours.
# ABOUTME: Font tokens are point values, so interpolating one into a px rule mis-sizes text off macOS.

import re
from pathlib import Path

import pytest

GUI_DIR = Path(__file__).parent.parent / "lunaNMR" / "gui"

# Colours the design system deliberately rejected: Apple's saturated system blue and
# red, in place of the softer PRIMARY_BUTTON_BG / DESTRUCTIVE_BUTTON_BG.
REJECTED_COLOURS = {
    "#007AFF": "PRIMARY_BUTTON_BG (#5B9EE5)",
    "#FF3B30": "DESTRUCTIVE_BUTTON_BG (#E8554E)",
}

FONT_TOKEN_AS_PX = re.compile(r"FONT_SIZE_[A-Z_]+\}px")


def gui_python_files():
    """Every Python source file under lunaNMR/gui, excluding caches."""
    return sorted(
        p for p in GUI_DIR.rglob("*.py")
        if "__pycache__" not in p.parts
    )


def test_gui_dir_is_where_we_think_it_is():
    """Guard against the scan silently covering nothing."""
    assert GUI_DIR.is_dir(), f"GUI directory not found at {GUI_DIR}"
    assert len(gui_python_files()) > 20


def test_font_size_tokens_are_never_interpolated_as_pixels():
    """FONT_SIZE_* tokens are point values; a px rule renders them ~33% small at 96 DPI.

    On macOS (logical DPI 72) pt and px coincide, so this mistake is invisible on the
    development platform and only shows on Windows/Linux.
    """
    offenders = []
    for path in gui_python_files():
        for lineno, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
            if FONT_TOKEN_AS_PX.search(line):
                offenders.append(f"{path.relative_to(GUI_DIR)}:{lineno}: {line.strip()}")

    assert not offenders, (
        f"{len(offenders)} font-size token(s) interpolated as px instead of pt:\n  "
        + "\n  ".join(offenders[:15])
    )


@pytest.mark.parametrize("colour,replacement", sorted(REJECTED_COLOURS.items()))
def test_rejected_system_colours_are_absent(colour, replacement):
    """The palette chose softer tones; the raw system colours must not creep back in."""
    pattern = re.compile(re.escape(colour), re.IGNORECASE)
    offenders = []
    for path in gui_python_files():
        if path.name.startswith("test_"):
            continue
        for lineno, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
            if pattern.search(line):
                offenders.append(f"{path.relative_to(GUI_DIR)}:{lineno}: {line.strip()}")

    assert not offenders, (
        f"{colour} used in {len(offenders)} place(s); use {replacement}:\n  "
        + "\n  ".join(offenders)
    )

# ABOUTME: Extracts series x-axis values from NMR filenames for time-series and titrations
# ABOUTME: Time mode reads _50ms/_1s; titration mode reads the _1o0/_0o5 suffix (o means .)

import os
import re
from typing import List, Optional, Tuple, Dict


class DelayExtractor:
    """
    Extract per-spectrum x-axis values from NMR filenames for series experiments.

    Two modes select what the trailing filename token means:

    Time mode (default) — relaxation delays:
    - T1_50ms.ft   -> 50 ms
    - T2_100ms.ft  -> 100 ms
    - T1_1s.ft     -> 1000 ms (converted from seconds)
    - T1_2.5s.ft   -> 2500 ms (fractional seconds)

    Titration mode — dimensionless titration points, with the filesystem-safe
    'o'-for-'.' convention:
    - sample_1o0.ft -> 1.0
    - sample_0o5.ft -> 0.5
    - sample_1.0.ft -> 1.0 (literal dot also accepted)
    - sample_2.ft   -> 2.0 (bare integer)

    The sorting / sequencing / column-naming machinery is value-agnostic and
    shared between modes; only value extraction differs.
    """

    # Patterns to match delay values (case-insensitive)
    # Support both with extension (_50ms.ft) and without (_50ms or T1_50ms)
    PATTERNS = [
        (r'_(\d+(?:\.\d+)?)ms(?:\.|$)', 1.0),      # _50ms. or _50ms (end of string)
        (r'_(\d+(?:\.\d+)?)s(?:\.|$)', 1000.0),    # _1s. or _1s (end of string)
    ]

    # Titration suffix: trailing _<value> where value uses 'o' or '.' for the
    # decimal point. An optional file extension may follow.
    TITRATION_PATTERN = re.compile(r'_(\d+(?:[o.]\d+)?)(?:\.\w+)?$')

    # NMR file extensions to look for when scanning folders
    NMR_EXTENSIONS = {'.ft', '.ucsf', '.pipe', '.2rr', '.2ii'}

    def __init__(self, mode: str = "time"):
        if mode not in ("time", "titration"):
            raise ValueError(f"mode must be 'time' or 'titration', got {mode!r}")
        self.mode = mode

    def extract_value(self, filename: str) -> Optional[float]:
        """
        Extract the series x-axis value for a filename according to the mode.

        Time mode returns a delay in milliseconds; titration mode returns the
        dimensionless titration point. Returns None if the filename has no
        value parseable in the current mode.
        """
        if self.mode == "titration":
            return self._extract_titration(filename)
        return self.extract_delay_ms(filename)

    def _extract_titration(self, filename: str) -> Optional[float]:
        """Extract a dimensionless titration point from the _<value> suffix."""
        basename = os.path.basename(filename)
        match = self.TITRATION_PATTERN.search(basename)
        if not match:
            return None
        token = match.group(1).replace('o', '.')
        try:
            return float(token)
        except ValueError:
            return None

    def extract_delay_ms(self, filename: str) -> Optional[float]:
        """
        Extract delay value in milliseconds from a filename.

        Args:
            filename: The filename (can include path) to parse

        Returns:
            Delay value in milliseconds, or None if no delay found
        """
        # Get just the filename if path is provided
        basename = os.path.basename(filename)

        for pattern, multiplier in self.PATTERNS:
            match = re.search(pattern, basename, re.IGNORECASE)
            if match:
                value = float(match.group(1))
                return value * multiplier

        return None

    def sort_files_by_delay(
        self, files: List[str]
    ) -> List[Tuple[str, float]]:
        """
        Sort a list of filenames by their extracted delay values.

        Files without valid delays are excluded from the result.

        Args:
            files: List of filenames to sort

        Returns:
            List of (filename, delay_ms) tuples, sorted by delay ascending
        """
        files_with_delays = []

        for filename in files:
            delay = self.extract_value(filename)
            if delay is not None:
                files_with_delays.append((filename, delay))

        # Sort by delay value
        files_with_delays.sort(key=lambda x: x[1])

        return files_with_delays

    def scan_folder(self, folder_path: str) -> List[Tuple[str, float]]:
        """
        Scan a folder for NMR files with delay values in their names.

        Args:
            folder_path: Path to the folder to scan

        Returns:
            List of (full_path, delay_ms) tuples, sorted by delay ascending.
            Returns empty list if folder doesn't exist or has no matching files.
        """
        if not os.path.isdir(folder_path):
            return []

        files_with_delays = []

        try:
            for filename in os.listdir(folder_path):
                # Check if it's likely an NMR file
                _, ext = os.path.splitext(filename.lower())
                if ext not in self.NMR_EXTENSIONS:
                    continue

                full_path = os.path.join(folder_path, filename)
                if not os.path.isfile(full_path):
                    continue

                delay = self.extract_value(filename)
                if delay is not None:
                    files_with_delays.append((full_path, delay))

        except PermissionError:
            return []

        # Sort by delay value
        files_with_delays.sort(key=lambda x: x[1])

        return files_with_delays

    def get_delay_summary(self, folder_path: str) -> dict:
        """
        Get a summary of delays found in a folder.

        Args:
            folder_path: Path to the folder to scan

        Returns:
            Dict with 'files': list of files, 'delays': list of delays,
            'min': minimum delay, 'max': maximum delay, 'count': number of files
        """
        result = self.scan_folder(folder_path)

        if not result:
            return {
                'files': [],
                'delays': [],
                'min': None,
                'max': None,
                'count': 0
            }

        files = [f for f, _ in result]
        delays = [d for _, d in result]

        return {
            'files': files,
            'delays': delays,
            'min': min(delays),
            'max': max(delays),
            'count': len(delays)
        }

    def sort_files_with_sequence(
        self, files: List[str]
    ) -> List[Tuple[str, float, int]]:
        """
        Sort files by delay and assign sequence numbers to duplicates.

        Files with the same delay value are grouped together. Within the same
        delay, files are sorted alphabetically so that replicates (marked with
        _b_, _c_, etc.) come after the original file.

        Args:
            files: List of filenames to sort

        Returns:
            List of (filename, delay_ms, sequence) tuples, sorted by delay
            then alphabetically within the same delay.
        """
        # Extract delays
        files_with_delays = []
        for filename in files:
            delay = self.extract_value(filename)
            if delay is not None:
                # Get basename for sorting (handles full paths)
                basename = os.path.basename(filename)
                files_with_delays.append((filename, delay, basename))

        # Sort by delay, then alphabetically by basename
        # This ensures _600ms comes before _b_600ms
        files_with_delays.sort(key=lambda x: (x[1], x[2]))

        # Assign sequence numbers within each delay group
        result = []
        delay_counts: Dict[float, int] = {}

        for filename, delay, _ in files_with_delays:
            delay_counts[delay] = delay_counts.get(delay, 0) + 1
            sequence = delay_counts[delay]
            result.append((filename, delay, sequence))

        return result

    def get_column_name(self, delay_ms: float, sequence: int) -> str:
        """
        Generate a unique column name for a delay value.

        First occurrence uses just the delay value (e.g., "50").
        Subsequent occurrences append the sequence number (e.g., "50_2", "50_3").

        Args:
            delay_ms: Extracted value (delay in ms, or titration point)
            sequence: Sequence number (1 for first occurrence, 2 for second, etc.)

        Returns:
            Unique column name string
        """
        # Format with minimal digits: whole numbers as ints ("50"), fractional
        # values at full precision ("0.15") so titration points are not truncated.
        if delay_ms == int(delay_ms):
            delay_str = str(int(delay_ms))
        else:
            delay_str = f"{delay_ms:g}"

        if sequence == 1:
            return delay_str
        else:
            return f"{delay_str}_{sequence}"

    def parse_column_name(self, name: str) -> Optional[Tuple[float, int]]:
        """
        Parse a column name to extract delay and sequence number.

        Inverse of get_column_name(). Handles formats like:
        - "50" -> (50.0, 1)
        - "50_2" -> (50.0, 2)
        - "100.5" -> (100.5, 1)
        - "100.5_3" -> (100.5, 3)

        Args:
            name: Column name string

        Returns:
            Tuple of (delay_ms, sequence) or None if not parseable
        """
        # Match: number (possibly with decimal) optionally followed by _N
        match = re.match(r'^(\d+(?:\.\d+)?)(?:_(\d+))?$', name)
        if match:
            try:
                delay = float(match.group(1))
                sequence = int(match.group(2)) if match.group(2) else 1
                return (delay, sequence)
            except ValueError:
                return None
        return None

    def scan_folder_with_sequence(
        self, folder_path: str
    ) -> List[Tuple[str, float, int]]:
        """
        Scan a folder for NMR files with delay values, including sequence numbers.

        Args:
            folder_path: Path to the folder to scan

        Returns:
            List of (full_path, delay_ms, sequence) tuples, sorted by delay
            with sequence numbers for duplicates.
        """
        if not os.path.isdir(folder_path):
            return []

        files_with_delays_and_idx = []

        try:
            for idx, filename in enumerate(sorted(os.listdir(folder_path))):
                # Check if it's likely an NMR file
                _, ext = os.path.splitext(filename.lower())
                if ext not in self.NMR_EXTENSIONS:
                    continue

                full_path = os.path.join(folder_path, filename)
                if not os.path.isfile(full_path):
                    continue

                delay = self.extract_value(filename)
                if delay is not None:
                    files_with_delays_and_idx.append((full_path, delay, idx))

        except PermissionError:
            return []

        # Sort by delay, then by original index
        files_with_delays_and_idx.sort(key=lambda x: (x[1], x[2]))

        # Assign sequence numbers
        result = []
        delay_counts: Dict[float, int] = {}

        for full_path, delay, _ in files_with_delays_and_idx:
            delay_counts[delay] = delay_counts.get(delay, 0) + 1
            sequence = delay_counts[delay]
            result.append((full_path, delay, sequence))

        return result

    def build_column_mapping(self, files: List[str]) -> Dict[str, str]:
        """
        Build a mapping from filenames to unique column names.

        This is the main entry point for processing series with duplicate delays.
        Returns a dict that maps each filename to its unique column name.

        Args:
            files: List of filenames to process

        Returns:
            Dict mapping filename -> unique column name (e.g., "50", "50_2")
        """
        sorted_files = self.sort_files_with_sequence(files)

        mapping = {}
        for filename, delay, sequence in sorted_files:
            mapping[filename] = self.get_column_name(delay, sequence)

        return mapping

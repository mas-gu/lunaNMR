# ABOUTME: Confines every path an MCP client names to one data folder, symlinks resolved.
# ABOUTME: Refuses an output folder that already exists, so a client cannot overwrite a run.

from pathlib import Path


class PathRefused(ValueError):
    """A path the client named is not one the server will touch."""


class DataRoot:
    def __init__(self, path):
        resolved = Path(path).expanduser().resolve()
        if not resolved.is_dir():
            raise PathRefused(f"data root is not an existing folder: {path}")
        self.path = resolved

    def _inside(self, name):
        candidate = Path(name).expanduser()
        if not candidate.is_absolute():
            candidate = self.path / candidate
        resolved = candidate.resolve()
        if resolved != self.path and self.path not in resolved.parents:
            raise PathRefused(f"outside the data root ({self.path}): {name}")
        return resolved

    def existing(self, name):
        """A file or folder that must already be there."""
        resolved = self._inside(name)
        if not resolved.exists():
            raise PathRefused(f"does not exist: {name}")
        return resolved

    def new_output(self, name):
        """A folder a run will create: it must not exist yet."""
        resolved = self._inside(name)
        if resolved.exists():
            raise PathRefused(f"already exists, choose a new output folder: {name}")
        return resolved

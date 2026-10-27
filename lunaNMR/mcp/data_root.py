# ABOUTME: Confines every path an MCP client names to one data folder, symlinks resolved.
# ABOUTME: Refuses an output folder that already exists, so a client cannot overwrite a run.

import os
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

    def new_file(self, name):
        """A file a run will create: inside the root, in a folder that exists, not yet there."""
        resolved = self._inside(name)
        if resolved.exists():
            raise PathRefused(f"already exists, choose a new file name: {name}")
        if not resolved.parent.is_dir():
            raise PathRefused(f"its folder does not exist: {name}")
        return resolved

    def writable_file(self, name):
        """A file path a run will write or merge into: inside the root, in a folder that exists."""
        resolved = self._inside(name)
        if resolved.is_dir():
            raise PathRefused(f"is a folder, not a file: {name}")
        if not resolved.parent.is_dir():
            raise PathRefused(f"its folder does not exist: {name}")
        return resolved

    def refuse_escaping_links(self, folder):
        """Scanning tools follow links on their own, so a link inside the folder that leads out
        of the root would read outside it. Refuse the folder, naming the first such link."""
        for directory, subfolders, files in os.walk(folder):
            for name in subfolders + files:
                link = Path(directory) / name
                if not link.is_symlink():
                    continue
                target = link.resolve()
                if target != self.path and self.path not in target.parents:
                    raise PathRefused(
                        f"{link.relative_to(self.path)} is a link to {target}, outside the "
                        f"data root ({self.path}); remove the link or copy the data in")

    def relative(self, path):
        return Path(path).resolve().relative_to(self.path)

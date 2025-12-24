# ABOUTME: Unified output routing system for print-to-signal migration.
# ABOUTME: Routes messages to callbacks (GUI) or print() (CLI) based on context.

from enum import Enum
from typing import Callable, Optional
import threading


class OutputLevel(Enum):
    """Output severity levels for message routing."""
    PROGRESS = 1   # Status updates: "Processing peak X..."
    INFO = 2       # Important user-facing information
    WARNING = 3    # Non-fatal issues
    ERROR = 4      # Failures


class OutputManager:
    """
    Thread-local output routing for GUI/CLI compatibility.

    When a callback is set for the current thread, messages are routed to it.
    When no callback is set, messages fall back to print().

    Usage:
        # In worker thread setup
        OutputManager.set_callback(my_progress_callback)

        # In processing code
        OutputManager.emit(OutputLevel.PROGRESS, "Processing...", progress=50)

        # Or use convenience functions
        from lunaNMR.utils.output_manager import log_progress
        log_progress("Processing...", progress=50, task="Fitting")
    """

    _local = threading.local()

    @classmethod
    def set_callback(cls, callback: Optional[Callable[[int, str, str, bool], None]]):
        """
        Set the progress callback for the current thread.

        Args:
            callback: Function with signature (progress, task, log_msg, failed)
                     - progress: Percentage (0-100) or -1 if indeterminate
                     - task: Short task description for status display
                     - log_msg: Detailed message for log area
                     - failed: True if this is a failure message
        """
        cls._local.callback = callback

    @classmethod
    def get_callback(cls) -> Optional[Callable]:
        """Get the callback for the current thread, or None."""
        return getattr(cls._local, 'callback', None)

    @classmethod
    def emit(
        cls,
        level: OutputLevel,
        message: str,
        progress: int = -1,
        task: str = "",
        failed: bool = False
    ):
        """
        Emit a message to callback or stdout.

        Args:
            level: Message severity (PROGRESS, INFO, WARNING, ERROR)
            message: The message content
            progress: Progress percentage (0-100) or -1 for indeterminate
            task: Short task description (defaults to message if empty)
            failed: Whether this indicates a failure
        """
        callback = cls.get_callback()

        if callback is not None:
            # Route to GUI callback
            callback(progress, task or message, message, failed)
        else:
            # Fallback to print for CLI usage
            print(message)


def log_progress(msg: str, progress: int = -1, task: str = ""):
    """Log a progress/status message."""
    OutputManager.emit(OutputLevel.PROGRESS, msg, progress=progress, task=task)


def log_info(msg: str):
    """Log an informational message."""
    OutputManager.emit(OutputLevel.INFO, msg)


def log_warning(msg: str):
    """Log a warning message."""
    OutputManager.emit(OutputLevel.WARNING, f"Warning: {msg}")


def log_error(msg: str, failed: bool = True):
    """Log an error message."""
    OutputManager.emit(OutputLevel.ERROR, msg, failed=failed)

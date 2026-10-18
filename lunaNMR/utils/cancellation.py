# ABOUTME: The stop event a caller passes to a fit or a series, and the exception a stopped run raises.
# ABOUTME: A run with no event (the CLI) never stops early; a stopped run returns nothing.


class Cancelled(Exception):
    """The run's stop event was set: it ended early and returns nothing."""


def raise_if_stopped(stop):
    """Raise Cancelled when the stop event (a threading.Event, or None) is set."""
    if stop is not None and stop.is_set():
        raise Cancelled()

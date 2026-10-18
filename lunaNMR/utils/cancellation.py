# ABOUTME: The stop event a caller passes to a fit or a series, and the exception a stopped run raises.
# ABOUTME: A run with no event (the CLI) never stops early; a stopped run returns nothing.

import time


class Cancelled(Exception):
    """The run's stop event was set: it ended early and returns nothing."""


def raise_if_stopped(stop):
    """Raise Cancelled when the stop event (a threading.Event, or None) is set."""
    if stop is not None and stop.is_set():
        raise Cancelled()


def result_unless_stopped(async_result, timeout, stop):
    """A pool task's result (timeout None: no limit). With a stop event the wait is
    polled, so that setting the event ends it (Cancelled) instead of the task's own end."""
    if stop is None:
        return async_result.get(timeout=timeout)
    deadline = time.monotonic() + (timeout if timeout is not None else float('inf'))
    while not async_result.ready() and time.monotonic() < deadline:
        raise_if_stopped(stop)
        async_result.wait(0.1)
    raise_if_stopped(stop)
    return async_result.get(timeout=0 if timeout is not None else None)

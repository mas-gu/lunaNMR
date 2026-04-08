# Fix: QThread Crash on Ubuntu After Multicore Analysis

## Issue

Application crashes on Ubuntu after T1/T2 multicore analysis completes:

```
Multicore analysis completed successfully!
QThread: Destroyed while thread '' is still running
Aborted (core dumped)
/home/masgu/anaconda3/lib/python3.13/multiprocessing/resource_tracker.py:301: UserWarning: resource_tracker: There appear to be 6 leaked semaphore objects to clean up at shutdown
```

## Root Cause

In `modules/dynamiXs_v2o0/dynamiXs_gui.py`, the `T1T2FittingWorker` (a QThread) is created and started but never properly cleaned up after completion. When the GUI closes or the worker reference changes, Qt destroys the QThread while it may still be running.

**Location**: `dynamiXs_gui.py` lines 1486-1491 (worker creation) and 1942-1983 (`_on_finished` without cleanup)

## Fix

### File: `lunaNMR_v1o0/modules/dynamiXs_v2o0/dynamiXs_gui.py`

### Edit 1: Add cleanup at end of `_on_finished` (after line 1983)

Change:
```python
        self.results_text.appendPlainText(f"\nSession updated: {session_key} stored")
        self.results_text.appendPlainText(f"Fitted experiments: {sorted(self.fitted_experiments)}")
```

To:
```python
        self.results_text.appendPlainText(f"\nSession updated: {session_key} stored")
        self.results_text.appendPlainText(f"Fitted experiments: {sorted(self.fitted_experiments)}")

        # Clean up worker thread to prevent crash on exit
        if hasattr(self, 'worker') and self.worker is not None:
            self.worker.wait()
            self.worker.deleteLater()
            self.worker = None
```

### Edit 2: Add cleanup in `_on_error` (after line 1989)

Change:
```python
    @Slot(str)
    def _on_error(self, error_msg: str):
        """Handle analysis error."""
        self.progress_bar.hide()
        self.results_text.appendPlainText(f"\nError during analysis:\n{error_msg}")
```

To:
```python
    @Slot(str)
    def _on_error(self, error_msg: str):
        """Handle analysis error."""
        self.progress_bar.hide()
        self.results_text.appendPlainText(f"\nError during analysis:\n{error_msg}")

        # Clean up worker thread to prevent crash on exit
        if hasattr(self, 'worker') and self.worker is not None:
            self.worker.wait()
            self.worker.deleteLater()
            self.worker = None
```

## Explanation

- `worker.wait()`: Blocks until the QThread's `run()` method has fully completed
- `worker.deleteLater()`: Schedules the QThread object for deletion by Qt's event loop (thread-safe)
- `self.worker = None`: Removes our reference so we don't accidentally reuse a deleted object

## Testing

1. Run a T1 or T2 multicore analysis
2. Wait for "Multicore analysis completed successfully!"
3. Close the application
4. Verify no crash or semaphore warnings appear

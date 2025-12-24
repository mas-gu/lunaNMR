"""ABOUTME: Qt progress dialog with thread-safe updates and detailed logging
ABOUTME: Implements the threading pattern for all background operations in lunaNMR
"""

from datetime import datetime, timedelta
from pathlib import Path
from typing import Optional
import logging

from PySide6.QtWidgets import (
    QVBoxLayout,
    QHBoxLayout,
    QLabel,
    QProgressBar,
    QTextEdit,
    QPushButton,
    QFrame,
    QWidget,
    QFileDialog,
    QMessageBox
)
from PySide6.QtCore import Signal, Slot, Qt
from PySide6.QtGui import QFont

from lunaNMR.gui.base.base_dialog import BaseDialog
from lunaNMR.gui.styles.design_system import (
    BG_COLOR,
    PANEL_BG_COLOR,
    FRAME_BG_COLOR,
    PRIMARY_TEXT,
    SECONDARY_TEXT,
    SUCCESS_GREEN,
    WARNING_ORANGE,
    ERROR_RED,
    PRIMARY_BUTTON_BG,
    PRIMARY_BUTTON_HOVER,
    SECONDARY_BUTTON_BG,
    SECONDARY_BUTTON_HOVER,
    DESTRUCTIVE_BUTTON_BG,
    DESTRUCTIVE_BUTTON_HOVER,
    SPACING_SM,
    SPACING_MD,
    SPACING_LG,
    FRAME_CORNER_RADIUS,
)

logger = logging.getLogger(__name__)


class ProgressDialog(BaseDialog):
    """
    Thread-safe progress dialog for long-running operations.

    This class establishes the threading pattern used throughout lunaNMR:
    - Worker threads emit signals
    - Signals are connected to slots in the main (GUI) thread
    - Slots update GUI elements safely
    - NO direct GUI updates from worker threads

    Features:
        - Thread-safe progress updates via Qt signals
        - Real-time status messages
        - Detailed logging with timestamps
        - ETA calculation and display
        - Pause/Resume/Cancel functionality
        - Save log to file
        - Non-modal operation (user can continue working)

    Signals:
        progress_updated(int, str): Progress value (0-100) and status message
        log_message(str, bool): Log message and whether it's an error
        operation_completed(str): Final completion message

    Example Usage:
        ```python
        from PySide6.QtCore import QThread, Signal
        from lunaNMR.gui.components.progress_dialog import ProgressDialog

        class WorkerThread(QThread):
            # Define signals for communication
            progress = Signal(int, str)  # progress, status
            log = Signal(str, bool)  # message, is_error
            finished_work = Signal(str)  # completion message

            def run(self):
                for i in range(100):
                    if self.isInterruptionRequested():
                        break

                    # Emit progress signal (thread-safe)
                    self.progress.emit(i, f"Processing item {i}")

                    # Emit log signal (thread-safe)
                    self.log.emit(f"Completed item {i}", False)

                    time.sleep(0.1)

                self.finished_work.emit("Processing complete!")

        # Create progress dialog
        dialog = ProgressDialog(parent=main_window, title="Processing Data")

        # Create worker thread
        worker = WorkerThread()

        # Connect worker signals to dialog slots
        worker.progress.connect(dialog.update_progress)
        worker.log.connect(dialog.add_log_message)
        worker.finished_work.connect(dialog.complete)

        # Connect dialog cancel to worker
        dialog.cancelled.connect(worker.requestInterruption)

        # Start work and show dialog
        worker.start()
        dialog.show()  # Non-modal, user can continue working
        ```

    Threading Best Practices:
        1. NEVER update GUI from worker thread directly
        2. ALWAYS use signals/slots for cross-thread communication
        3. Use QThread.requestInterruption() for cancellation
        4. Check isInterruptionRequested() in worker loop
        5. Emit progress signals frequently for responsive UI
        6. Keep worker logic separate from GUI logic
    """

    # Signals for external communication
    cancelled = Signal()  # Emitted when user clicks Cancel
    paused = Signal(bool)  # Emitted when user pauses/resumes (True = paused)

    def __init__(
        self,
        parent=None,
        title: str = "Processing",
        show_details: bool = True,
        modal: bool = False  # Non-modal by default
    ):
        """
        Initialize progress dialog.

        Args:
            parent: Parent widget
            title: Dialog title
            show_details: Whether to show detailed logging section
            modal: Whether dialog should be modal (blocks parent)
        """
        super().__init__(
            parent=parent,
            title=title,
            default_size=(600, 550),
            min_size=(600, 550),
            modal=modal
        )

        # Progress tracking
        self.start_time = datetime.now()
        self.total_tasks = 0
        self.completed_tasks = 0
        self.failed_tasks = 0

        # State flags
        self._is_cancelled = False
        self._is_completed = False
        self._is_paused = False

        # Store whether to show details
        self._show_details = show_details

        # Build UI
        self._create_widgets()

        # Center on parent or screen
        if parent:
            self.center_on_parent()
        else:
            self.center_on_screen()

    def _create_widgets(self):
        """Create and layout all dialog widgets."""
        # Main layout
        main_layout = QVBoxLayout()
        main_layout.setSpacing(SPACING_MD)
        main_layout.setContentsMargins(SPACING_LG, SPACING_LG, SPACING_LG, SPACING_LG)

        # Create sections
        self._create_progress_section(main_layout)
        self._create_task_section(main_layout)
        self._create_stats_section(main_layout)

        if self._show_details:
            self._create_log_section(main_layout)

        # Add stretch before buttons
        main_layout.addStretch()

        # Create control buttons at bottom
        self._create_button_section(main_layout)

        self.setLayout(main_layout)

    def _create_progress_section(self, parent_layout: QVBoxLayout):
        """Create progress bar section."""
        # Progress bar with styled container
        progress_frame = QFrame()
        progress_frame.setObjectName("progressFrame")
        progress_frame.setStyleSheet(f"""
            QFrame#progressFrame {{
                background-color: {PANEL_BG_COLOR};
                border-radius: {FRAME_CORNER_RADIUS}px;
            }}
        """)

        progress_layout = QVBoxLayout(progress_frame)
        progress_layout.setContentsMargins(4, 4, 4, 4)

        # Create progress bar
        self.progress_bar = QProgressBar()
        self.progress_bar.setRange(0, 100)
        self.progress_bar.setValue(0)
        self.progress_bar.setTextVisible(True)
        self.progress_bar.setFormat("%v%")
        self.progress_bar.setMinimumHeight(20)

        # Style the progress bar
        self.progress_bar.setStyleSheet(f"""
            QProgressBar {{
                border: none;
                border-radius: 4px;
                background-color: {PANEL_BG_COLOR};
                text-align: center;
                color: {PRIMARY_TEXT};
            }}
            QProgressBar::chunk {{
                background-color: {SUCCESS_GREEN};
                border-radius: 4px;
            }}
        """)

        progress_layout.addWidget(self.progress_bar)
        parent_layout.addWidget(progress_frame)

    def _create_task_section(self, parent_layout: QVBoxLayout):
        """Create current task display section."""
        task_frame = QFrame()
        task_frame.setObjectName("taskFrame")
        task_frame.setStyleSheet(f"""
            QFrame#taskFrame {{
                background-color: {PANEL_BG_COLOR};
                border-radius: {FRAME_CORNER_RADIUS}px;
                padding: {SPACING_SM}px;
            }}
        """)
        task_layout = QVBoxLayout(task_frame)
        task_layout.setSpacing(SPACING_SM)
        task_layout.setContentsMargins(SPACING_SM, SPACING_SM, SPACING_SM, SPACING_SM)

        # Title label
        title_label = QLabel("Current Task:")
        title_font = QFont()
        title_font.setPointSize(11)
        title_font.setBold(True)
        title_label.setFont(title_font)
        title_label.setStyleSheet(f"color: {PRIMARY_TEXT};")

        # Task content label (sized for 3 lines of text) - use full width
        self.task_label = QLabel("Starting...")
        task_font = QFont()
        task_font.setPointSize(10)
        self.task_label.setFont(task_font)
        self.task_label.setStyleSheet(f"color: {PRIMARY_TEXT};")
        self.task_label.setWordWrap(True)
        self.task_label.setMinimumHeight(60)
        self.task_label.setSizePolicy(
            self.task_label.sizePolicy().horizontalPolicy(),
            self.task_label.sizePolicy().verticalPolicy()
        )

        task_layout.addWidget(title_label)
        task_layout.addWidget(self.task_label)  # No alignment constraint - use full width

        parent_layout.addWidget(task_frame)

    def _create_stats_section(self, parent_layout: QVBoxLayout):
        """Create statistics display section."""
        stats_frame = QFrame()
        stats_frame.setObjectName("statsFrame")
        stats_frame.setStyleSheet(f"""
            QFrame#statsFrame {{
                background-color: {PANEL_BG_COLOR};
                border-radius: {FRAME_CORNER_RADIUS}px;
                padding: {SPACING_SM}px;
            }}
        """)

        stats_layout = QVBoxLayout(stats_frame)
        stats_layout.setSpacing(SPACING_SM)

        # Title
        stats_title = QLabel("📊 Statistics")
        stats_title_font = QFont()
        stats_title_font.setPointSize(10)
        stats_title_font.setBold(True)
        stats_title.setFont(stats_title_font)
        stats_title.setStyleSheet(f"color: {PRIMARY_TEXT};")

        # Stats content
        self.stats_label = QLabel("Starting...")
        stats_font = QFont()
        stats_font.setPointSize(9)
        self.stats_label.setFont(stats_font)
        self.stats_label.setStyleSheet(f"color: {PRIMARY_TEXT};")
        self.stats_label.setWordWrap(True)

        stats_layout.addWidget(stats_title)
        stats_layout.addWidget(self.stats_label)

        parent_layout.addWidget(stats_frame)

    def _create_log_section(self, parent_layout: QVBoxLayout):
        """Create detailed log section."""
        log_frame = QFrame()
        log_frame.setObjectName("logFrame")
        log_frame.setStyleSheet(f"""
            QFrame#logFrame {{
                background-color: {PANEL_BG_COLOR};
                border-radius: {FRAME_CORNER_RADIUS}px;
                padding: {SPACING_SM}px;
            }}
        """)

        log_layout = QVBoxLayout(log_frame)
        log_layout.setSpacing(SPACING_SM)

        # Title
        log_title = QLabel("📋 Detailed Log")
        log_title_font = QFont()
        log_title_font.setPointSize(10)
        log_title_font.setBold(True)
        log_title.setFont(log_title_font)
        log_title.setStyleSheet(f"color: {PRIMARY_TEXT};")

        # Log text area
        self.log_text = QTextEdit()
        self.log_text.setReadOnly(True)
        self.log_text.setMinimumHeight(150)

        # Monospace font for log
        log_font = QFont("Courier")
        log_font.setPointSize(9)
        self.log_text.setFont(log_font)

        self.log_text.setStyleSheet(f"""
            QTextEdit {{
                background-color: {FRAME_BG_COLOR};
                color: {PRIMARY_TEXT};
                border: 1px solid {PANEL_BG_COLOR};
                border-radius: 4px;
            }}
        """)

        log_layout.addWidget(log_title)
        log_layout.addWidget(self.log_text, 1)  # Stretch factor 1

        parent_layout.addWidget(log_frame, 1)  # Stretch factor 1

    def _create_button_section(self, parent_layout: QVBoxLayout):
        """Create control buttons section."""
        button_layout = QHBoxLayout()
        button_layout.setSpacing(SPACING_SM)

        # Left side buttons (Pause, Cancel)
        self.pause_button = QPushButton("⏸️ Pause")
        self.pause_button.setMinimumWidth(100)
        self.pause_button.clicked.connect(self._on_pause_clicked)
        self.pause_button.setStyleSheet(f"""
            QPushButton {{
                background-color: {WARNING_ORANGE};
                color: white;
                border: none;
                border-radius: 10px;
                padding: 6px 12px;
                font-size: 10pt;
            }}
            QPushButton:hover {{
                background-color: #DD9043;
            }}
        """)

        self.cancel_button = QPushButton("❌ Cancel")
        self.cancel_button.setMinimumWidth(100)
        self.cancel_button.clicked.connect(self._on_cancel_clicked)
        self.cancel_button.setStyleSheet(f"""
            QPushButton {{
                background-color: {DESTRUCTIVE_BUTTON_BG};
                color: white;
                border: none;
                border-radius: 10px;
                padding: 6px 12px;
                font-size: 10pt;
            }}
            QPushButton:hover {{
                background-color: {DESTRUCTIVE_BUTTON_HOVER};
            }}
        """)

        button_layout.addWidget(self.pause_button)
        button_layout.addWidget(self.cancel_button)
        button_layout.addStretch()

        # Right side buttons (Save Log, Close)
        self.save_log_button = QPushButton("💾 Save Log")
        self.save_log_button.setMinimumWidth(120)
        self.save_log_button.setEnabled(False)
        self.save_log_button.clicked.connect(self._on_save_log_clicked)
        self.save_log_button.setStyleSheet(f"""
            QPushButton {{
                background-color: {SECONDARY_BUTTON_BG};
                color: {PRIMARY_TEXT};
                border: 1px solid #C8C8CD;
                border-radius: 10px;
                padding: 6px 12px;
                font-size: 10pt;
            }}
            QPushButton:hover {{
                background-color: {SECONDARY_BUTTON_HOVER};
            }}
            QPushButton:disabled {{
                background-color: #F0F0F0;
                color: #C0C0C0;
            }}
        """)

        self.close_button = QPushButton("✅ Close")
        self.close_button.setMinimumWidth(100)
        self.close_button.setEnabled(False)
        self.close_button.clicked.connect(self.accept)
        self.close_button.setStyleSheet(f"""
            QPushButton {{
                background-color: {SUCCESS_GREEN};
                color: white;
                border: none;
                border-radius: 10px;
                padding: 6px 12px;
                font-size: 10pt;
            }}
            QPushButton:hover {{
                background-color: #2AA64A;
            }}
            QPushButton:disabled {{
                background-color: #A0D0A0;
                color: white;
            }}
        """)

        button_layout.addWidget(self.save_log_button)
        button_layout.addWidget(self.close_button)

        parent_layout.addLayout(button_layout)

    @Slot(int, str)
    def update_progress(self, value: int, status: str = ""):
        """
        Update progress bar and status (THREAD-SAFE).

        This slot can be connected to worker thread signals.

        Args:
            value: Progress value (0-100)
            status: Optional status message
        """
        # Update progress bar
        self.progress_bar.setValue(value)

        # Update status if provided
        if status:
            self.task_label.setText(status)

        # Update statistics
        self._update_statistics(value)

    @Slot(str, bool)
    def add_log_message(self, message: str, is_error: bool = False, count_in_stats: bool = True):
        """
        Add a message to the log (THREAD-SAFE).

        This slot can be connected to worker thread signals.

        Args:
            message: Log message
            is_error: Whether this is an error message
            count_in_stats: Whether to count this message in success/failure statistics
                           Set to False for optimization/calibration steps
        """
        if not self._show_details:
            return

        # Format with timestamp
        timestamp = datetime.now().strftime('%H:%M:%S')
        icon = "❌" if is_error else "✅"
        log_entry = f"[{timestamp}] {icon} {message}"

        # Append to log
        self.log_text.append(log_entry)

        # Scroll to bottom
        self.log_text.verticalScrollBar().setValue(
            self.log_text.verticalScrollBar().maximum()
        )

        # Track success/failure only if counting is enabled
        if count_in_stats:
            if is_error:
                self.failed_tasks += 1
            else:
                self.completed_tasks += 1

    def reset_statistics(self, reset_timer: bool = False):
        """
        Reset success/failure counters.

        Call this before starting the main fitting phase to exclude
        optimization/calibration steps from the final statistics.

        Args:
            reset_timer: If True, also reset the elapsed time counter.
                        Default False to keep cumulative timing from dialog start.
        """
        self.completed_tasks = 0
        self.failed_tasks = 0
        if reset_timer:
            self.start_time = datetime.now()

    @Slot(str)
    def complete(self, message: str = "Processing completed"):
        """
        Mark processing as completed (THREAD-SAFE).

        This slot can be connected to worker thread signals.

        Args:
            message: Completion message
        """
        self._is_completed = True

        # Update progress to 100%
        self.progress_bar.setValue(100)
        self.task_label.setText(message)

        # Final statistics
        total_time = datetime.now() - self.start_time
        total_time_str = str(total_time).split('.')[0]  # Remove microseconds

        final_stats = (
            f"🏁 Completed in: {total_time_str}\n"
            f"✅ Successful: {self.completed_tasks} | ❌ Failed: {self.failed_tasks}"
        )
        self.stats_label.setText(final_stats)

        # Enable/disable buttons
        self.cancel_button.setEnabled(False)
        self.pause_button.setEnabled(False)
        self.close_button.setEnabled(True)

        if self._show_details:
            self.save_log_button.setEnabled(True)

    def _update_statistics(self, progress_value: int):
        """Update statistics display."""
        elapsed = datetime.now() - self.start_time
        elapsed_str = str(elapsed).split('.')[0]  # Remove microseconds

        # Calculate ETA
        if progress_value > 0:
            eta_total = elapsed.total_seconds() * (100 / progress_value)
            eta_remaining = eta_total - elapsed.total_seconds()
            eta_str = str(timedelta(seconds=int(eta_remaining)))
        else:
            eta_str = "Calculating..."

        # Format statistics
        stats_text = (
            f"⏱️  Elapsed: {elapsed_str} | 🔮 ETA: {eta_str}\n"
            f"✅ Completed: {self.completed_tasks} | ❌ Failed: {self.failed_tasks}"
        )
        self.stats_label.setText(stats_text)

    def _on_pause_clicked(self):
        """Handle pause button click."""
        self._is_paused = not self._is_paused

        if self._is_paused:
            self.pause_button.setText("▶️ Resume")
            self.task_label.setText("⏸️ Processing paused...")
        else:
            self.pause_button.setText("⏸️ Pause")

        # Emit signal for worker thread to handle
        self.paused.emit(self._is_paused)

    def _on_cancel_clicked(self):
        """Handle cancel button click."""
        self._is_cancelled = True

        # Emit signal for worker thread
        self.cancelled.emit()

        # Close dialog
        self.reject()

    def _on_save_log_clicked(self):
        """Handle save log button click."""
        if not self._show_details:
            return

        # Open file dialog
        filename, _ = QFileDialog.getSaveFileName(
            self,
            "Save Processing Log",
            str(Path.home() / "processing_log.txt"),
            "Text files (*.txt);;All files (*.*)"
        )

        if filename:
            try:
                with open(filename, 'w', encoding='utf-8') as f:
                    f.write(self.log_text.toPlainText())

                QMessageBox.information(
                    self,
                    "Success",
                    f"Log saved to:\n{filename}"
                )
            except Exception as e:
                QMessageBox.critical(
                    self,
                    "Error",
                    f"Failed to save log:\n{str(e)}"
                )

    # Public properties for external access
    @property
    def is_cancelled(self) -> bool:
        """Check if operation was cancelled."""
        return self._is_cancelled

    @property
    def is_completed(self) -> bool:
        """Check if operation completed."""
        return self._is_completed

    @property
    def is_paused(self) -> bool:
        """Check if operation is paused."""
        return self._is_paused

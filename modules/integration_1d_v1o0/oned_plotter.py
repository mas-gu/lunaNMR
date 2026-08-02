# ABOUTME: Interactive 1D spectrum plot - zoom, axis scaling, box select and click select.
# ABOUTME: Emits selected peaks to the page; holds no measurement logic of its own.

from PySide6.QtCore import Signal
from PySide6.QtWidgets import (QCheckBox, QComboBox, QDoubleSpinBox, QHBoxLayout,
                               QLabel, QPushButton, QVBoxLayout, QWidget)

from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
from matplotlib.backends.backend_qtagg import NavigationToolbar2QT
from matplotlib.figure import Figure
from matplotlib.widgets import RectangleSelector

# Vertical zoom range. A hydrolysis series spans orders of magnitude between the
# growing and decaying peak, so the y-axis needs far more than a full-scale view.
Y_SCALE_MIN = 0.001
Y_SCALE_MAX = 1000.0

PEAK_COLOUR = '#d62728'
SELECTED_COLOUR = '#1f77b4'


class OneDPlotter(QWidget):
    """Spectrum trace with interactive peak selection.

    Two ways to pick a peak, both requested by the user:
      - drag a rectangle over a region, everything detectable inside becomes a peak
      - click a position directly, which snaps to the local maximum
    """

    peaks_selected = Signal(list)      # [{'ppm', 'height', ...}] from a box drag
    position_clicked = Signal(float)   # a bare ppm from a click

    def __init__(self, parent=None):
        super().__init__(parent)

        self.spectrum = None
        self.peaks = []
        self._select_mode = 'box'

        self.figure = Figure(figsize=(10, 4.5))
        self.canvas = FigureCanvasQTAgg(self.figure)
        self.axes = self.figure.add_subplot(111)
        self.toolbar = NavigationToolbar2QT(self.canvas, self)

        self._build_controls()

        layout = QVBoxLayout(self)
        layout.addWidget(self.toolbar)
        layout.addLayout(self.controls)
        layout.addWidget(self.canvas)

        self._selector = RectangleSelector(
            self.axes, self._on_box, useblit=True, button=[1],
            minspanx=1e-6, spancoords='data', interactive=False)
        self.canvas.mpl_connect('button_press_event', self._on_click)

        self._apply_mode()

    def _build_controls(self):
        self.controls = QHBoxLayout()

        self.mode_box = QComboBox()
        self.mode_box.addItems(['Box select (auto-detect)', 'Click position'])
        self.mode_box.currentIndexChanged.connect(self._apply_mode)

        self.y_scale = QDoubleSpinBox()
        self.y_scale.setRange(Y_SCALE_MIN, Y_SCALE_MAX)
        self.y_scale.setDecimals(3)
        self.y_scale.setValue(1.0)
        self.y_scale.setSingleStep(0.1)
        self.y_scale.valueChanged.connect(self.refresh)

        self.show_peaks = QCheckBox("Show peaks")
        self.show_peaks.setChecked(True)
        self.show_peaks.stateChanged.connect(self.refresh)

        self.reset_button = QPushButton("Reset view")
        self.reset_button.clicked.connect(self.reset_view)

        self.controls.addWidget(QLabel("Mode:"))
        self.controls.addWidget(self.mode_box)
        self.controls.addSpacing(12)
        self.controls.addWidget(QLabel("Y scale:"))
        self.controls.addWidget(self.y_scale)
        self.controls.addSpacing(12)
        self.controls.addWidget(self.show_peaks)
        self.controls.addStretch()
        self.controls.addWidget(self.reset_button)

    def _apply_mode(self):
        self._select_mode = 'box' if self.mode_box.currentIndex() == 0 else 'click'
        self._selector.set_active(self._select_mode == 'box')

    def set_spectrum(self, spectrum, keep_view=False):
        bounds = (self.axes.get_xlim(), self.axes.get_ylim()) if keep_view else None
        self.spectrum = spectrum
        self.refresh(bounds=bounds)

    def set_peaks(self, peaks):
        self.peaks = list(peaks or [])
        self.refresh(keep_view=True)

    def reset_view(self):
        self.y_scale.setValue(1.0)
        self.refresh()

    def refresh(self, *args, bounds=None, keep_view=False):
        if keep_view and self.spectrum is not None:
            bounds = (self.axes.get_xlim(), self.axes.get_ylim())

        self.axes.clear()

        if self.spectrum is None:
            self.canvas.draw_idle()
            return

        ppm, data = self.spectrum.ppm_axis, self.spectrum.data
        self.axes.plot(ppm, data, lw=0.6, color='#333333')

        if self.show_peaks.isChecked():
            for peak in self.peaks:
                colour = SELECTED_COLOUR if peak.get('selected') else PEAK_COLOUR
                self.axes.axvline(peak['ppm'], color=colour, lw=0.8, alpha=0.7)
                label = peak.get('assignment') or f"{peak['ppm']:.4f}"
                self.axes.annotate(label, (peak['ppm'], peak.get('height', 0)),
                                   textcoords='offset points', xytext=(0, 6),
                                   ha='center', fontsize=7, color=colour)

        self.axes.set_xlabel('ppm')
        self.axes.invert_xaxis()

        if bounds is not None:
            self.axes.set_xlim(bounds[0])
            self.axes.set_ylim(bounds[1])
        else:
            span = float(data.max()) or 1.0
            self.axes.set_ylim(-0.1 * span / self.y_scale.value(),
                               span / self.y_scale.value())

        self.figure.tight_layout()
        self.canvas.draw_idle()

    def _on_box(self, press, release):
        if self.spectrum is None or press.xdata is None or release.xdata is None:
            return
        if abs(release.xdata - press.xdata) <= 0:
            return

        from oned_series import peaks_in_range

        found = peaks_in_range(self.spectrum, press.xdata, release.xdata)
        if found:
            self.peaks_selected.emit(found)

    def _on_click(self, event):
        if (self._select_mode != 'click' or self.spectrum is None
                or event.inaxes is not self.axes or event.button != 1
                or event.xdata is None or self.toolbar.mode):
            return

        self.position_clicked.emit(float(event.xdata))

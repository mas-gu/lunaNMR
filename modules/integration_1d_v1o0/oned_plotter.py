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

# One scroll notch shrinks the visible span to this fraction of itself.
ZOOM_FACTOR = 0.8

# Movement below this many pixels between press and release counts as a click, not a
# drag, so a pick with a twitchy hand still registers.
DRAG_THRESHOLD_PIXELS = 3


class OneDPlotter(QWidget):
    """Spectrum trace with interactive peak selection.

    Two ways to pick a peak:
      - drag a rectangle over a region, everything detectable inside becomes a peak
      - click a position directly, which snaps to the local maximum

    Navigation follows the convention used elsewhere in lunaNMR: a drag pans, a press
    and release in the same place is a pick. Scroll zooms about the pointer, and the
    middle button pans or picks in either mode so box selection keeps the left drag.
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

        self._pan = None
        self.canvas.mpl_connect('scroll_event', self._on_scroll)
        self.canvas.mpl_connect('button_press_event', self._on_press)
        self.canvas.mpl_connect('motion_notify_event', self._on_motion)
        self.canvas.mpl_connect('button_release_event', self._on_release)

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
        self.y_scale.valueChanged.connect(self._apply_y_scale)

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

        hint = QLabel("scroll: zoom  ·  scroll on Y axis: scale trace  ·  "
                      "middle-drag: pan  ·  middle-click: pick peak")
        hint.setStyleSheet("color: #666; font-size: 10px;")
        self.controls.addWidget(hint)
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

    # ------------------------------------------------------------ navigation

    def _on_y_axis(self, event):
        """True when the pointer is over the y-axis strip beside the plotting area."""
        if event.x is None or event.y is None:
            return False

        box = self.axes.get_window_extent()

        return event.x < box.x0 and box.y0 <= event.y <= box.y1

    def _on_scroll(self, event):
        """Scroll zooms about the pointer; shift-scroll zooms the intensity axis.

        Scrolling on the y axis itself rescales the spectrum amplitude instead, which
        is the control the y-scale box holds - so the two stay in step rather than
        offering two independent notions of vertical size.
        """
        if self.spectrum is None:
            return

        direction = getattr(event, 'button', None)
        step = getattr(event, 'step', 0)
        zooming_in = direction == 'up' or step > 0
        factor = ZOOM_FACTOR if zooming_in else 1.0 / ZOOM_FACTOR

        if event.inaxes is not self.axes:
            if self._on_y_axis(event):
                # a larger y scale draws the trace taller, so zooming in raises it
                self._scale_amplitude(1.0 / factor)
            return

        if event.key and 'shift' in str(event.key):
            self._zoom_axis(self.axes.get_ylim(), self.axes.set_ylim,
                            getattr(event, 'ydata', None), factor)
        else:
            if event.xdata is None:
                return
            self._zoom_axis(self.axes.get_xlim(), self.axes.set_xlim,
                            event.xdata, factor)

        self.canvas.draw_idle()

    def _apply_y_scale(self, *args):
        """Set the intensity limits from the current scale, leaving the ppm view alone.

        Redrawing the whole plot instead would reset the ppm zoom, so rescaling the
        trace would throw away the region the user had navigated to.
        """
        if self.spectrum is None:
            return

        span = float(self.spectrum.data.max()) or 1.0
        scale = self.y_scale.value()
        self.axes.set_ylim(-0.1 * span / scale, span / scale)
        self.canvas.draw_idle()

    def _scale_amplitude(self, multiplier):
        """Multiply the spectrum's vertical scale, clamped to the spin box range."""
        target = self.y_scale.value() * multiplier

        self.y_scale.setValue(min(max(target, self.y_scale.minimum()),
                                  self.y_scale.maximum()))

    @staticmethod
    def _zoom_axis(limits, apply_limits, anchor, factor):
        """Scale `limits` by `factor` while holding `anchor` at the same screen fraction."""
        low, high = limits
        if anchor is None:
            anchor = (low + high) / 2.0

        apply_limits(anchor + (low - anchor) * factor,
                     anchor + (high - anchor) * factor)

    def _pan_allowed(self, button):
        """Middle button always pans. The left button pans only when it is not
        reserved for drawing the box-selection rectangle."""
        if button == 2:
            return True
        return button == 1 and self._select_mode == 'click'

    def _on_press(self, event):
        if (event.inaxes is not self.axes or self.spectrum is None
                or self.toolbar.mode or event.button not in (1, 2)):
            return

        self._pan = {
            'button': event.button,
            'x': event.x,
            'xlim': self.axes.get_xlim(),
            'ylim': self.axes.get_ylim(),
            'y': event.y,
            'xdata': event.xdata,
            'moved': False,
        }

    def _on_motion(self, event):
        if not self._pan or event.inaxes is not self.axes:
            return

        dx_pixels = event.x - self._pan['x']
        dy_pixels = event.y - self._pan['y']

        if abs(dx_pixels) > DRAG_THRESHOLD_PIXELS or abs(dy_pixels) > DRAG_THRESHOLD_PIXELS:
            self._pan['moved'] = True

        if not (self._pan['moved'] and self._pan_allowed(self._pan['button'])):
            return

        inverse = self.axes.transData.inverted()
        origin_x, origin_y = inverse.transform((0.0, 0.0))
        shifted_x, shifted_y = inverse.transform((dx_pixels, dy_pixels))

        low, high = self._pan['xlim']
        self.axes.set_xlim(low - (shifted_x - origin_x), high - (shifted_x - origin_x))
        low, high = self._pan['ylim']
        self.axes.set_ylim(low - (shifted_y - origin_y), high - (shifted_y - origin_y))

        self.canvas.draw_idle()

    def _on_release(self, event):
        pan, self._pan = self._pan, None

        if not pan or pan['moved'] or event.xdata is None:
            return

        # A press and release in the same place is a pick, not a pan: the middle
        # button picks in either mode, the left button only in click mode.
        if pan['button'] == 2 or self._select_mode == 'click':
            self.position_clicked.emit(float(event.xdata))

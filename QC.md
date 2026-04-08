# Series QC Module - Complete Implementation Specification

## Overview

Create a QC module accessible under "Modules" menu that analyzes series integration results from `training_data.json`. Supports single-series analysis and multi-series comparison.

---

## Location & Access

- **Menu**: Modules → Series QC (after DynamiXs entry)
- **Dialog**: Non-modal, reusable pattern from DynamiXs
- **Output**: Auto-saves `QC_report.json` in same directory as `training_data.json`

---

## Files to Create

### 1. `lunaNMR/utils/series_qc_analyzer.py`

Core analysis module with two classes:

```python
# ABOUTME: Analyzes series integration results for quality control.
# ABOUTME: Computes per-peak and per-spectrum statistics, flags anomalies.

from dataclasses import dataclass, field
from typing import List, Dict, Tuple, Optional
from pathlib import Path
import json
import numpy as np

@dataclass
class PeakQCStats:
    """QC statistics for a single peak across all spectra in series."""
    peak_id: str
    assignment: str  # SPARKY assignment (e.g., "A.2.ASP.H")

    # Time series (per-spectrum values, indexed by spectrum order)
    volumes: List[Optional[float]] = field(default_factory=list)
    heights: List[Optional[float]] = field(default_factory=list)
    r2_values: List[Optional[float]] = field(default_factory=list)
    lw_f1_values: List[Optional[float]] = field(default_factory=list)  # 15N or 13C
    lw_f2_values: List[Optional[float]] = field(default_factory=list)  # 1H
    cluster_sizes: List[int] = field(default_factory=list)
    pos_f1_values: List[Optional[float]] = field(default_factory=list)
    pos_f2_values: List[Optional[float]] = field(default_factory=list)

    # Aggregates (computed after loading)
    volume_cv: float = 0.0
    height_cv: float = 0.0
    r2_mean: float = 0.0
    r2_min: float = 0.0
    lw_f1_cv: float = 0.0
    lw_f2_cv: float = 0.0
    pos_f1_drift: float = 0.0  # max - min
    pos_f2_drift: float = 0.0

    # Cluster stability
    is_cluster_stable: bool = True
    cluster_transitions: List[Tuple[int, int, int]] = field(default_factory=list)

    # Outlier detection
    outlier_spectra: List[int] = field(default_factory=list)

    # Flags
    flags: List[str] = field(default_factory=list)


@dataclass
class SeriesQCReport:
    """Complete QC report for a series integration run."""
    source_path: str
    n_peaks: int
    n_spectra: int
    spectrum_names: List[str] = field(default_factory=list)

    # Processing context (affects expected CV values)
    processing_mode: str = "unknown"  # "parallel" or "sequential"
    nucleus_type: str = "unknown"     # "15N", "13C", etc.
    fix_positions_used: bool = False
    fix_linewidths_used: bool = False

    # Aggregates
    volume_cv_median: float = 0.0
    volume_cv_mean: float = 0.0
    r2_mean: float = 0.0
    r2_median: float = 0.0
    n_cluster_stable: int = 0
    n_lw_stable: int = 0
    fit_success_rate: float = 0.0

    # Per-spectrum quality (identify bad spectra)
    per_spectrum_r2_mean: List[float] = field(default_factory=list)
    per_spectrum_n_failed: List[int] = field(default_factory=list)
    outlier_spectra: List[int] = field(default_factory=list)

    # Per-peak
    peak_stats: Dict[str, PeakQCStats] = field(default_factory=dict)
    flagged_peaks: List[str] = field(default_factory=list)

    # Special categories
    dummy_peak_ids: List[str] = field(default_factory=list)
    failed_peak_ids: List[str] = field(default_factory=list)


class SeriesQCAnalyzer:
    """Analyzes training_data.json for quality control metrics."""

    # Default thresholds (may be overridden by quality_categories.py)
    VOLUME_CV_WARN = 0.3
    R2_POOR = 0.7
    LW_CV_WARN = 0.3
    LW_JUMP_THRESHOLD = 0.5  # 50% change between consecutive spectra

    def __init__(self):
        self.data = None
        self.report = None

    def load_training_data(self, path: str) -> bool:
        """Load and parse training_data.json."""
        pass

    def compute_per_peak_stats(self) -> None:
        """Compute CV, means, and flags for each peak."""
        pass

    def analyze_cluster_stability(self) -> None:
        """Track cluster size changes across spectra."""
        pass

    def detect_anomalies(self) -> None:
        """Detect LW jumps, boundary hits, outliers."""
        pass

    def identify_outlier_spectra(self) -> None:
        """Find spectra with many problems (low R², many failures)."""
        pass

    def generate_report(self) -> SeriesQCReport:
        """Generate complete QC report."""
        pass

    def save_report(self, output_dir: str) -> str:
        """Save QC_report.json to output directory. Returns path."""
        pass


class SeriesComparator:
    """Compares multiple series integration runs."""

    def __init__(self):
        self.series: Dict[str, SeriesQCReport] = {}

    def add_series(self, name: str, training_data_path: str) -> None:
        """Load and analyze a series."""
        pass

    def get_common_peaks(self) -> List[str]:
        """Get peak IDs present in ALL series (intersection)."""
        pass

    def get_series_only_peaks(self, series_name: str) -> List[str]:
        """Get peaks only in this series (not in others)."""
        pass

    def compare(self) -> Dict:
        """Generate comparison report with winner per metric."""
        pass

    def get_per_peak_winner(self) -> Dict[str, str]:
        """For each peak, which series gave best result."""
        pass
```

---

### 2. `lunaNMR/gui/dialogs/series_qc_dialog.py`

```python
# ABOUTME: GUI dialog for Series QC analysis and comparison.
# ABOUTME: Provides drag-drop interface and visual QC reports.

from PySide6.QtWidgets import (
    QDialog, QVBoxLayout, QHBoxLayout, QLabel, QPushButton,
    QTableWidget, QTableWidgetItem, QProgressBar, QListWidget,
    QAbstractItemView, QFrame, QSplitter, QTabWidget, QGroupBox
)
from PySide6.QtCore import Qt, Signal, QMimeData
from PySide6.QtGui import QDrag

class DraggableSeriesList(QListWidget):
    """List widget that supports dragging series results."""
    # Copy pattern from modules/dynamiXs_v2o0/dynamiXs_gui.py:463-501
    pass

class DropTargetLabel(QLabel):
    """Label that accepts dropped series results."""
    series_dropped = Signal(str, str, str)  # field_name, series_name, path
    # Copy pattern from modules/dynamiXs_v2o0/dynamiXs_gui.py:504-559
    pass

class SeriesQCDialog(QDialog):
    """Main QC dialog with drag-drop, analysis, and comparison."""

    def __init__(self, parent=None, main_window=None):
        super().__init__(parent)
        self.main_window = main_window
        self.setWindowTitle("Series QC Analysis")
        self.setMinimumSize(900, 700)
        self.setup_ui()

    def setup_ui(self):
        """Build the UI with 3 main sections."""
        layout = QVBoxLayout(self)

        # Section 1: Drop zone for series selection
        self._create_drop_zone(layout)

        # Section 2: Action buttons
        self._create_action_buttons(layout)

        # Section 3: Results area (tabs for single/compare)
        self._create_results_area(layout)

    def _create_drop_zone(self, parent_layout):
        """Create drag-drop area for series results."""
        # Use QListWidget with drag-drop enabled
        # Show: folder name, spectrum count, [x] remove button
        pass

    def _create_action_buttons(self, parent_layout):
        """Create Browse / Analyze Single / Compare Selected buttons."""
        pass

    def browse_for_series(self):
        """Open folder dialog to select series results.

        Browse Behavior:
        - Opens folder selection dialog (multi-select enabled)
        - For each selected folder, searches for training_data.json:
          1. Check folder itself
          2. Check data_collection/ subfolder
          3. Check any subfolder containing training_data.json
        - Immediately validates each found file (valid JSON, has required fields)
        - Invalid selections show error message and are not added
        - Valid selections are added to the series list
        """
        pass

    def _create_results_area(self, parent_layout):
        """Create tabbed results display."""
        # Tab 1: Summary Dashboard
        # Tab 2: Flagged Peaks Table
        # Tab 3: Comparison (when 2+ series)
        # Tab 4: Distribution Plots (optional)
        pass

    def _create_summary_dashboard(self) -> QWidget:
        """Create summary panel with progress bars."""
        # Show:
        # - Peaks: N, Spectra: M, Fit Success: X%
        # - Volume CV: median [progress bar] Good/Warning
        # - R² Quality: mean [progress bar] Good/Warning
        # - Cluster Stable: N/M [progress bar] Good/Warning
        # - LW Stable: N/M [progress bar] Good/Warning
        # - Flagged peaks count
        pass

    def _create_flagged_table(self) -> QTableWidget:
        """Create sortable table of flagged peaks."""
        # Columns: Peak, Issue, Volume CV, R² mean, LW CV, Details
        pass

    def _create_comparison_table(self) -> QTableWidget:
        """Create comparison table for multiple series."""
        # Columns: Metric, Series A, Series B, ..., Winner
        pass
```

---

### 3. Modifications to Existing Files

#### `lunaNMR/gui/main_window.py` (~line 4767)

Add after DynamiXs menu entry:

```python
# In menu creation section (around line 4767):
series_qc_action = QAction("Series QC Analysis", self)
series_qc_action.triggered.connect(self.launch_series_qc)
modules_menu.addAction(series_qc_action)

# Add method (after launch_dynamixs):
def launch_series_qc(self):
    """Launch Series QC Analysis dialog."""
    from lunaNMR.gui.dialogs import SeriesQCDialog

    if hasattr(self, 'series_qc_dialog') and self.series_qc_dialog and self.series_qc_dialog.isVisible():
        self.series_qc_dialog.raise_()
        self.series_qc_dialog.activateWindow()
        return

    self.series_qc_dialog = SeriesQCDialog(parent=self, main_window=self)
    self.series_qc_dialog.show()
```

#### `lunaNMR/gui/dialogs/__init__.py`

Add export:
```python
from .series_qc_dialog import SeriesQCDialog
```

---

## training_data.json Field Mapping

### Per-Peak Fields to Extract:

| Field in JSON | Use in QC | Description |
|---------------|-----------|-------------|
| `peak_id` | Identification | Peak identifier |
| `assignment` | Display | SPARKY assignment |
| `volume` | Volume CV | Integrated volume |
| `height` | Height CV | Peak height |
| `r_squared` | Fit quality | R² of fit |
| `lw_total_f1` | LW stability | Total linewidth F1 |
| `lw_total_f2` | LW stability | Total linewidth F2 |
| `pos_f1` | Position drift | Fitted position F1 |
| `pos_f2` | Position drift | Fitted position F2 |
| `cluster_size` | Cluster stability | Size of cluster (1=isolated) |
| `is_clustered` | Cluster flag | Boolean |
| `success` | Fit success | Boolean |
| `was_detected` | Dummy detection | False = dummy peak |
| `detection_quality` | Dummy detection | "Missed" = dummy |
| `fix_positions_applied` | Context | Positions were fixed |
| `fix_linewidths_applied` | Context | Linewidths were fixed |

### Per-Spectrum Fields to Extract:

| Field in JSON | Use in QC | Description |
|---------------|-----------|-------------|
| `spectrum_name` | Identification | Filename |
| `processing_mode` | Context | "parallel" or "sequential" |
| `nucleus_type` | Thresholds | "15N", "13C", etc. |
| `peak_count` | Summary | Total peaks fitted |
| `is_reference_spectrum` | Context | Boolean |

---

## Flagging Criteria

### Default Thresholds:

| Flag | Condition | Severity |
|------|-----------|----------|
| `high_vol_cv` | Volume CV > 0.3 | Warning |
| `poor_r2` | R² mean < 0.7 | Warning |
| `very_poor_r2` | R² mean < 0.5 | Error |
| `lw_unstable` | LW CV > 0.3 | Warning |
| `lw_jump` | LW change > 50% between consecutive spectra | Warning |
| `cluster_changed` | cluster_size varies across spectra | Info |
| `position_drift` | pos drift > 0.02 ppm (1H) or > 0.2 ppm (15N) | Warning |
| `fit_failed` | success=False in any spectrum | Error |
| `dummy_peak` | was_detected=False | Info |

### Context-Aware Adjustments:

- If `fix_linewidths_applied=True`: LW CV should be ~0, flag if > 0.05
- If `fix_positions_applied=True`: Position drift should be ~0
- Nucleus-specific: 13C may have different typical LW range than 15N

---

## GUI Layout Mockups

### Main Dialog:

```
┌─────────────────────────────────────────────────────────────────┐
│  Series QC Analysis                                        [X]  │
├─────────────────────────────────────────────────────────────────┤
│                                                                 │
│  ┌─────────────────────────────────────────────────────────┐   │
│  │  Drop series results folders here                        │   │
│  │  (or click Browse)                                       │   │
│  │                                                          │   │
│  │  📁 series_results_param_A       12 spectra    [x]      │   │
│  │  📁 series_results_param_B       12 spectra    [x]      │   │
│  │                                                          │   │
│  └─────────────────────────────────────────────────────────┘   │
│                                                                 │
│  [Browse...]  [Analyze Single]  [Compare Selected]             │
│                                                                 │
├─────────────────────────────────────────────────────────────────┤
│  ┌─────────┬─────────────┬────────────┬────────────┐           │
│  │ Summary │ Flagged (4) │ Comparison │ Plots      │           │
│  └─────────┴─────────────┴────────────┴────────────┘           │
│                                                                 │
│  ┌─────────────────────────────────────────────────────────┐   │
│  │  SERIES QC SUMMARY                                       │   │
│  │  ─────────────────                                       │   │
│  │  Peaks: 45    Spectra: 12    Fit Success: 98.2%         │   │
│  │                                                          │   │
│  │  Volume CV:     0.12 median  [████████░░] Good          │   │
│  │  R² Quality:    0.89 mean    [█████████░] Good          │   │
│  │  Cluster Stable: 43/45       [█████████░] Good          │   │
│  │  LW Stable:     41/45        [████████░░] Warning       │   │
│  │                                                          │   │
│  │  ⚠ 4 peaks flagged for review                           │   │
│  │                                                          │   │
│  │  Outlier Spectra: spectrum_05.ft (low R²)               │   │
│  └─────────────────────────────────────────────────────────┘   │
│                                                                 │
└─────────────────────────────────────────────────────────────────┘
```

### Flagged Peaks Tab:

```
┌─────────────────────────────────────────────────────────────────┐
│ Peak       │ Issue        │ Vol CV │ R² mean │ LW CV │ Details │
├────────────┼──────────────┼────────┼─────────┼───────┼─────────┤
│ Peak_023   │ high_vol_cv  │ 0.42   │ 0.85    │ 0.15  │ size=3  │
│ Peak_067   │ lw_jump      │ 0.18   │ 0.91    │ 0.55  │ @spec 5 │
│ Peak_089   │ poor_r2      │ 0.25   │ 0.52    │ 0.22  │ 4 fails │
│ Peak_112   │ cluster_chg  │ 0.31   │ 0.78    │ 0.28  │ 1→3→1   │
└─────────────────────────────────────────────────────────────────┘
```

### Comparison Tab:

```
┌─────────────────────────────────────────────────────────────────┐
│ Metric             │ param_A  │ param_B  │ Winner │
├────────────────────┼──────────┼──────────┼────────┤
│ Volume CV (median) │ 0.15     │ 0.12     │ B      │
│ R² (mean)          │ 0.87     │ 0.91     │ B      │
│ Cluster Stable     │ 40/45    │ 43/45    │ B      │
│ Flagged Peaks      │ 6        │ 4        │ B      │
│ Fit Success Rate   │ 96.5%    │ 98.2%    │ B      │
├────────────────────┴──────────┴──────────┴────────┤
│ Common peaks: 43   │   A only: 2   │   B only: 0  │
└─────────────────────────────────────────────────────────────────┘
```

---

## Edge Cases to Handle

1. **Empty/failed spectrum**: `peak_count=0` or all peaks have `success=False`
   - Flag spectrum as outlier, exclude from per-peak stats

2. **Dummy peaks**: `was_detected=False` or `detection_quality="Missed"`
   - List separately, don't include in aggregate stats

3. **NaN/None values**: Missing volumes, linewidths
   - Skip in CV calculation, note in flags

4. **fix_linewidths=True**: LW should be constant
   - LW CV threshold drops to 0.05 (expect ~0)

5. **Cascade mode position drift**: Positions may drift systematically
   - Track trend, not just CV

6. **Different peak lists in comparison**: Series A has 45 peaks, Series B has 43
   - Report intersection (43 common) + exclusives (2 A-only, 0 B-only)

---

## Implementation Order

### Phase 1: Core Analyzer (TDD)
1. Write tests for `SeriesQCAnalyzer.load_training_data()`
2. Implement loading and parsing
3. Write tests for `compute_per_peak_stats()` (CV calculation)
4. Implement per-peak statistics
5. Write tests for `analyze_cluster_stability()`
6. Implement cluster tracking
7. Write tests for `generate_report()`
8. Implement report generation and save

### Phase 2: Comparator (TDD)
1. Write tests for `SeriesComparator.add_series()`
2. Implement series loading
3. Write tests for `get_common_peaks()` / `get_series_only_peaks()`
4. Implement intersection logic
5. Write tests for `compare()`
6. Implement comparison with winner selection

### Phase 3: GUI Dialog
1. Create dialog skeleton with layout
2. Implement drag-drop zone (copy from DynamiXs)
3. Implement summary dashboard
4. Implement flagged peaks table
5. Implement comparison table
6. Add Browse button for file selection
7. Connect to analyzer

### Phase 4: Menu Integration
1. Add menu entry in main_window.py
2. Add dialog export in dialogs/__init__.py
3. Test full workflow

---

## Success Criteria

1. Can load any valid `training_data.json` and generate report
2. Summary dashboard shows all key metrics with visual indicators
3. Flagged peaks table is sortable and shows relevant details
4. Comparison mode correctly identifies winner per metric
5. Handles edge cases (empty spectra, dummy peaks, NaN) gracefully
6. Auto-saves `QC_report.json` after analysis
7. Accessible from Modules menu, follows DynamiXs dialog pattern

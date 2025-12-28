#!/usr/bin/env python3
"""
Example Usage of LunaNMR Qt Design System

This script demonstrates how to use the design system and stylesheet
in a Qt/PySide6 application.

Run this script to see a demo window with styled widgets:
    python3 example_usage.py

Author: Guillaume Mas
Date: 2025
"""

import sys

try:
    from PySide6.QtWidgets import (
        QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
        QPushButton, QLabel, QLineEdit, QGroupBox, QCheckBox, QRadioButton,
        QProgressBar, QSlider, QComboBox, QSpinBox
    )
    from PySide6.QtCore import Qt
except ImportError:
    print("PySide6 not installed. Install with: pip install PySide6")
    sys.exit(1)

# Import design system constants and stylesheet loader
from lunaNMR.gui.styles import (
    # Core functions
    load_stylesheet,
    
    # Colors for dynamic styling
    SPACING_SM,
    SPACING_MD,
    SPACING_LG,
)


class DesignSystemDemo(QMainWindow):
    """Demo window showing all styled widgets from the design system."""
    
    def __init__(self):
        super().__init__()
        self.setWindowTitle("LunaNMR Design System Demo")
        self.setMinimumSize(800, 600)
        
        # Create central widget and main layout
        central_widget = QWidget()
        self.setCentralWidget(central_widget)
        
        main_layout = QVBoxLayout(central_widget)
        main_layout.setSpacing(SPACING_MD)
        main_layout.setContentsMargins(SPACING_LG, SPACING_LG, SPACING_LG, SPACING_LG)
        
        # Add demo sections
        self._add_header(main_layout)
        self._add_buttons_section(main_layout)
        self._add_inputs_section(main_layout)
        self._add_controls_section(main_layout)
        
    def _add_header(self, layout):
        """Add header with title and description."""
        title = QLabel("LunaNMR Design System Demo")
        title.setProperty("class", "header-large")
        layout.addWidget(title)
        
        description = QLabel("Demonstration of Qt/PySide6 styled widgets using main.qss")
        description.setProperty("class", "secondary")
        layout.addWidget(description)
        
    def _add_buttons_section(self, layout):
        """Add button examples."""
        group = QGroupBox("Buttons")
        group_layout = QHBoxLayout(group)
        group_layout.setSpacing(SPACING_SM)
        
        # Primary button
        primary_btn = QPushButton("Primary Action")
        primary_btn.setProperty("class", "primary")
        group_layout.addWidget(primary_btn)
        
        # Secondary button (default)
        secondary_btn = QPushButton("Secondary")
        group_layout.addWidget(secondary_btn)
        
        # Success button
        success_btn = QPushButton("Success")
        success_btn.setProperty("class", "success")
        group_layout.addWidget(success_btn)
        
        # Destructive button
        destructive_btn = QPushButton("Destructive")
        destructive_btn.setProperty("class", "destructive")
        group_layout.addWidget(destructive_btn)
        
        # Icon button
        icon_btn = QPushButton("🔄")
        icon_btn.setProperty("class", "icon")
        group_layout.addWidget(icon_btn)
        
        group_layout.addStretch()
        layout.addWidget(group)
        
    def _add_inputs_section(self, layout):
        """Add input widget examples."""
        group = QGroupBox("Input Fields")
        group_layout = QVBoxLayout(group)
        group_layout.setSpacing(SPACING_SM)
        
        # Text input
        line_edit = QLineEdit()
        line_edit.setPlaceholderText("Enter text here...")
        group_layout.addWidget(QLabel("Text Input:"))
        group_layout.addWidget(line_edit)
        
        # Combo box
        combo = QComboBox()
        combo.addItems(["Option 1", "Option 2", "Option 3"])
        group_layout.addWidget(QLabel("Dropdown:"))
        group_layout.addWidget(combo)
        
        # Spin box
        spin = QSpinBox()
        spin.setRange(0, 100)
        spin.setValue(50)
        group_layout.addWidget(QLabel("Numeric Input:"))
        group_layout.addWidget(spin)
        
        layout.addWidget(group)
        
    def _add_controls_section(self, layout):
        """Add control widget examples."""
        group = QGroupBox("Controls")
        group_layout = QVBoxLayout(group)
        group_layout.setSpacing(SPACING_SM)
        
        # Checkboxes
        checkbox1 = QCheckBox("Enable Option 1")
        checkbox1.setChecked(True)
        checkbox2 = QCheckBox("Enable Option 2")
        group_layout.addWidget(checkbox1)
        group_layout.addWidget(checkbox2)
        
        # Radio buttons
        radio_layout = QHBoxLayout()
        radio1 = QRadioButton("Choice A")
        radio1.setChecked(True)
        radio2 = QRadioButton("Choice B")
        radio_layout.addWidget(radio1)
        radio_layout.addWidget(radio2)
        group_layout.addLayout(radio_layout)
        
        # Slider
        slider = QSlider(Qt.Horizontal)
        slider.setRange(0, 100)
        slider.setValue(75)
        group_layout.addWidget(QLabel("Slider:"))
        group_layout.addWidget(slider)
        
        # Progress bar
        progress = QProgressBar()
        progress.setValue(65)
        group_layout.addWidget(QLabel("Progress:"))
        group_layout.addWidget(progress)
        
        layout.addWidget(group)
        layout.addStretch()


def main():
    """Run the demo application."""
    app = QApplication(sys.argv)
    
    # Load and apply the stylesheet
    try:
        stylesheet = load_stylesheet()
        app.setStyleSheet(stylesheet)
        print("✓ Stylesheet loaded successfully")
    except Exception as e:
        print(f"⚠ Warning: Could not load stylesheet: {e}")
    
    # Create and show demo window
    window = DesignSystemDemo()
    window.show()
    
    print("\nDesign System Demo")
    print("=" * 50)
    print("This window demonstrates the LunaNMR Qt design system")
    print("All widgets are styled using main.qss")
    print("\nClose the window to exit.")
    
    sys.exit(app.exec())


if __name__ == "__main__":
    main()

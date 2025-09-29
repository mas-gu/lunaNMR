"""
Simplified Parameter GUI Component
==================================

This module provides a clean GUI interface for the simplified 3-5 parameter system
that implements Priority 1 improvements by reducing parameter complexity.

Author: Guillaume Mas
Date: 2025
"""

import tkinter as tk
from tkinter import ttk, messagebox
import threading
from typing import Optional, Callable

class SimplifiedParameterFrame:
    """
    GUI component for simplified parameter control.

    Provides clean interface for 3-5 core parameters instead of 25+ legacy parameters.
    """

    def __init__(self, parent_frame, parameter_manager, on_parameter_change: Optional[Callable] = None):
        """
        Initialize simplified parameter frame.

        Parameters:
        -----------
        parent_frame : tk.Widget
            Parent GUI frame
        parameter_manager : NMRParameterManager
            Parameter manager instance
        on_parameter_change : callable, optional
            Callback function called when parameters change
        """
        self.parent_frame = parent_frame
        self.parameter_manager = parameter_manager
        self.on_parameter_change = on_parameter_change

        # Check if simplified mode is available
        self.simplified_available = (
            hasattr(parameter_manager, 'simplified_manager') and
            parameter_manager.simplified_manager is not None
        )

        # Create GUI variables
        self.create_variables()

        # Create the GUI
        self.create_gui()

        # Initialize values
        self.update_from_manager()

    def create_variables(self):
        """Create tkinter variables for simplified parameters"""

        # Mode selection
        self.use_simplified_mode = tk.BooleanVar(value=False)

        if self.simplified_available:
            # Core simplified parameters
            self.sensitivity = tk.DoubleVar(value=0.5)
            self.window_scale = tk.DoubleVar(value=1.0)
            self.quality_target = tk.DoubleVar(value=0.85)

            # Advanced parameters
            self.noise_method = tk.StringVar(value='auto')
            self.baseline_method = tk.StringVar(value='auto')

            # Status variables
            self.adaptive_threshold = tk.DoubleVar(value=0.85)
            self.parameter_count = tk.StringVar(value="Core: 3, Total: 5")

    def create_gui(self):
        """Create the simplified parameter GUI"""

        # Main frame
        self.main_frame = ttk.LabelFrame(self.parent_frame, text="🚀 Automated Fitting Parameters",
                                        padding="10")
        self.main_frame.pack(fill="x", padx=5, pady=5)

        if not self.simplified_available:
            # Show unavailable message
            unavailable_label = ttk.Label(
                self.main_frame,
                text="⚠️ Simplified parameters not available\nUsing legacy parameter system",
                foreground="orange"
            )
            unavailable_label.pack(pady=10)
            return

        # Mode selection frame
        mode_frame = ttk.Frame(self.main_frame)
        mode_frame.pack(fill="x", pady=(0, 10))

        mode_check = ttk.Checkbutton(
            mode_frame,
            text="🎯 Enable Simplified Mode (3-5 parameters instead of 25+)",
            variable=self.use_simplified_mode,
            command=self.on_mode_change
        )
        mode_check.pack(side="left")

        # Parameter count display
        count_label = ttk.Label(mode_frame, textvariable=self.parameter_count,
                               foreground="blue", font=("TkDefaultFont", 9))
        count_label.pack(side="right", padx=(10, 0))

        # Core parameters frame (only visible in simplified mode)
        self.core_frame = ttk.LabelFrame(self.main_frame, text="🎛️ Core Parameters",
                                        padding="5")

        # Sensitivity
        sens_frame = ttk.Frame(self.core_frame)
        sens_frame.pack(fill="x", pady=2)

        ttk.Label(sens_frame, text="Detection Sensitivity:", width=20).pack(side="left")
        sens_scale = ttk.Scale(sens_frame, from_=0.1, to=0.9, variable=self.sensitivity,
                              orient="horizontal", length=200, command=self.on_parameter_change_wrapper)
        sens_scale.pack(side="left", padx=(5, 5))

        sens_value = ttk.Label(sens_frame, width=8)
        sens_value.pack(side="left")

        # Update sensitivity label
        def update_sens_label(*args):
            sens_value.config(text=f"{self.sensitivity.get():.2f}")
        self.sensitivity.trace_add("write", update_sens_label)
        update_sens_label()

        # Window Scale
        window_frame = ttk.Frame(self.core_frame)
        window_frame.pack(fill="x", pady=2)

        ttk.Label(window_frame, text="Window Scale:", width=20).pack(side="left")
        window_scale = ttk.Scale(window_frame, from_=0.5, to=3.0, variable=self.window_scale,
                                orient="horizontal", length=200, command=self.on_parameter_change_wrapper)
        window_scale.pack(side="left", padx=(5, 5))

        window_value = ttk.Label(window_frame, width=8)
        window_value.pack(side="left")

        def update_window_label(*args):
            window_value.config(text=f"{self.window_scale.get():.2f}")
        self.window_scale.trace_add("write", update_window_label)
        update_window_label()

        # Quality Target
        quality_frame = ttk.Frame(self.core_frame)
        quality_frame.pack(fill="x", pady=2)

        ttk.Label(quality_frame, text="Quality Target:", width=20).pack(side="left")
        quality_scale = ttk.Scale(quality_frame, from_=0.5, to=0.95, variable=self.quality_target,
                                 orient="horizontal", length=200, command=self.on_parameter_change_wrapper)
        quality_scale.pack(side="left", padx=(5, 5))

        quality_value = ttk.Label(quality_frame, width=8)
        quality_value.pack(side="left")

        def update_quality_label(*args):
            quality_value.config(text=f"{self.quality_target.get():.2f}")
        self.quality_target.trace_add("write", update_quality_label)
        update_quality_label()

        # Advanced parameters frame (collapsible)
        self.advanced_frame = ttk.LabelFrame(self.main_frame, text="🔧 Advanced Options",
                                           padding="5")

        # Noise estimation method
        noise_frame = ttk.Frame(self.advanced_frame)
        noise_frame.pack(fill="x", pady=2)

        ttk.Label(noise_frame, text="Noise Estimation:", width=20).pack(side="left")
        noise_combo = ttk.Combobox(noise_frame, textvariable=self.noise_method,
                                  values=['auto', 'robust', 'percentile'], state="readonly", width=15)
        noise_combo.pack(side="left", padx=(5, 0))
        noise_combo.bind('<<ComboboxSelected>>', self.on_parameter_change_wrapper)

        # Baseline method
        baseline_frame = ttk.Frame(self.advanced_frame)
        baseline_frame.pack(fill="x", pady=2)

        ttk.Label(baseline_frame, text="Baseline Method:", width=20).pack(side="left")
        baseline_combo = ttk.Combobox(baseline_frame, textvariable=self.baseline_method,
                                     values=['auto', 'polynomial', 'iterative'], state="readonly", width=15)
        baseline_combo.pack(side="left", padx=(5, 0))
        baseline_combo.bind('<<ComboboxSelected>>', self.on_parameter_change_wrapper)

        # Status frame
        self.status_frame = ttk.LabelFrame(self.main_frame, text="📊 Adaptive Status",
                                          padding="5")

        # Adaptive threshold display
        threshold_frame = ttk.Frame(self.status_frame)
        threshold_frame.pack(fill="x", pady=2)

        ttk.Label(threshold_frame, text="Adaptive Threshold:", width=20).pack(side="left")
        threshold_label = ttk.Label(threshold_frame, textvariable=self.adaptive_threshold,
                                   foreground="darkgreen", width=8)
        threshold_label.pack(side="left", padx=(5, 0))

        ttk.Label(threshold_frame, text="(SNR-based)", foreground="gray").pack(side="left", padx=(5, 0))

        # Control buttons
        button_frame = ttk.Frame(self.main_frame)
        button_frame.pack(fill="x", pady=(10, 0))

        ttk.Button(button_frame, text="Reset to Defaults", command=self.reset_to_defaults).pack(side="left")
        ttk.Button(button_frame, text="Apply Parameters", command=self.apply_parameters).pack(side="left", padx=(5, 0))

        # Help button
        ttk.Button(button_frame, text="Help", command=self.show_help).pack(side="right")

        # Initially hide frames
        self.update_frame_visibility()

    def on_mode_change(self):
        """Handle simplified mode toggle"""

        if self.use_simplified_mode.get():
            # Enable simplified mode
            self.parameter_manager.enable_simplified_mode(
                sensitivity=self.sensitivity.get(),
                window_scale=self.window_scale.get(),
                quality_target=self.quality_target.get(),
                noise_estimation_method=self.noise_method.get(),
                baseline_method=self.baseline_method.get()
            )
        else:
            # Disable simplified mode
            self.parameter_manager.disable_simplified_mode()

        self.update_frame_visibility()
        self.update_parameter_count()

        # Notify parent of parameter change
        if self.on_parameter_change:
            self.on_parameter_change()

    def update_frame_visibility(self):
        """Update visibility of parameter frames based on mode"""

        if not self.simplified_available:
            return

        if self.use_simplified_mode.get():
            self.core_frame.pack(fill="x", pady=5)
            self.advanced_frame.pack(fill="x", pady=5)
            self.status_frame.pack(fill="x", pady=5)
        else:
            self.core_frame.pack_forget()
            self.advanced_frame.pack_forget()
            self.status_frame.pack_forget()

    def on_parameter_change_wrapper(self, *args):
        """Wrapper for parameter change events"""

        if not self.use_simplified_mode.get():
            return

        # Update parameter manager
        self.parameter_manager.update_simplified_parameters(
            sensitivity=self.sensitivity.get(),
            window_scale=self.window_scale.get(),
            quality_target=self.quality_target.get(),
            noise_estimation_method=self.noise_method.get(),
            baseline_method=self.baseline_method.get()
        )

        # Update adaptive threshold (placeholder - would need spectrum data)
        # For now, just show the quality target
        self.adaptive_threshold.set(self.quality_target.get())

        # Notify parent
        if self.on_parameter_change:
            threading.Thread(target=self.on_parameter_change, daemon=True).start()

    def update_parameter_count(self):
        """Update parameter count display"""

        if self.use_simplified_mode.get():
            self.parameter_count.set("Core: 3, Total: 5")
        else:
            self.parameter_count.set("Legacy: 25+ parameters")

    def update_from_manager(self):
        """Update GUI from parameter manager state"""

        if not self.simplified_available:
            return

        # Update mode
        self.use_simplified_mode.set(self.parameter_manager.use_simplified_mode)

        if self.parameter_manager.use_simplified_mode and self.parameter_manager.simplified_manager:
            # Update values from simplified manager
            params = self.parameter_manager.simplified_manager.simplified_params

            self.sensitivity.set(params.sensitivity)
            self.window_scale.set(params.window_scale)
            self.quality_target.set(params.quality_target)
            self.noise_method.set(params.noise_estimation_method)
            self.baseline_method.set(params.baseline_method)

        self.update_frame_visibility()
        self.update_parameter_count()

    def reset_to_defaults(self):
        """Reset parameters to defaults"""

        if not self.simplified_available:
            return

        if self.use_simplified_mode.get():
            self.parameter_manager.simplified_manager.reset_to_defaults()
            self.update_from_manager()
        else:
            self.parameter_manager.reset_to_defaults()

        messagebox.showinfo("Reset Complete", "Parameters reset to defaults")

    def apply_parameters(self):
        """Apply current parameters"""

        if not self.simplified_available:
            messagebox.showinfo("Applied", "Using legacy parameters")
            return

        if self.use_simplified_mode.get():
            # Validate simplified parameters
            is_valid, errors = self.parameter_manager.validate_simplified_parameters()

            if not is_valid:
                messagebox.showerror("Validation Error", "Parameter validation failed:\n" + "\n".join(errors))
                return

            messagebox.showinfo("Applied", "Simplified parameters applied successfully")
        else:
            # Validate legacy parameters
            errors = self.parameter_manager.validate_all_parameters()

            if errors:
                messagebox.showerror("Validation Error", "Parameter validation failed:\n" + "\n".join(errors))
                return

            messagebox.showinfo("Applied", "Legacy parameters applied successfully")

        # Notify parent
        if self.on_parameter_change:
            self.on_parameter_change()

    def show_help(self):
        """Show help dialog for simplified parameters"""

        help_text = """
🚀 Automated Fitting Parameters Help

SIMPLIFIED MODE (Recommended):
Uses only 3-5 core parameters instead of 25+ legacy parameters.

CORE PARAMETERS:
• Detection Sensitivity (0.1-0.9): Higher = more peaks detected
• Window Scale (0.5-3.0): Adjusts fitting window sizes
• Quality Target (0.5-0.95): Target R² value for fits

ADVANCED OPTIONS:
• Noise Estimation: Method for noise level calculation
• Baseline Method: Baseline correction approach

BENEFITS:
✅ Dramatically simplified parameter tuning
✅ Adaptive quality thresholds based on SNR
✅ Natural parameter space prevents bounds violations
✅ Automatic parameter derivation for all legacy systems
✅ Full backward compatibility maintained

LEGACY MODE:
Uses traditional 25+ parameter system for maximum control.

The simplified system automatically calculates optimal values for all
legacy parameters based on your 3-5 core settings and spectrum characteristics.
"""

        help_window = tk.Toplevel(self.parent_frame)
        help_window.title("Simplified Parameters Help")
        help_window.geometry("600x500")

        text_widget = tk.Text(help_window, wrap="word", padx=10, pady=10)
        text_widget.pack(fill="both", expand=True)
        text_widget.insert("1.0", help_text)
        text_widget.config(state="disabled")

        # Add scrollbar
        scrollbar = ttk.Scrollbar(help_window, orient="vertical", command=text_widget.yview)
        scrollbar.pack(side="right", fill="y")
        text_widget.config(yscrollcommand=scrollbar.set)

    def get_current_mode(self):
        """Get current parameter mode"""
        return "simplified" if self.use_simplified_mode.get() else "legacy"

    def get_simplified_parameters(self):
        """Get current simplified parameter values"""
        if not self.simplified_available or not self.use_simplified_mode.get():
            return None

        return {
            'sensitivity': self.sensitivity.get(),
            'window_scale': self.window_scale.get(),
            'quality_target': self.quality_target.get(),
            'noise_estimation_method': self.noise_method.get(),
            'baseline_method': self.baseline_method.get()
        }

    def set_adaptive_threshold(self, threshold: float):
        """Update adaptive threshold display (called from spectrum analysis)"""
        if self.simplified_available:
            self.adaptive_threshold.set(threshold)
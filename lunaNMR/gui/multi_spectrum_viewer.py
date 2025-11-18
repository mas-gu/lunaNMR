# ABOUTME: Multi-spectrum overlay viewer for simultaneous visualization of multiple NMR spectra
# ABOUTME: Allows toggling spectrum visibility, color customization, and per-spectrum Voigt analysis

import tkinter as tk
from tkinter import ttk, colorchooser, messagebox
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure
import numpy as np
import pandas as pd
import os
from pathlib import Path

# Import LunaNMR components
from lunaNMR.gui.gui_components import PeakNavigator, natural_sort_key
from lunaNMR.gui.visualization import VoigtAnalysisPlotter


class SpectrumData:
    """Data container for individual spectrum in multi-spectrum viewer"""

    def __init__(self, name, file_path, result_data):
        self.name = name
        self.file_path = file_path
        self.result_data = result_data
        self.color = None  # Assigned display color
        self.visible = False  # Visibility toggle state
        self.data = None  # NMR data (loaded on demand)
        self.integrator = None  # EnhancedVoigtIntegrator instance
        self.fitted_peaks = None  # Fitted peak results
        self.detected_peaks = None  # Detected peak list

    def load_data(self, file_manager):
        """Load spectrum data on demand"""
        if self.data is None:
            print(f"📂 Loading data for {self.name}...")
            try:
                # Check file exists
                if not os.path.exists(self.file_path):
                    print(f"❌ File not found: {self.file_path}")
                    return False

                # Create integrator (takes no parameters)
                from lunaNMR.core.core_integrator import EnhancedVoigtIntegrator
                self.integrator = EnhancedVoigtIntegrator()

                # Load NMR file using integrator's load method
                success = self.integrator.load_nmr_file(self.file_path)
                if not success:
                    print(f"❌ Failed to load NMR file: {self.file_path}")
                    return False

                # Verify data loaded
                if not hasattr(self.integrator, 'nmr_data') or self.integrator.nmr_data is None:
                    print(f"❌ NMR data is None after loading")
                    return False

                # Verify PPM axes are available
                if not hasattr(self.integrator, 'ppm_x_axis') or not hasattr(self.integrator, 'ppm_y_axis'):
                    print(f"❌ Missing PPM axes after loading")
                    return False

                self.data = self.integrator.nmr_data

                # Load fitted peaks from result_data
                # multi_spectrum_processor stores results in 'integration_results'
                self.fitted_peaks = self.result_data.get('integration_results', [])
                if not self.fitted_peaks:
                    # Fallback to 'fitted_peaks' for compatibility with legacy formats
                    self.fitted_peaks = self.result_data.get('fitted_peaks', [])

                # Ensure fitted_peaks is a list, not None
                if self.fitted_peaks is None:
                    self.fitted_peaks = []

                self.detected_peaks = self.result_data.get('detected_peaks', [])
                print(f"✅ Loaded {self.name}: {len(self.fitted_peaks)} fitted peaks")
            except Exception as e:
                print(f"❌ Error loading {self.name}: {e}")
                import traceback
                traceback.print_exc()
                return False
        return self.data is not None


class MultiSpectrumViewer:
    """Multi-spectrum overlay viewer for simultaneous visualization"""

    def __init__(self, parent, file_manager, all_results):
        self.parent = parent
        self.file_manager = file_manager
        self.all_results = all_results

        # Check for empty results
        if not all_results or len(all_results) == 0:
            messagebox.showwarning("No Data", "No spectra found in series results")
            raise ValueError("No spectra in all_results")

        # Spectrum data management
        self.spectra = []  # List of SpectrumData objects
        self.reference_peaks = None  # Peak list from first spectrum

        # Create viewer window with dynamic screen scaling
        self.window = tk.Toplevel(parent)
        self.window.title("Multi-Spectrum Overlay Viewer")

        # Dynamic sizing: Adaptive to screen size
        screen_width = parent.winfo_screenwidth()
        screen_height = parent.winfo_screenheight()

        # Use 85% of screen dimensions (slightly larger than spectrum_browser's 80%×75%)
        # But cap at reasonable maximums for very large screens
        window_width = min(int(screen_width * 0.85), 1600)
        window_height = min(int(screen_height * 0.85), 1000)

        # Set minimum size based on screen size (more adaptive for small screens)
        # For screens <1400px wide, reduce minimum to 900px to ensure buttons visible
        min_width = 900 if screen_width < 1400 else 1000
        min_height = 600 if screen_height < 900 else 700

        self.window.geometry(f"{window_width}x{window_height}")
        self.window.minsize(min_width, min_height)

        print(f"📐 Multi-spectrum viewer window size: {window_width}×{window_height} (screen: {screen_width}×{screen_height})")
        print(f"   Minimum size: {min_width}×{min_height}")

        # Register cleanup handler
        self.window.protocol("WM_DELETE_WINDOW", self.on_closing)

        # Contour control variables (with master window for proper cleanup)
        # Use low detail preset as default (same as spectrum_browser.py low preset)
        self.contour_levels = tk.IntVar(master=self.window, value=10)
        self.contour_min = tk.DoubleVar(master=self.window, value=0.2)
        self.contour_increment = tk.DoubleVar(master=self.window, value=2.0)

        # Voigt analysis state
        self.selected_spectrum_index = 0  # Currently selected spectrum for Voigt analysis
        self.selected_peak_index = None  # Currently selected peak for Voigt analysis

        # Color palette for spectra (matplotlib default cycle)
        self.default_colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd',
                               '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']

        # Initialize spectrum data
        self.initialize_spectra()

        # Setup UI
        self.setup_viewer()

        # Load first spectrum by default
        self.window.after(100, self.initialize_first_spectrum)

    def initialize_spectra(self):
        """Initialize SpectrumData objects for all spectra in series"""
        print(f"🔧 Initializing multi-spectrum viewer with {len(self.all_results)} spectra")

        # Create SpectrumData objects
        temp_spectra = []
        for idx, result in enumerate(self.all_results):
            spectrum_name = result.get('spectrum_name', f'spectrum_{idx+1:03d}')
            file_path = result.get('spectrum_file', '')

            # Create SpectrumData object
            spec_data = SpectrumData(spectrum_name, file_path, result)
            temp_spectra.append(spec_data)

        # Sort spectra using natural sort by name
        # Then REVERSE so reference (first in time) is at bottom, newest at top
        temp_spectra.sort(key=lambda s: natural_sort_key(s.name), reverse=True)

        print(f"📊 Spectrum order (newest → oldest):")
        for i, spec in enumerate(temp_spectra):
            print(f"   {i+1}. {spec.name}")

        # Assign colors and visibility after sorting
        for idx, spec_data in enumerate(temp_spectra):
            # Assign color from palette (cycle through colors)
            spec_data.color = self.default_colors[idx % len(self.default_colors)]

            # Last spectrum (reference, oldest) visible by default
            spec_data.visible = (idx == len(temp_spectra) - 1)

            self.spectra.append(spec_data)

        print(f"✅ Initialized {len(self.spectra)} spectra (reference at bottom: {self.spectra[-1].name})")

    def initialize_first_spectrum(self):
        """Load reference spectrum (last in list) and populate reference peak list"""
        if len(self.spectra) > 0:
            # Reference spectrum is at the bottom (last in list after reverse sort)
            reference_spectrum = self.spectra[-1]

            # Load reference spectrum data
            if reference_spectrum.load_data(self.file_manager):
                # Extract reference peaks from reference spectrum
                self.reference_peaks = reference_spectrum.fitted_peaks

                # Check if peaks loaded
                if not self.reference_peaks:
                    print(f"⚠️ No fitted peaks found in {reference_spectrum.name}")
                    self.reference_peaks = []  # Ensure it's a list, not None
                    messagebox.showwarning("No Peaks",
                        f"No fitted peaks found in {reference_spectrum.name}.\n"
                        "Peak list will be empty.")
                else:
                    print(f"📋 Reference peak list: {len(self.reference_peaks)} peaks from {reference_spectrum.name}")

                # Populate peak list in main tab
                self.populate_main_peak_list()

                # Plot reference spectrum
                self.update_overlay_plot()
            else:
                messagebox.showerror("Load Error",
                    f"Failed to load reference spectrum: {reference_spectrum.name}")

    def on_closing(self):
        """Cleanup when window closes"""
        try:
            print("🧹 Cleaning up multi-spectrum viewer...")
            # Close matplotlib figures
            plt.close(self.fig_overlay)
            plt.close(self.fig_voigt)
            print("✅ Matplotlib figures closed")
        except Exception as e:
            print(f"⚠️ Error during cleanup: {e}")
        finally:
            # Destroy window
            self.window.destroy()

    def setup_viewer(self):
        """Setup the multi-spectrum viewer layout with tabbed interface"""
        # Create notebook with tabs
        self.notebook = ttk.Notebook(self.window)
        self.notebook.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)

        # Tab 1: Main Spectrum Overlay
        self.main_tab = ttk.Frame(self.notebook)
        self.notebook.add(self.main_tab, text="📊 Spectrum Overlay")
        self.setup_main_overlay_tab()

        # Tab 2: Voigt Analysis
        self.voigt_tab = ttk.Frame(self.notebook)
        self.notebook.add(self.voigt_tab, text="📈 Voigt Analysis")
        self.setup_voigt_analysis_tab()

        # Tab 3: 3D Voigt Analysis (supplementary visualization)
        self.voigt_3d_tab = ttk.Frame(self.notebook)
        self.notebook.add(self.voigt_3d_tab, text="🎨 3D Voigt Analysis")
        self.setup_voigt_3d_analysis_tab()

    def setup_main_overlay_tab(self):
        """Setup main spectrum overlay tab with grid-based layout for guaranteed visibility"""
        # Configure grid weights for main_tab (3 rows: buttons, content, controls)
        # Row 0: Spectrum buttons - fixed height (weight=0)
        # Row 1: Plot + peak list - expandable (weight=1)
        # Row 2: Contour controls - fixed height (weight=0)
        self.main_tab.grid_rowconfigure(0, weight=0, minsize=80)   # Fixed ~80px for buttons
        self.main_tab.grid_rowconfigure(1, weight=1)                # Expandable plot area
        self.main_tab.grid_rowconfigure(2, weight=0, minsize=100)  # Fixed ~100px for controls
        self.main_tab.grid_columnconfigure(0, weight=1)

        print("📐 Using grid-based layout: Row 0 (buttons, fixed 80px), Row 1 (plot, expandable), Row 2 (controls, fixed 100px)")

        # Row 0: Spectrum control bar (fixed height)
        control_bar = ttk.Frame(self.main_tab)
        control_bar.grid(row=0, column=0, sticky='ew', padx=5, pady=5)

        ttk.Label(control_bar, text="Spectra:", font=('TkDefaultFont', 10, 'bold')).pack(side=tk.LEFT, padx=5)

        # Scrollable spectrum button container with adaptive height
        # Height adjusted to accommodate buttons with checkboxes (50px ensures all content visible)
        self.button_canvas = tk.Canvas(control_bar, height=50, highlightthickness=0)
        button_scrollbar = ttk.Scrollbar(control_bar, orient="horizontal", command=self.button_canvas.xview)
        button_frame = ttk.Frame(self.button_canvas)

        button_frame.bind(
            "<Configure>",
            lambda e: self.button_canvas.configure(scrollregion=self.button_canvas.bbox("all"))
        )

        self.button_canvas.create_window((0, 0), window=button_frame, anchor="nw")
        self.button_canvas.configure(xscrollcommand=button_scrollbar.set)

        # Enable horizontal mouse wheel scrolling for button container
        def scroll_buttons(event):
            self.button_canvas.xview_scroll(int(-1 * (event.delta / 120)), "units")

        # Bind mouse wheel to canvas and all child widgets
        self.button_canvas.bind("<MouseWheel>", scroll_buttons)  # Windows/Linux
        self.button_canvas.bind("<Button-4>", lambda e: self.button_canvas.xview_scroll(-1, "units"))  # Linux scroll up
        self.button_canvas.bind("<Button-5>", lambda e: self.button_canvas.xview_scroll(1, "units"))   # Linux scroll down

        self.button_canvas.pack(side=tk.TOP, fill=tk.X, expand=True)
        button_scrollbar.pack(side=tk.TOP, fill=tk.X)

        # Create spectrum toggle buttons
        self.spectrum_buttons = []
        for idx, spec_data in enumerate(self.spectra):
            btn_frame = self.create_spectrum_button(button_frame, idx, spec_data)
            btn_frame.pack(side=tk.LEFT, padx=2)

        # Row 1: Horizontal paned window (spectrum plot left, peak list right) - expandable
        content_paned = ttk.PanedWindow(self.main_tab, orient='horizontal')
        content_paned.grid(row=1, column=0, sticky='nsew', padx=5, pady=(5, 0))

        # Left panel (75%): Spectrum overlay plot
        left_panel = ttk.Frame(content_paned)
        content_paned.add(left_panel, weight=3)

        plot_container = ttk.LabelFrame(left_panel, text="📈 Spectrum Overlay", padding=5)
        plot_container.pack(fill=tk.BOTH, expand=True)

        # Create matplotlib figure for overlay plot with compact size
        # Reduced from (10, 8) to (6, 4) to ensure bottom controls always visible
        # Figure will resize with window but starts at reasonable size
        self.fig_overlay = Figure(figsize=(6, 4))
        self.ax_overlay = self.fig_overlay.add_subplot(111)

        self.canvas_overlay = FigureCanvasTkAgg(self.fig_overlay, plot_container)
        self.canvas_overlay.get_tk_widget().pack(fill=tk.BOTH, expand=True, pady=(0, 5))

        # Add toolbar
        toolbar = NavigationToolbar2Tk(self.canvas_overlay, plot_container)
        toolbar.update()

        # Right panel (25%): Peak list
        right_panel = ttk.Frame(content_paned)
        content_paned.add(right_panel, weight=1)

        peak_list_container = ttk.LabelFrame(right_panel, text="📋 Peak List (Reference)", padding=5)
        peak_list_container.pack(fill=tk.BOTH, expand=True)

        # Create peak list treeview
        self.setup_main_peak_list(peak_list_container)

        # Row 2: Contour controls (fixed height, always visible)
        contour_frame = ttk.LabelFrame(self.main_tab, text="🎨 Contour Controls (Applied to All Spectra)", padding=5)
        contour_frame.grid(row=2, column=0, sticky='ew', padx=5, pady=5)
        self.setup_contour_controls(contour_frame)

    def create_spectrum_button(self, parent, idx, spec_data):
        """Create spectrum toggle button with color box and checkbox"""
        btn_frame = ttk.Frame(parent, relief='raised', borderwidth=1)

        # Color box (clickable to change color)
        color_canvas = tk.Canvas(btn_frame, width=20, height=20, bg=spec_data.color,
                                highlightthickness=1, highlightbackground='gray')
        color_canvas.pack(side=tk.LEFT, padx=2)

        # Right-click to change color
        color_canvas.bind('<Button-3>', lambda e: self.change_spectrum_color(idx))
        color_canvas.bind('<Button-2>', lambda e: self.change_spectrum_color(idx))  # Mac middle click

        # Spectrum name label (2x larger to show full names)
        name_label = ttk.Label(btn_frame, text=spec_data.name, width=30)
        name_label.pack(side=tk.LEFT, padx=2)

        # Visibility checkbox (with master for proper cleanup)
        var = tk.BooleanVar(master=self.window, value=spec_data.visible)
        checkbox = ttk.Checkbutton(btn_frame, variable=var,
                                   command=lambda: self.toggle_spectrum_visibility(idx, var))
        checkbox.pack(side=tk.LEFT, padx=2)

        # Store references
        self.spectrum_buttons.append({
            'frame': btn_frame,
            'color_canvas': color_canvas,
            'name_label': name_label,
            'checkbox': checkbox,
            'var': var
        })

        return btn_frame

    def change_spectrum_color(self, idx):
        """Open color picker to change spectrum color"""
        spec_data = self.spectra[idx]
        color = colorchooser.askcolor(title=f"Choose color for {spec_data.name}",
                                      initialcolor=spec_data.color)
        # Check color is not None (user cancelled) and has hex value
        if color and color[1]:  # color[1] is hex string
            spec_data.color = color[1]
            # Update color box
            self.spectrum_buttons[idx]['color_canvas'].config(bg=spec_data.color)
            # Redraw plot if spectrum is visible
            if spec_data.visible:
                self.update_overlay_plot()

    def toggle_spectrum_visibility(self, idx, var):
        """Toggle spectrum visibility and update plot"""
        spec_data = self.spectra[idx]
        spec_data.visible = var.get()

        print(f"{'✅' if spec_data.visible else '❌'} {spec_data.name} visibility: {spec_data.visible}")

        # Load data if needed
        if spec_data.visible and spec_data.data is None:
            if not spec_data.load_data(self.file_manager):
                messagebox.showerror("Error", f"Failed to load {spec_data.name}")
                var.set(False)
                spec_data.visible = False
                return

        # Update overlay plot
        self.update_overlay_plot()

    def setup_main_peak_list(self, parent):
        """Setup peak list treeview in main tab"""
        # Peak list selector dropdown (which spectrum's peaks to display)
        selector_frame = ttk.Frame(parent)
        selector_frame.pack(fill=tk.X, pady=(0, 5))

        ttk.Label(selector_frame, text="Show peaks from:").pack(side=tk.LEFT, padx=5)
        self.peak_list_selector = ttk.Combobox(selector_frame, state='readonly', width=20)
        self.peak_list_selector.pack(side=tk.LEFT, padx=5)
        self.peak_list_selector.bind('<<ComboboxSelected>>', self.on_peak_list_changed)

        # Populate selector with spectrum names
        spectrum_names = [spec.name for spec in self.spectra]
        self.peak_list_selector['values'] = spectrum_names
        if spectrum_names:
            self.peak_list_selector.current(0)

        # Create treeview with scrollbars
        tree_frame = ttk.Frame(parent)
        tree_frame.pack(fill=tk.BOTH, expand=True)

        # Scrollbars
        vsb = ttk.Scrollbar(tree_frame, orient="vertical")
        hsb = ttk.Scrollbar(tree_frame, orient="horizontal")

        # Treeview
        self.main_peak_tree = ttk.Treeview(tree_frame,
                                          columns=('Assignment', 'X_ppm', 'Y_ppm', 'Height'),
                                          show='headings',
                                          yscrollcommand=vsb.set,
                                          xscrollcommand=hsb.set)

        vsb.config(command=self.main_peak_tree.yview)
        hsb.config(command=self.main_peak_tree.xview)

        # Configure columns
        self.main_peak_tree.heading('Assignment', text='Assignment')
        self.main_peak_tree.heading('X_ppm', text='X (ppm)')
        self.main_peak_tree.heading('Y_ppm', text='Y (ppm)')
        self.main_peak_tree.heading('Height', text='Height')

        self.main_peak_tree.column('Assignment', width=80, anchor='center')
        self.main_peak_tree.column('X_ppm', width=80, anchor='center')
        self.main_peak_tree.column('Y_ppm', width=80, anchor='center')
        self.main_peak_tree.column('Height', width=100, anchor='center')

        # Pack
        self.main_peak_tree.grid(row=0, column=0, sticky='nsew')
        vsb.grid(row=0, column=1, sticky='ns')
        hsb.grid(row=1, column=0, sticky='ew')

        tree_frame.grid_rowconfigure(0, weight=1)
        tree_frame.grid_columnconfigure(0, weight=1)

        # Double-click to center on peak
        self.main_peak_tree.bind('<Double-1>', self.on_main_peak_double_click)

    def on_peak_list_changed(self, event=None):
        """Handle peak list selector change"""
        selected_idx = self.peak_list_selector.current()
        if selected_idx < 0:
            return

        spec_data = self.spectra[selected_idx]
        print(f"📋 Changed peak list to: {spec_data.name}")

        # Load spectrum data if not loaded yet
        if spec_data.data is None:
            if not spec_data.load_data(self.file_manager):
                messagebox.showerror("Error", f"Failed to load {spec_data.name}")
                return

        # Update reference peaks to this spectrum's peaks
        self.reference_peaks = spec_data.fitted_peaks if spec_data.fitted_peaks else []

        # Refresh peak list display
        self.populate_main_peak_list()

        # Update overlay plot to show peaks from selected spectrum
        self.update_overlay_plot()

    def populate_main_peak_list(self):
        """Populate main peak list with reference peaks"""
        # Clear existing items
        for item in self.main_peak_tree.get_children():
            self.main_peak_tree.delete(item)

        if not self.reference_peaks:
            return

        # Debug: Print first peak structure to see available keys
        if len(self.reference_peaks) > 0:
            print(f"🔍 DEBUG: First peak keys: {list(self.reference_peaks[0].keys())}")
            print(f"🔍 DEBUG: First peak sample: {self.reference_peaks[0]}")

        # Populate with reference peaks
        for peak in self.reference_peaks:
            assignment = peak.get('assignment', 'Unknown')
            # Try multiple possible keys for coordinates
            x_ppm = peak.get('center_x', peak.get('peak_x', peak.get('ppm_x', peak.get('x', 0))))
            y_ppm = peak.get('center_y', peak.get('peak_y', peak.get('ppm_y', peak.get('y', 0))))
            height = peak.get('height', peak.get('amplitude', peak.get('intensity', 0)))

            self.main_peak_tree.insert('', 'end', values=(
                assignment,
                f"{x_ppm:.3f}",
                f"{y_ppm:.1f}",
                f"{height:.2e}" if height > 0 else "N/A"
            ))

        print(f"📋 Populated main peak list: {len(self.reference_peaks)} peaks")

    def on_main_peak_double_click(self, event):
        """Handle double-click on peak in main peak list"""
        selection = self.main_peak_tree.selection()
        if not selection:
            return

        # Get peak index
        item = selection[0]
        peak_idx = self.main_peak_tree.index(item)

        if peak_idx < len(self.reference_peaks):
            peak = self.reference_peaks[peak_idx]
            # Try multiple possible keys for coordinates
            x_ppm = peak.get('center_x', peak.get('peak_x', peak.get('ppm_x', peak.get('x', 0))))
            y_ppm = peak.get('center_y', peak.get('peak_y', peak.get('ppm_y', peak.get('y', 0))))

            # Center plot on peak
            self.center_overlay_on_peak(x_ppm, y_ppm)

    def center_overlay_on_peak(self, x_ppm, y_ppm, window=0.5):
        """Center overlay plot on specified peak position"""
        # Set axis limits centered on peak
        self.ax_overlay.set_xlim(x_ppm - window, x_ppm + window)
        self.ax_overlay.set_ylim(y_ppm - window*10, y_ppm + window*10)  # Wider Y window

        self.canvas_overlay.draw()
        print(f"🎯 Centered on peak at ({x_ppm:.3f}, {y_ppm:.1f})")

    def setup_contour_controls(self, parent):
        """Setup contour control widgets"""
        # Levels
        ttk.Label(parent, text="Levels:").grid(row=0, column=0, padx=5, sticky='w')
        levels_spin = ttk.Spinbox(parent, from_=5, to=100, width=8,
                                 textvariable=self.contour_levels)
        levels_spin.grid(row=0, column=1, padx=5, sticky='w')

        # Min level
        ttk.Label(parent, text="Min Level:").grid(row=0, column=2, padx=5, sticky='w')
        min_spin = ttk.Spinbox(parent, from_=0.01, to=10.0, increment=0.01, width=8,
                              textvariable=self.contour_min)
        min_spin.grid(row=0, column=3, padx=5, sticky='w')

        # Increment
        ttk.Label(parent, text="Increment:").grid(row=0, column=4, padx=5, sticky='w')
        inc_spin = ttk.Spinbox(parent, from_=1.01, to=2.0, increment=0.01, width=8,
                              textvariable=self.contour_increment)
        inc_spin.grid(row=0, column=5, padx=5, sticky='w')

        # Update button
        update_btn = ttk.Button(parent, text="Update Plot", command=self.update_overlay_plot)
        update_btn.grid(row=0, column=6, padx=20)

        # Export button
        export_btn = ttk.Button(parent, text="Export PNG", command=self.export_overlay_plot)
        export_btn.grid(row=0, column=7, padx=5)

    def update_overlay_plot(self):
        """Redraw overlay plot with all visible spectra"""
        # Store current axis limits to preserve zoom
        xlim_before = self.ax_overlay.get_xlim() if self.ax_overlay.get_xlim() != (0.0, 1.0) else None
        ylim_before = self.ax_overlay.get_ylim() if self.ax_overlay.get_ylim() != (0.0, 1.0) else None

        # Clear axes (resets inversion state)
        self.ax_overlay.clear()

        # Get visible spectra
        visible_spectra = [s for s in self.spectra if s.visible and s.data is not None]

        if not visible_spectra:
            self.ax_overlay.text(0.5, 0.5, 'No spectra selected\nCheck spectrum boxes to display',
                               ha='center', va='center', transform=self.ax_overlay.transAxes)
            self.canvas_overlay.draw()
            return

        print(f"🎨 Updating overlay plot with {len(visible_spectra)} spectra")

        # Get contour parameters
        num_levels = self.contour_levels.get()
        min_level = self.contour_min.get()
        increment = self.contour_increment.get()  # Not used in geomspace method, kept for UI compatibility

        # Plot spectra in reverse order (reference first = background, newest last = foreground)
        # This ensures n, n+1, n+2... plotting order with reference at bottom
        plotting_order = list(reversed(visible_spectra))

        # Plot each visible spectrum
        for spec_data in plotting_order:
            if spec_data.integrator is None:
                continue

            # Verify PPM axes are available
            if not hasattr(spec_data.integrator, 'ppm_x_axis') or not hasattr(spec_data.integrator, 'ppm_y_axis'):
                print(f"⚠️ {spec_data.name}: Missing PPM axes, skipping")
                continue

            # Get axes from integrator
            ppm_x = spec_data.integrator.ppm_x_axis
            ppm_y = spec_data.integrator.ppm_y_axis
            data = spec_data.data

            # Calculate contour levels using same method as spectrum_browser.py
            # Use geometrically spaced levels based on data intensity range
            data_abs = np.abs(data)
            max_intensity = np.max(data_abs)
            min_intensity = min_level * max_intensity

            # Use geomspace for better contour distribution (same as spectrum_browser.py)
            level_values = np.geomspace(min_intensity, max_intensity * 0.8, num_levels)

            # Plot contours with spectrum color
            self.ax_overlay.contour(ppm_x, ppm_y, data,
                                   levels=level_values,
                                   colors=spec_data.color,
                                   linewidths=1.0,
                                   alpha=0.8)

            print(f"   ✓ Plotted {spec_data.name} (color: {spec_data.color})")

        # Detect nucleus type for axis labels
        y_label = '¹⁵N/¹³C (ppm)'  # Generic label
        if visible_spectra and visible_spectra[0].integrator:
            nucleus = getattr(visible_spectra[0].integrator, 'detected_nucleus_type', None)
            if nucleus == '15N-HSQC':
                y_label = '¹⁵N (ppm)'
            elif nucleus == '13C-HSQC':
                y_label = '¹³C (ppm)'

        # Formatting
        self.ax_overlay.set_xlabel('¹H (ppm)', fontsize=10)
        self.ax_overlay.set_ylabel(y_label, fontsize=10)
        # Invert axes for NMR convention (high ppm on left/top)
        self.ax_overlay.invert_xaxis()
        self.ax_overlay.invert_yaxis()
        self.ax_overlay.set_title(f'Multi-Spectrum Overlay ({len(visible_spectra)} spectra)', fontsize=11)

        # Add legend
        legend_labels = [spec_data.name for spec_data in visible_spectra]
        legend_colors = [spec_data.color for spec_data in visible_spectra]
        legend_handles = [plt.Line2D([0], [0], color=color, linewidth=2)
                         for color in legend_colors]
        self.ax_overlay.legend(legend_handles, legend_labels, loc='best', fontsize=8)

        # Plot peak markers from selected spectrum (same as spectrum_browser.py)
        self.plot_peak_markers()

        # Restore zoom if it was set before (preserves user's zoom level)
        if xlim_before is not None and ylim_before is not None:
            # Check if limits need to be inverted (after invert_xaxis/yaxis calls)
            # xlim_before and ylim_before are from BEFORE inversion, so apply them directly
            self.ax_overlay.set_xlim(xlim_before)
            self.ax_overlay.set_ylim(ylim_before)
            print(f"   🔍 Restored zoom: X={xlim_before}, Y={ylim_before}")

        self.fig_overlay.tight_layout()
        self.canvas_overlay.draw()

        print(f"✅ Overlay plot updated")

    def export_overlay_plot(self):
        """Export overlay plot as PNG"""
        from tkinter import filedialog

        filename = filedialog.asksaveasfilename(
            defaultextension=".png",
            filetypes=[("PNG files", "*.png"), ("All files", "*.*")],
            title="Export Overlay Plot"
        )

        if filename:
            self.fig_overlay.savefig(filename, dpi=300, bbox_inches='tight')
            print(f"💾 Exported overlay plot to {filename}")
            messagebox.showinfo("Success", f"Overlay plot exported to:\n{filename}")

    def plot_peak_markers(self):
        """Plot peak markers on overlay spectrum (same as spectrum_browser.py)"""
        if not self.reference_peaks:
            return

        try:
            x_coords = []
            y_coords = []
            labels = []
            colors = []

            # Extract peak coordinates and quality info
            for peak in self.reference_peaks:
                # Try multiple possible keys for coordinates
                x_ppm = peak.get('center_x', peak.get('peak_x', peak.get('ppm_x', peak.get('Position_X', 0))))
                y_ppm = peak.get('center_y', peak.get('peak_y', peak.get('ppm_y', peak.get('Position_Y', 0))))
                assignment = peak.get('assignment', peak.get('Assignment', ''))
                quality = peak.get('quality', peak.get('fitting_quality', peak.get('Quality', 'Unknown')))

                if x_ppm and y_ppm and x_ppm != 0 and y_ppm != 0:
                    x_coords.append(float(x_ppm))
                    y_coords.append(float(y_ppm))
                    labels.append(assignment)

                    # Color code by quality (same as spectrum_browser.py)
                    if quality in ['Excellent', 'excellent']:
                        colors.append('lime')
                    elif quality in ['Good', 'good']:
                        colors.append('green')
                    elif quality in ['Fair', 'fair']:
                        colors.append('orange')
                    elif quality in ['Poor', 'poor']:
                        colors.append('red')
                    else:
                        colors.append('yellow')  # Default for unknown

            # Plot peak markers
            if x_coords and y_coords:
                self.ax_overlay.scatter(x_coords, y_coords, c=colors, s=40, marker='o',
                                      alpha=0.8, edgecolors='black', linewidths=1,
                                      label='Peaks', zorder=10)

                # Add peak labels
                for x, y, label in zip(x_coords, y_coords, labels):
                    if label and label != 'Unknown':
                        self.ax_overlay.annotate(label, (x, y), xytext=(5, 5),
                                               textcoords='offset points', fontsize=8,
                                               color='black', weight='bold',
                                               bbox=dict(boxstyle='round,pad=0.2', facecolor='white', alpha=0.8))

                print(f"✅ Plotted {len(x_coords)} peak markers")

        except Exception as e:
            print(f"⚠️ Error plotting peak markers: {e}")

    def setup_voigt_analysis_tab(self):
        """Setup Voigt analysis tab with spectrum selector and peak navigator"""
        # Create horizontal paned window
        voigt_paned = ttk.PanedWindow(self.voigt_tab, orient='horizontal')
        voigt_paned.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)

        # Left panel (75%): Voigt analysis plots
        left_panel = ttk.Frame(voigt_paned)
        voigt_paned.add(left_panel, weight=3)

        # Right panel (25%): Spectrum selector + Peak navigator
        right_panel = ttk.Frame(voigt_paned)
        voigt_paned.add(right_panel, weight=1)

        # Left panel: Voigt plots
        plot_container = ttk.LabelFrame(left_panel, text="📈 Voigt Analysis", padding=5)
        plot_container.pack(fill=tk.BOTH, expand=True)

        # Create 2x1 grid for Voigt analysis plots with compact size
        # Reduced from (8, 6) to (6, 4.5) for better fit
        self.fig_voigt, self.axes_voigt = plt.subplots(2, 1, figsize=(6, 4.5))

        self.canvas_voigt = FigureCanvasTkAgg(self.fig_voigt, plot_container)
        self.canvas_voigt.get_tk_widget().pack(fill=tk.BOTH, expand=True)

        # Add toolbar
        toolbar_voigt = NavigationToolbar2Tk(self.canvas_voigt, plot_container)
        toolbar_voigt.update()

        # Create VoigtAnalysisPlotter instance
        self.voigt_plotter = VoigtAnalysisPlotter(self.fig_voigt, self.axes_voigt)

        # Initialize with placeholder (axes_voigt is 1D array from subplots(2, 1))
        for ax in self.axes_voigt:
            ax.text(0.5, 0.5, 'Select spectrum and peak to view Voigt analysis',
                   ha='center', va='center', transform=ax.transAxes, fontsize='small')
            ax.set_title('Voigt Analysis - No Selection')

        self.canvas_voigt.draw()

        # Right panel: Spectrum selector + Peak navigator
        # Spectrum selector at top
        selector_frame = ttk.LabelFrame(right_panel, text="🔬 Select Spectrum", padding=5)
        selector_frame.pack(fill=tk.X, pady=(0, 10))

        ttk.Label(selector_frame, text="Spectrum:").pack(anchor='w')

        self.spectrum_selector = ttk.Combobox(selector_frame, state='readonly', width=20)
        self.spectrum_selector.pack(fill=tk.X, pady=5)
        self.spectrum_selector.bind('<<ComboboxSelected>>', self.on_spectrum_selected)

        # Populate spectrum selector
        spectrum_names = [spec.name for spec in self.spectra]
        self.spectrum_selector['values'] = spectrum_names
        if spectrum_names:
            self.spectrum_selector.current(0)

        # Peak navigator below
        navigator_frame = ttk.LabelFrame(right_panel, text="📋 Peak Navigator", padding=5)
        navigator_frame.pack(fill=tk.BOTH, expand=True)

        self.voigt_peak_navigator = PeakNavigator(navigator_frame)
        self.voigt_peak_navigator.pack(fill=tk.BOTH, expand=True)

        # Set ourselves as spectrum controller for navigator callbacks
        # Use set_spectrum_controller method (same as spectrum_browser.py)
        self.voigt_peak_navigator.set_spectrum_controller(self)

        # Populate peak navigator for first spectrum (after widgets are created)
        # Schedule this to run after the UI is fully initialized
        if spectrum_names:
            self.window.after(200, self.populate_voigt_for_first_spectrum)

    def populate_voigt_for_first_spectrum(self):
        """Populate Voigt peak navigator for first spectrum on initialization"""
        try:
            print("🔧 Initializing Voigt tab with first spectrum's peaks...")
            # Trigger on_spectrum_selected to load first spectrum's peaks
            self.on_spectrum_selected()
        except Exception as e:
            print(f"⚠️ Failed to populate Voigt tab: {e}")
            import traceback
            traceback.print_exc()

    def setup_voigt_3d_analysis_tab(self):
        """Setup 3D Voigt analysis tab with Peak Navigator (like 2D tab and main GUI)"""
        # Create horizontal paned window (like 2D Voigt tab)
        voigt_3d_paned = ttk.PanedWindow(self.voigt_3d_tab, orient='horizontal')
        voigt_3d_paned.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)

        # Left panel (75%): 3D plot with controls
        left_panel_3d = ttk.Frame(voigt_3d_paned)
        voigt_3d_paned.add(left_panel_3d, weight=3)

        # Right panel (25%): Peak Navigator
        right_panel_3d = ttk.Frame(voigt_3d_paned)
        voigt_3d_paned.add(right_panel_3d, weight=1)

        # Left panel: 3D plot container
        plot_container = ttk.LabelFrame(left_panel_3d, text="🎨 3D Voigt Surface Analysis", padding=5)
        plot_container.pack(fill=tk.BOTH, expand=True)

        # Create control frame at top
        control_frame_3d = ttk.Frame(plot_container)
        control_frame_3d.pack(side=tk.TOP, fill=tk.X, padx=5, pady=3)

        # Row 1: Layer toggling checkboxes
        layer_frame = ttk.LabelFrame(control_frame_3d, text="Layer Visibility", padding=3)
        layer_frame.pack(side=tk.LEFT, padx=3)

        self.show_exp_3d_var = tk.BooleanVar(value=True)
        self.show_fit_3d_var = tk.BooleanVar(value=True)
        self.show_resid_3d_var = tk.BooleanVar(value=True)
        # self.show_cross_3d_var = tk.BooleanVar(value=True)  # Disabled - code kept for future use

        ttk.Checkbutton(layer_frame, text="Experimental", variable=self.show_exp_3d_var,
                        command=lambda: self.voigt_plotter_3d.toggle_experimental(self.show_exp_3d_var.get())
                        ).pack(side=tk.LEFT, padx=2)
        ttk.Checkbutton(layer_frame, text="Fitted", variable=self.show_fit_3d_var,
                        command=lambda: self.voigt_plotter_3d.toggle_fitted(self.show_fit_3d_var.get())
                        ).pack(side=tk.LEFT, padx=2)
        ttk.Checkbutton(layer_frame, text="Residuals", variable=self.show_resid_3d_var,
                        command=lambda: self.voigt_plotter_3d.toggle_residuals(self.show_resid_3d_var.get())
                        ).pack(side=tk.LEFT, padx=2)
        # ttk.Checkbutton(layer_frame, text="Cross-Sections", variable=self.show_cross_3d_var,
        #                 command=lambda: self.voigt_plotter_3d.toggle_cross_sections(self.show_cross_3d_var.get())
        #                 ).pack(side=tk.LEFT, padx=2)  # Disabled - code kept for future use

        # Row 2: Residual mode radio buttons
        residual_frame = ttk.LabelFrame(control_frame_3d, text="Residual Mode", padding=3)
        residual_frame.pack(side=tk.LEFT, padx=3)

        self.residual_mode_3d_var = tk.StringVar(value='overlay')
        ttk.Radiobutton(residual_frame, text="Separate Panel", variable=self.residual_mode_3d_var,
                        value='separate',
                        command=lambda: self.voigt_plotter_3d.set_residual_mode('separate')
                        ).pack(side=tk.LEFT, padx=2)
        ttk.Radiobutton(residual_frame, text="Overlay", variable=self.residual_mode_3d_var,
                        value='overlay',
                        command=lambda: self.voigt_plotter_3d.set_residual_mode('overlay')
                        ).pack(side=tk.LEFT, padx=2)

        # Row 3: Intensity scaling slider
        intensity_frame = ttk.LabelFrame(control_frame_3d, text="Intensity Scale", padding=3)
        intensity_frame.pack(side=tk.LEFT, padx=3, fill=tk.X, expand=True)

        ttk.Label(intensity_frame, text="50%").pack(side=tk.LEFT, padx=2)

        self.intensity_scale_3d_var = tk.DoubleVar(value=100.0)
        intensity_slider_3d = tk.Scale(
            intensity_frame,
            from_=50,
            to=200,
            orient=tk.HORIZONTAL,
            variable=self.intensity_scale_3d_var,
            command=self._on_intensity_scale_change_3d,
            resolution=5,
            showvalue=0,
            length=200
        )
        intensity_slider_3d.pack(side=tk.LEFT, fill=tk.X, expand=True, padx=2)

        self.intensity_scale_label_3d = ttk.Label(intensity_frame, text="100%", width=5)
        self.intensity_scale_label_3d.pack(side=tk.LEFT, padx=2)

        ttk.Label(intensity_frame, text="200%").pack(side=tk.LEFT, padx=2)

        # Create 3D Voigt analysis figure - will be dynamically resized by plotter
        # Reduced from (15, 5) to (10, 4) for better fit in window
        self.fig_voigt_3d = plt.figure(figsize=(10, 4))

        self.canvas_voigt_3d = FigureCanvasTkAgg(self.fig_voigt_3d, plot_container)
        self.canvas_voigt_3d.get_tk_widget().pack(fill=tk.BOTH, expand=True)

        # Add toolbar
        toolbar_voigt_3d = NavigationToolbar2Tk(self.canvas_voigt_3d, plot_container)
        toolbar_voigt_3d.update()

        # Create VoigtAnalysisPlotter instance for 3D (SHARED with main GUI)
        self.voigt_plotter_3d = VoigtAnalysisPlotter(self.fig_voigt_3d, None)

        # Initialize with placeholder (axes_voigt_3d is dynamically created by plotter)
        ax = self.fig_voigt_3d.add_subplot(111)
        ax.axis('off')
        ax.text(0.5, 0.5, 'Select spectrum and peak to view 3D Voigt analysis',
               ha='center', va='center', transform=ax.transAxes, fontsize=12,
               bbox=dict(boxstyle='round,pad=1.0', facecolor='lightgray', alpha=0.9))
        ax.set_title('3D Voigt Analysis - No Peak Selected', fontsize=12, fontweight='bold')

        self.canvas_voigt_3d.draw()

        # Right panel: Spectrum selector + Peak Navigator (same as 2D Voigt tab)
        # Spectrum selector at top
        selector_frame_3d = ttk.LabelFrame(right_panel_3d, text="🔬 Select Spectrum", padding=5)
        selector_frame_3d.pack(fill=tk.X, pady=(0, 10))

        ttk.Label(selector_frame_3d, text="Spectrum:").pack(anchor='w')

        self.spectrum_selector_3d = ttk.Combobox(selector_frame_3d, state='readonly', width=20)
        self.spectrum_selector_3d.pack(fill=tk.X, pady=5)
        self.spectrum_selector_3d.bind('<<ComboboxSelected>>', self.on_spectrum_selected_3d)

        # Populate spectrum selector
        spectrum_names = [spec.name for spec in self.spectra]
        self.spectrum_selector_3d['values'] = spectrum_names
        if spectrum_names:
            self.spectrum_selector_3d.current(0)

        # Peak Navigator below
        navigator_frame_3d = ttk.LabelFrame(right_panel_3d, text="📋 Peak Navigator", padding=5)
        navigator_frame_3d.pack(fill=tk.BOTH, expand=True)

        # Create Peak Navigator instance (share same navigator as 2D tab)
        self.voigt_3d_peak_navigator = PeakNavigator(navigator_frame_3d)
        self.voigt_3d_peak_navigator.pack(fill=tk.BOTH, expand=True)

        # Set ourselves as spectrum controller for navigator callbacks
        self.voigt_3d_peak_navigator.set_spectrum_controller(self)

        # Populate peak navigator for first spectrum (after widgets are created)
        # Schedule this to run after the UI is fully initialized
        if spectrum_names:
            self.window.after(250, self.populate_voigt_3d_for_first_spectrum)

    def populate_voigt_3d_for_first_spectrum(self):
        """Populate 3D Voigt peak navigator for first spectrum on initialization"""
        try:
            print("🔧 Initializing 3D Voigt tab with first spectrum's peaks...")
            # Trigger on_spectrum_selected_3d to load first spectrum's peaks
            self.on_spectrum_selected_3d()
        except Exception as e:
            print(f"⚠️ Failed to populate 3D Voigt tab: {e}")
            import traceback
            traceback.print_exc()

    def _on_intensity_scale_change_3d(self, value):
        """Handle intensity scale slider change for 3D Voigt plot"""
        # Update label
        self.intensity_scale_label_3d.config(text=f"{int(float(value))}%")

        # Cancel pending update if exists (debouncing)
        if hasattr(self, '_scale_update_id_3d'):
            self.window.after_cancel(self._scale_update_id_3d)

        # Schedule update after 100ms of no slider movement
        self._scale_update_id_3d = self.window.after(100,
            lambda: self.voigt_plotter_3d.set_intensity_scale(float(value)))

    def on_spectrum_selected(self, event=None):
        """Handle spectrum selection in 2D Voigt tab"""
        selected_idx = self.spectrum_selector.current()
        if selected_idx < 0:
            return

        self.selected_spectrum_index = selected_idx
        spec_data = self.spectra[selected_idx]

        print(f"🔬 Selected spectrum for 2D Voigt analysis: {spec_data.name}")

        # Load spectrum data if not loaded
        if spec_data.data is None:
            if not spec_data.load_data(self.file_manager):
                messagebox.showerror("Error", f"Failed to load {spec_data.name}")
                return

        # Update peak navigator with this spectrum's peaks
        self.populate_voigt_peak_navigator(spec_data)

    def on_spectrum_selected_3d(self, event=None):
        """Handle spectrum selection in 3D Voigt tab"""
        selected_idx = self.spectrum_selector_3d.current()
        if selected_idx < 0:
            return

        self.selected_spectrum_index = selected_idx
        spec_data = self.spectra[selected_idx]

        print(f"🔬 Selected spectrum for 3D Voigt analysis: {spec_data.name}")

        # Load spectrum data if not loaded
        if spec_data.data is None:
            if not spec_data.load_data(self.file_manager):
                messagebox.showerror("Error", f"Failed to load {spec_data.name}")
                return

        # Update peak navigator with this spectrum's peaks
        self.populate_voigt_peak_navigator(spec_data)

    def populate_voigt_peak_navigator(self, spec_data):
        """Populate peak navigators (2D and 3D tabs) with peaks from selected spectrum"""
        if not spec_data.fitted_peaks:
            print(f"⚠️ No fitted peaks for {spec_data.name}")
            # Load empty list to clear navigators
            self.voigt_peak_navigator.load_detected_peaks([])
            self.voigt_peak_navigator.refresh_peak_list()
            if hasattr(self, 'voigt_3d_peak_navigator'):
                self.voigt_3d_peak_navigator.load_detected_peaks([])
                self.voigt_3d_peak_navigator.refresh_peak_list()
            return

        # Debug: Print first peak structure to see available keys
        if len(spec_data.fitted_peaks) > 0:
            print(f"🔍 DEBUG Voigt: First peak keys: {list(spec_data.fitted_peaks[0].keys())}")

        # Update both 2D and 3D peak navigators
        self.voigt_peak_navigator.load_detected_peaks(spec_data.fitted_peaks)
        self.voigt_peak_navigator.refresh_peak_list()

        if hasattr(self, 'voigt_3d_peak_navigator'):
            self.voigt_3d_peak_navigator.load_detected_peaks(spec_data.fitted_peaks)
            self.voigt_3d_peak_navigator.refresh_peak_list()

        print(f"📋 Populated peak navigators: {len(spec_data.fitted_peaks)} peaks from {spec_data.name}")

    def center_on_coordinates(self, x, y):
        """
        Controller interface method called by PeakNavigator when peak selected.
        For multi-spectrum viewer, we don't center the main overlay plot.
        """
        # Store selected coordinates
        self.selected_peak_index = None  # Will be updated by set_selected_peak
        print(f"📍 Peak selected at ({x:.3f}, {y:.1f})")

    def set_selected_peak(self, peak_index, peak_type, source="navigator"):
        """
        Controller interface method called by PeakNavigator.
        Update our selected peak index for Voigt analysis.
        """
        self.selected_peak_index = peak_index
        print(f"📍 Set selected peak index: {peak_index}")

        # Auto-analyze when peak selected
        if peak_index is not None:
            self.plot_selected_peak_voigt()

    def navigator_show_peak_analysis(self, peak_type, peak_index):
        """Show Voigt analysis for selected peak (called by peak navigator analyze button)"""
        if peak_index is not None:
            self.set_selected_peak(peak_index, peak_type, source="navigator")

    def plot_selected_peak_voigt(self):
        """Plot Voigt analysis for currently selected peak"""
        if self.selected_peak_index is None:
            print("⚠️ No peak selected")
            return

        spec_data = self.spectra[self.selected_spectrum_index]

        # Check if fitted_peaks exists and is not None
        if not spec_data.fitted_peaks:
            messagebox.showwarning("No Peaks", "No fitted peaks available for selected spectrum")
            return

        # Check peak index is valid
        if self.selected_peak_index >= len(spec_data.fitted_peaks):
            messagebox.showerror("Error", "Invalid peak index")
            return

        peak = spec_data.fitted_peaks[self.selected_peak_index]

        print(f"🔬 Analyzing peak {peak.get('assignment', 'Unknown')} from {spec_data.name}")

        # Plot Voigt analysis
        self.plot_voigt_analysis(spec_data, peak)

    def plot_voigt_analysis(self, spec_data, peak):
        """Plot Voigt analysis for selected peak from selected spectrum (same as spectrum_browser.py)"""
        if spec_data.integrator is None:
            messagebox.showerror("Error", "Spectrum data not loaded")
            return

        # Extract peak assignment
        assignment = peak.get('assignment', peak.get('Assignment', 'Unknown'))

        # CRITICAL: Find stored Voigt result with full analysis data (same as spectrum_browser.py)
        stored_result = self.find_stored_voigt_result(spec_data, assignment)

        # Plot using VoigtAnalysisPlotter
        try:
            if stored_result:
                # VoigtAnalysisPlotter.plot_voigt_analysis() takes complete peak result dict
                self.voigt_plotter.plot_voigt_analysis(stored_result)

                # Also update 3D view
                self.voigt_plotter_3d.plot_voigt_analysis_3d(stored_result)
                self.canvas_voigt_3d.draw()

                # Update overall figure title to show spectrum name
                # (works for both 1D and 2D layouts)
                self.fig_voigt.suptitle(
                    f"Voigt Analysis - {spec_data.name} - Peak {assignment}",
                    fontsize=11, fontweight='bold'
                )

                print(f"✅ Plotted Voigt analysis for peak {assignment} from {spec_data.name}")
            else:
                # No stored result - show placeholder
                for ax in self.axes_voigt:
                    ax.clear()
                    ax.text(0.5, 0.5, "❌ No Voigt fitting results available",
                           ha='center', va='center', transform=ax.transAxes, fontsize='small')
                    ax.set_title('Voigt Analysis - No Results')

                print(f"⚠️ No stored Voigt result found for peak {assignment}")

            self.canvas_voigt.draw()

        except Exception as e:
            print(f"❌ Error plotting Voigt analysis: {e}")
            import traceback
            traceback.print_exc()

            # Show error in plots
            for ax in self.axes_voigt:
                ax.clear()
                ax.text(0.5, 0.5, f"Error: {str(e)[:50]}",
                       ha='center', va='center', transform=ax.transAxes, fontsize='small')
            self.canvas_voigt.draw()

    def find_stored_voigt_result(self, spec_data, assignment):
        """
        Find stored Voigt fitting results for a peak by assignment
        Returns the detailed fitted_results data from series processing
        (Same logic as spectrum_browser.py)
        """
        if not assignment:
            return None

        # First check if we have fitted_results from result_data (detailed Voigt analysis)
        # This is where the complete region_2d and fit data is stored
        fitted_results = spec_data.result_data.get('fitted_results', [])
        if fitted_results:
            for result in fitted_results:
                # Try both assignment key variations
                result_assignment = result.get('assignment') or result.get('Assignment')
                if result_assignment == assignment:
                    r_squared = result.get('r_squared') or result.get('R_Squared') or result.get('avg_r_squared') or 0
                    print(f"✅ Found stored Voigt result for {assignment}: R²={r_squared:.3f}")
                    return result

        # Fallback: check if integrator has fitted_peaks
        if hasattr(spec_data.integrator, 'fitted_peaks') and spec_data.integrator.fitted_peaks:
            for result in spec_data.integrator.fitted_peaks:
                if result.get('assignment') == assignment:
                    print(f"✅ Found integrator Voigt result for {assignment}")
                    return result

        print(f"⚠️ No stored Voigt result found for {assignment}")
        return None


def open_multi_spectrum_viewer(parent, file_manager, all_results):
    """
    Open multi-spectrum overlay viewer window

    Parameters:
    -----------
    parent : tk.Tk
        Parent window
    file_manager : NMRFileManager
        File manager for loading NMR spectra
    all_results : list
        List of result dictionaries for all spectra
    """
    try:
        viewer = MultiSpectrumViewer(parent, file_manager, all_results)
        return viewer
    except Exception as e:
        print(f"❌ Error opening multi-spectrum viewer: {e}")
        import traceback
        traceback.print_exc()
        messagebox.showerror("Error", f"Failed to open multi-spectrum viewer:\n{str(e)}")
        return None

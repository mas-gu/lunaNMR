#!/usr/bin/env python3
"""
Spectrum Browser Module

This module provides functionality for browsing and selecting individual spectra
from series integration results for quality control and detailed analysis.

Classes:
- SpectrumBrowserDialog: Main spectrum selection dialog
- SpectrumViewer: Individual spectrum analysis window
- PeakAnalysisPanel: Peak properties and Voigt analysis display

Author: Guillaume Mas
Date: 2025
"""

import os
import sys
import tkinter as tk
from tkinter import ttk, messagebox
import pandas as pd
import numpy as np
from pathlib import Path

# Add current directory to path for imports
current_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, current_dir)

try:
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
    from matplotlib.figure import Figure
    MATPLOTLIB_AVAILABLE = True
except ImportError:
    MATPLOTLIB_AVAILABLE = False
    print("Warning: matplotlib not available for spectrum browser")

try:
    from lunaNMR.core.core_integrator import EnhancedVoigtIntegrator
    from lunaNMR.gui.visualization import SpectrumPlotter, VoigtAnalysisPlotter
    from lunaNMR.utils.file_manager import NMRFileManager
    CORE_MODULES_AVAILABLE = True
except ImportError:
    print("Warning: Some modules not available for spectrum browser")
    CORE_MODULES_AVAILABLE = False

    # Create fallback classes for testing
    class EnhancedVoigtIntegrator:
        def __init__(self):
            self.nmr_data = None
            self.fitted_peaks = []
            self.integration_results = []
            self.enhanced_fitter = None
            self.fitting_parameters = {}
            print("⚠️ Using fallback integrator - enhanced features not available")

        def load_nmr_file(self, file_path):
            return False

    class SpectrumPlotter:
        def __init__(self, fig, ax):
            self.fig = fig
            self.ax = ax

        def plot_spectrum(self, integrator):
            pass

        def plot_peaks(self, integrator, show_detected=True, show_assigned=True):
            pass

class SpectrumBrowserDialog:
    """Dialog for browsing and selecting individual spectra from series results"""

    def __init__(self, parent, batch_results, series_processor, original_data_folder=None):
        self.parent = parent
        self.batch_results = batch_results
        self.original_data_folder = original_data_folder
        self.series_processor = series_processor
        self.selected_spectrum = None

        # Create dialog window
        self.dialog = tk.Toplevel(parent)
        self.dialog.title("Browse Individual Spectra - Quality Control")
        self.dialog.geometry("900x600")
        self.dialog.transient(parent)
        # Remove grab_set() to allow interaction with spectrum viewer windows
        # self.dialog.grab_set()  # This was blocking spectrum viewer buttons

        self.setup_dialog()
        self.populate_spectrum_list()

    def setup_dialog(self):
        """Setup the dialog layout"""
        # Main frame
        main_frame = ttk.Frame(self.dialog)
        main_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

        # Title and info
        title_frame = ttk.Frame(main_frame)
        title_frame.pack(fill=tk.X, pady=(0, 10))

        ttk.Label(title_frame, text="📊 Individual Spectrum Browser",
                 font=('TkDefaultFont', 14, 'bold')).pack(side=tk.LEFT)

        if self.batch_results:
            summary = self.batch_results.get_summary()
            info_text = f"Series: {summary['total_spectra']} spectra, {summary['successful']} successful ({summary['success_rate']:.1f}%)"
            ttk.Label(title_frame, text=info_text,
                     font=('TkDefaultFont', 10), foreground='blue').pack(side=tk.RIGHT)

        # Search and filter frame
        filter_frame = ttk.LabelFrame(main_frame, text="🔍 Filter & Search", padding=10)
        filter_frame.pack(fill=tk.X, pady=(0, 10))

        # Search box
        search_frame = ttk.Frame(filter_frame)
        search_frame.pack(fill=tk.X)

        ttk.Label(search_frame, text="Search:").pack(side=tk.LEFT)
        self.search_var = tk.StringVar()
        self.search_var.trace('w', self.filter_spectra)
        search_entry = ttk.Entry(search_frame, textvariable=self.search_var, width=30)
        search_entry.pack(side=tk.LEFT, padx=(5, 20))

        # Quality filter
        ttk.Label(search_frame, text="Quality:").pack(side=tk.LEFT)
        self.quality_filter = tk.StringVar(value="All")
        quality_combo = ttk.Combobox(search_frame, textvariable=self.quality_filter,
                                   values=["All", "Excellent", "Good", "Fair", "Poor", "Failed"],
                                   width=12, state="readonly")
        quality_combo.pack(side=tk.LEFT, padx=(5, 20))
        quality_combo.bind('<<ComboboxSelected>>', self.filter_spectra)

        # Status filter
        ttk.Label(search_frame, text="Status:").pack(side=tk.LEFT)
        self.status_filter = tk.StringVar(value="All")
        status_combo = ttk.Combobox(search_frame, textvariable=self.status_filter,
                                  values=["All", "Success", "Failed"],
                                  width=10, state="readonly")
        status_combo.pack(side=tk.LEFT, padx=(5, 0))
        status_combo.bind('<<ComboboxSelected>>', self.filter_spectra)

        # Spectrum list frame
        list_frame = ttk.LabelFrame(main_frame, text="📋 Spectrum List", padding=10)
        list_frame.pack(fill=tk.BOTH, expand=True, pady=(0, 10))

        # Create treeview for spectrum list
        columns = ('Status', 'Spectrum', 'Type', 'Detected', 'Total', 'Detection Rate', 'Quality', 'Processing Time')
        self.tree = ttk.Treeview(list_frame, columns=columns, show='headings', height=15)

        # Configure column headings and widths
        column_widths = {'Status': 60, 'Spectrum': 200, 'Type': 80, 'Detected': 80,
                        'Total': 60, 'Detection Rate': 100, 'Quality': 80, 'Processing Time': 100}

        for col in columns:
            self.tree.heading(col, text=col, command=lambda c=col: self.sort_by_column(c))
            self.tree.column(col, width=column_widths.get(col, 100))

        # Add scrollbars
        v_scrollbar = ttk.Scrollbar(list_frame, orient=tk.VERTICAL, command=self.tree.yview)
        h_scrollbar = ttk.Scrollbar(list_frame, orient=tk.HORIZONTAL, command=self.tree.xview)
        self.tree.configure(yscrollcommand=v_scrollbar.set, xscrollcommand=h_scrollbar.set)

        # Pack treeview and scrollbars
        self.tree.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        v_scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
        h_scrollbar.pack(side=tk.BOTTOM, fill=tk.X)

        # Bind double-click event
        self.tree.bind('<Double-1>', self.on_spectrum_double_click)
        self.tree.bind('<Button-1>', self.on_spectrum_select)

        # Buttons frame
        button_frame = ttk.Frame(main_frame)
        button_frame.pack(fill=tk.X)

        ttk.Button(button_frame, text="📊 Open Spectrum Viewer",
                  command=self.open_spectrum_viewer).pack(side=tk.LEFT, padx=(0, 10))
        ttk.Button(button_frame, text="📈 Quick Analysis",
                  command=self.show_quick_analysis).pack(side=tk.LEFT, padx=(0, 10))
        ttk.Button(button_frame, text="🔄 Refresh List",
                  command=self.refresh_list).pack(side=tk.LEFT, padx=(0, 20))

        ttk.Button(button_frame, text="Close",
                  command=self.dialog.destroy).pack(side=tk.RIGHT)

    def populate_spectrum_list(self):
        """Populate the spectrum list with data from batch results"""
        if not self.batch_results or not hasattr(self.batch_results, 'results'):
            return

        # Clear existing items
        for item in self.tree.get_children():
            self.tree.delete(item)

        # Add spectrum data
        for spectrum_name, result in self.batch_results.results.items():
            # Calculate quality score
            quality_score = self.calculate_quality_score(result)
            quality_text = self.get_quality_text(quality_score)

            # Status icon
            status_icon = "✅" if result.get('status') == 'success' else "❌"

            # Spectrum type
            spectrum_type = "📍 REF" if result.get('is_reference', False) else "📄"

            # Processing time
            proc_time = result.get('processing_time', 0)
            proc_time_str = f"{proc_time:.2f}s" if proc_time else "N/A"

            # Insert row
            item_id = self.tree.insert('', tk.END, values=(
                status_icon,
                spectrum_name,
                spectrum_type,
                result.get('detected_peaks', 0),
                result.get('total_peaks', 0),
                f"{result.get('detection_rate', 0):.1f}%",
                quality_text,
                proc_time_str
            ))

            # Color coding based on quality and status
            if result.get('status') == 'failed':
                self.tree.set(item_id, 'Status', '❌')
                tags = ['failed']
            elif quality_score >= 90:
                tags = ['excellent']
            elif quality_score >= 70:
                tags = ['good']
            elif quality_score >= 50:
                tags = ['fair']
            else:
                tags = ['poor']

            # Apply tags for color coding
            for tag in tags:
                self.tree.item(item_id, tags=(tag,))

        # Configure tag colors
        self.tree.tag_configure('failed', background='#ffebee')
        self.tree.tag_configure('poor', background='#fff3e0')
        self.tree.tag_configure('fair', background='#fffde7')
        self.tree.tag_configure('good', background='#f1f8e9')
        self.tree.tag_configure('excellent', background='#e8f5e8')

    def calculate_quality_score(self, result):
        """Calculate quality score for a spectrum (0-100)"""
        if result.get('status') == 'failed':
            return 0

        detection_rate = result.get('detection_rate', 0)
        detected_peaks = result.get('detected_peaks', 0)

        # Base score from detection rate
        score = detection_rate

        # Bonus for having detected peaks
        if detected_peaks > 0:
            score += min(10, detected_peaks / 10)  # Up to 10 bonus points

        # Penalty for very low detection
        if detection_rate < 50:
            score *= 0.8

        return min(100, max(0, score))

    def get_quality_text(self, score):
        """Convert quality score to text"""
        if score >= 90:
            return "Excellent"
        elif score >= 70:
            return "Good"
        elif score >= 50:
            return "Fair"
        elif score > 0:
            return "Poor"
        else:
            return "Failed"

    def filter_spectra(self, *args):
        """Filter spectrum list based on search and filter criteria"""
        search_text = self.search_var.get().lower()
        quality_filter = self.quality_filter.get()
        status_filter = self.status_filter.get()

        # Clear and repopulate with filtered items
        for item in self.tree.get_children():
            self.tree.delete(item)

        if not self.batch_results:
            return

        for spectrum_name, result in self.batch_results.results.items():
            # Apply search filter
            if search_text and search_text not in spectrum_name.lower():
                continue

            # Apply status filter
            if status_filter != "All":
                if status_filter == "Success" and result.get('status') != 'success':
                    continue
                if status_filter == "Failed" and result.get('status') != 'failed':
                    continue

            # Apply quality filter
            if quality_filter != "All":
                quality_score = self.calculate_quality_score(result)
                quality_text = self.get_quality_text(quality_score)
                if quality_filter != quality_text:
                    continue

            # Add item (same logic as populate_spectrum_list)
            quality_score = self.calculate_quality_score(result)
            quality_text = self.get_quality_text(quality_score)
            status_icon = "✅" if result.get('status') == 'success' else "❌"
            spectrum_type = "📍 REF" if result.get('is_reference', False) else "📄"
            proc_time = result.get('processing_time', 0)
            proc_time_str = f"{proc_time:.2f}s" if proc_time else "N/A"

            item_id = self.tree.insert('', tk.END, values=(
                status_icon, spectrum_name, spectrum_type,
                result.get('detected_peaks', 0), result.get('total_peaks', 0),
                f"{result.get('detection_rate', 0):.1f}%", quality_text, proc_time_str
            ))

            # Apply color coding
            if result.get('status') == 'failed':
                tags = ['failed']
            elif quality_score >= 90:
                tags = ['excellent']
            elif quality_score >= 70:
                tags = ['good']
            elif quality_score >= 50:
                tags = ['fair']
            else:
                tags = ['poor']

            for tag in tags:
                self.tree.item(item_id, tags=(tag,))

    def sort_by_column(self, column):
        """Sort the tree by the selected column"""
        items = [(self.tree.set(item, column), item) for item in self.tree.get_children('')]

        # Sort numerically for numeric columns
        if column in ['Detected', 'Total', 'Detection Rate', 'Processing Time']:
            try:
                items.sort(key=lambda x: float(x[0].replace('%', '').replace('s', '')) if x[0] not in ['N/A', ''] else 0)
            except:
                items.sort()
        else:
            items.sort()

        # Rearrange items
        for index, (val, item) in enumerate(items):
            self.tree.move(item, '', index)

    def on_spectrum_select(self, event):
        """Handle spectrum selection"""
        selection = self.tree.selection()
        if selection:
            item = selection[0]
            spectrum_name = self.tree.item(item, 'values')[1]  # Spectrum name is in column 1
            self.selected_spectrum = spectrum_name

    def on_spectrum_double_click(self, event):
        """Handle double-click on spectrum to open viewer"""
        self.open_spectrum_viewer()

    def open_spectrum_viewer(self):
        """Open the spectrum viewer for the selected spectrum"""
        if not self.selected_spectrum:
            messagebox.showwarning("No Selection", "Please select a spectrum first.")
            return

        try:
            # Get spectrum data
            if self.selected_spectrum not in self.batch_results.results:
                messagebox.showerror("Error", f"Spectrum data not found for {self.selected_spectrum}")
                return

            raw_result_data = self.batch_results.results[self.selected_spectrum]

            # Find the spectrum file path FIRST
            spectrum_file_path = self.find_spectrum_file_path(self.selected_spectrum)
            if not spectrum_file_path:
                messagebox.showerror("Error", f"Could not find spectrum file for {self.selected_spectrum}")
                return

            # Create spectrum viewer with correct data
            viewer = SpectrumViewer(self.dialog, self.selected_spectrum, spectrum_file_path,
                                  raw_result_data, self.series_processor, self.original_data_folder)

        except Exception as e:
            messagebox.showerror("Error", f"Failed to open spectrum viewer:\n{str(e)}")


    def _prepare_spectrum_viewer_data(self, raw_data):
        """Prepare enhanced data for spectrum viewer with proper coordinates"""
        enhanced_data = raw_data.copy()

        # Get integration results
        integration_results = enhanced_data.get('integration_results', [])

        print(f"🔧 Preparing spectrum viewer data:")
        print(f"   Raw data keys: {list(raw_data.keys())}")
        print(f"   Integration results count: {len(integration_results)}")

        if integration_results:
            # Verify and enhance coordinate data
            for i, peak in enumerate(integration_results):
                if i < 3:  # Debug first 3 peaks
                    print(f"   Peak {i+1} data keys: {list(peak.keys())}")
                    pos_x = peak.get('Position_X', peak.get('ppm_x', 0))
                    pos_y = peak.get('Position_Y', peak.get('ppm_y', 0))
                    print(f"   Peak {i+1} coordinates: X={pos_x}, Y={pos_y}")

                # Ensure all required fields exist
                if 'Position_X' not in peak and 'ppm_x' in peak:
                    peak['Position_X'] = peak['ppm_x']
                if 'Position_Y' not in peak and 'ppm_y' in peak:
                    peak['Position_Y'] = peak['ppm_y']

                # Add visualization fields if missing
                if 'detected' not in peak:
                    peak['detected'] = True
                if 'fitted' not in peak:
                    peak['fitted'] = True
                if 'success' not in peak:
                    peak['success'] = True

        # Also check fitted_results as fallback
        fitted_results = enhanced_data.get('fitted_results', [])
        if fitted_results and not integration_results:
            print(f"   Using fitted_results as fallback ({len(fitted_results)} peaks)")
            # Convert fitted_results to integration format on-the-fly
            enhanced_integration_results = []
            for result in fitted_results:
                if isinstance(result, dict):
                    # STANDARDIZED dual-format data structure (primary and legacy)
                    enhanced_peak = {
                        # Primary fields (from standard format)
                        'assignment': result.get('assignment') or result.get('Assignment') or 'Unknown',
                        'peak_number': result.get('peak_number') or result.get('Peak_Number') or 0,
                        'ppm_x': result.get('ppm_x') or result.get('Position_X') or 0,
                        'ppm_y': result.get('ppm_y') or result.get('Position_Y') or 0,
                        'height': result.get('height') or result.get('Height') or 0,
                        'volume': result.get('volume') or result.get('Volume') or 0,
                        'snr': result.get('snr') or result.get('SNR') or 0,
                        'quality': result.get('quality') or result.get('Quality') or result.get('fitting_quality') or 'Unknown',
                        'r_squared': result.get('r_squared') or result.get('R_Squared') or result.get('avg_r_squared') or 0,

                        # Legacy fields (for backward compatibility)
                        'Assignment': result.get('assignment') or result.get('Assignment') or 'Unknown',
                        'Position_X': result.get('ppm_x') or result.get('Position_X') or 0,
                        'Position_Y': result.get('ppm_y') or result.get('Position_Y') or 0,
                        'Height': result.get('height') or result.get('Height') or 0,
                        'Volume': result.get('volume') or result.get('Volume') or 0,
                        'SNR': result.get('snr') or result.get('SNR') or 0,
                        'Quality': result.get('quality') or result.get('Quality') or result.get('fitting_quality') or 'Unknown',
                        'R_Squared': result.get('r_squared') or result.get('R_Squared') or result.get('avg_r_squared') or 0,
                        'Peak_Number': result.get('peak_number') or result.get('Peak_Number') or 0,
                    #enhanced_peak = {
                    #    'Assignment': result.get('assignment', 'Unknown'),
                    #    'Position_X': result.get('ppm_x', result.get('Position_X', 0)),
                    #    'Position_Y': result.get('ppm_y', result.get('Position_Y', 0)),
                    #    'Height': result.get('height', result.get('Height', 0)),
                    #    'Volume': result.get('integration_volume', result.get('Volume', 0)),
                    #    'SNR': result.get('snr', result.get('SNR', 0)),
                    #    'Quality': result.get('fitting_quality', result.get('Quality', 'Unknown')),
                    #    'R_Squared': result.get('avg_r_squared', result.get('R_Squared', 0)),
                    #    'Peak_Number': result.get('peak_number', result.get('Peak_Number', 0)),
                        'detected': True,
                        'fitted': True,
                        'success': True
                    }
                    enhanced_integration_results.append(enhanced_peak)

            enhanced_data['integration_results'] = enhanced_integration_results
            print(f"   ✅ Converted {len(enhanced_integration_results)} fitted_results to integration format")

        return enhanced_data



    def find_spectrum_file_path(self, spectrum_name):
        """Find the full file path for a spectrum"""
        print(f"🔍 Looking for spectrum: {spectrum_name}")

        # PRIORITY 1: Use the original data folder if available
        if self.original_data_folder:
            direct_path = os.path.join(self.original_data_folder, spectrum_name)
            print(f"🔍 Checking original data folder: {direct_path}")
            if os.path.exists(direct_path):
                print(f"✅ Found spectrum in original data folder: {direct_path}")
                return direct_path

            # Also search subdirectories of the original data folder
            for root, dirs, files in os.walk(self.original_data_folder):
                if spectrum_name in files:
                    found_path = os.path.join(root, spectrum_name)
                    print(f"✅ Found spectrum in subdirectory: {found_path}")
                    return found_path

        # PRIORITY 2: Check batch results metadata for stored file paths
        if hasattr(self.batch_results, 'metadata') and 'processing_files' in self.batch_results.metadata:
            print(f"🔍 Checking metadata processing_files...")
            for file_path in self.batch_results.metadata['processing_files']:
                if os.path.basename(file_path) == spectrum_name:
                    if os.path.exists(file_path):
                        print(f"✅ Found spectrum from metadata: {file_path}")
                        return file_path
                    else:
                        print(f"⚠️ File in metadata doesn't exist: {file_path}")

        # PRIORITY 3: Search backwards from series processor output folder (fallback)
        if self.series_processor and hasattr(self.series_processor, 'output_folder'):
            output_folder = self.series_processor.output_folder
            if output_folder:
                base_path = os.path.dirname(output_folder)
                search_paths = [
                    base_path,  # Direct parent
                    os.path.dirname(base_path),  # Grandparent
                    os.path.dirname(os.path.dirname(base_path))  # Great-grandparent
                ]

                print(f"🔍 Fallback search in output folder hierarchy...")
                for search_path in search_paths:
                    for root, dirs, files in os.walk(search_path):
                        if spectrum_name in files:
                            found_path = os.path.join(root, spectrum_name)
                            print(f"✅ Found spectrum via fallback search: {found_path}")
                            return found_path

        print(f"❌ Could not find spectrum file: {spectrum_name}")
        return None

    def show_quick_analysis(self):
        """Show quick analysis for selected spectrum"""
        if not self.selected_spectrum:
            messagebox.showwarning("No Selection", "Please select a spectrum first.")
            return

        result_data = self.batch_results.results.get(self.selected_spectrum)
        if not result_data:
            return

        # Create quick analysis dialog
        analysis_window = tk.Toplevel(self.dialog)
        analysis_window.title(f"Quick Analysis - {self.selected_spectrum}")
        analysis_window.geometry("500x400")
        # Allow analysis window to be independent for better interaction
        # analysis_window.transient(self.dialog)

        # Analysis content
        text_widget = tk.Text(analysis_window, wrap=tk.WORD, font=('Courier', 10))
        scrollbar = ttk.Scrollbar(analysis_window, orient=tk.VERTICAL, command=text_widget.yview)
        text_widget.configure(yscrollcommand=scrollbar.set)

        # Generate analysis text
        analysis_text = self.generate_quick_analysis_text(result_data)
        text_widget.insert(tk.END, analysis_text)
        text_widget.config(state=tk.DISABLED)

        text_widget.pack(side=tk.LEFT, fill=tk.BOTH, expand=True, padx=10, pady=10)
        scrollbar.pack(side=tk.RIGHT, fill=tk.Y, pady=10)

    def generate_quick_analysis_text(self, result_data):
        """Generate quick analysis text for a spectrum"""
        quality_score = self.calculate_quality_score(result_data)
        quality_text = self.get_quality_text(quality_score)

        text = f"""QUICK ANALYSIS REPORT
{'='*50}

Spectrum: {self.selected_spectrum}
Status: {result_data.get('status', 'Unknown').title()}
Overall Quality: {quality_text} ({quality_score:.1f}/100)

DETECTION SUMMARY:
• Total Reference Peaks: {result_data.get('total_peaks', 0)}
• Detected Peaks: {result_data.get('detected_peaks', 0)}
• Detection Rate: {result_data.get('detection_rate', 0):.1f}%
• Processing Time: {result_data.get('processing_time', 0):.2f} seconds

QUALITY METRICS:
• Spectrum Type: {'Reference' if result_data.get('is_reference', False) else 'Sample'}
• Noise Level: {result_data.get('noise_level', 0):.2e}

"""

        if result_data.get('status') == 'failed':
            text += f"ERROR DETAILS:\n• {result_data.get('error', 'Unknown error')}\n\n"

        if result_data.get('integration_results'):
            integrations = result_data['integration_results']
            if integrations:
                text += f"INTEGRATION RESULTS:\n"
                text += f"• Total Integrations: {len(integrations)}\n"

                # Quality breakdown
                #good_quality = sum(1 for i in integrations if i.get('SNR', 0) >= 5)
                #text += f"• Good Quality (SNR ≥ 5): {good_quality}\n"
                #text += f"• Quality Rate: {good_quality/len(integrations)*100:.1f}%\n\n"

                # Top 5 peaks by intensity
                #sorted_peaks = sorted(integrations, key=lambda x: x.get('Height', 0), reverse=True)
                #text += "TOP PEAKS BY INTENSITY:\n"
                #for i, peak in enumerate(sorted_peaks[:5], 1):
                #    text += f"  {i}. {peak.get('Assignment', 'Unknown')}: "
                #    text += f"Height={peak.get('Height', 0):.2e}, SNR={peak.get('SNR', 0):.1f}\n"

                # STANDARDIZED quality breakdown (primary fields first)
                good_quality = sum(1 for i in integrations if (i.get('snr') or i.get('SNR') or 0) >= 5)
                text += f"• Good Quality (SNR ≥ 5): {good_quality}\n"
                text += f"• Quality Rate: {good_quality/len(integrations)*100:.1f}%\n\n"

                # STANDARDIZED top 5 peaks by intensity (primary fields first)
                sorted_peaks = sorted(integrations, key=lambda x: (x.get('height') or x.get('Height') or 0), reverse=True)
                text += "TOP PEAKS BY INTENSITY:\n"
                for i, peak in enumerate(sorted_peaks[:5], 1):
                    assignment = peak.get('assignment') or peak.get('Assignment') or 'Unknown'
                    height = peak.get('height') or peak.get('Height') or 0
                    snr = peak.get('snr') or peak.get('SNR') or 0
                    text += f"  {i}. {assignment}: Height={height:.2e}, SNR={snr:.1f}\n"

        return text

    def refresh_list(self):
        """Refresh the spectrum list"""
        self.populate_spectrum_list()


class SpectrumViewer:
    """Individual spectrum analysis window with interactive peak analysis"""

    def __init__(self, parent, spectrum_name, spectrum_file_path, result_data, series_processor, original_data_folder=None):
        self.parent = parent
        self.spectrum_name = spectrum_name
        self.original_data_folder = original_data_folder
        self.spectrum_file_path = spectrum_file_path
        self.result_data = result_data
        self.series_processor = series_processor
        self.selected_peak = None
        self.integrator = None

        # Contour control variables (user-requested default values)
        self.contour_levels = tk.IntVar(value=40)  # User requested: levels 40
        self.contour_min = tk.DoubleVar(value=0.15)  # User requested: min level 0.15
        self.contour_increment = tk.DoubleVar(value=1.1)  # User requested: increment 1.1
        self.colormap_var = tk.StringVar(value='viridis')

        # Performance and interaction state
        self.plot_needs_update = True
        self.zoom_history = []
        self.current_zoom = None
        self.plot_cache = None

        # Create viewer window
        self.window = tk.Toplevel(parent)
        self.window.title(f"Spectrum Viewer - {spectrum_name}")
        self.window.geometry("1400x900")
        self.window.minsize(1200, 800)
        # Make spectrum viewer independent - not transient to allow button interaction
        # self.window.transient(parent)  # This was preventing button clicks

        # Ensure window is properly displayed and can receive focus
        self.window.update_idletasks()
        self.window.lift()  # Bring to front
        self.window.focus_force()  # Ensure it can receive input

        self.setup_viewer()

        # Force window update to be responsive first
        self.window.update_idletasks()
        self.window.update()

        # Load data with minimal delay for UI responsiveness
        self.window.after(50, self.load_spectrum_data_proper)

    def setup_viewer(self):
        """Setup the spectrum viewer layout with tabbed interface like main GUI"""
        # TABBED INTERFACE: Create notebook with tabs like main GUI
        self.notebook = ttk.Notebook(self.window)
        self.notebook.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)

        # Tab 1: Main Spectrum (current functionality preserved)
        self.main_tab = ttk.Frame(self.notebook)
        self.notebook.add(self.main_tab, text="📊 Main Spectrum")
        self.setup_main_spectrum_tab()

        # Tab 2: Voigt Analysis (integrated directly instead of popup window)
        self.voigt_tab = ttk.Frame(self.notebook)
        self.notebook.add(self.voigt_tab, text="📈 Voigt Analysis")
        self.setup_voigt_analysis_tab()

    def setup_main_spectrum_tab(self):
        """Setup the main spectrum tab (preserves current functionality)"""
        # ENHANCED: Create scrollable main container for small screen compatibility
        canvas_container = tk.Canvas(self.main_tab, highlightthickness=0)
        scrollbar = ttk.Scrollbar(self.main_tab, orient="vertical", command=canvas_container.yview)
        scrollable_frame = ttk.Frame(canvas_container)

        # Configure scrolling
        scrollable_frame.bind(
            "<Configure>",
            lambda e: canvas_container.configure(scrollregion=canvas_container.bbox("all"))
        )

        canvas_container.create_window((0, 0), window=scrollable_frame, anchor="nw")
        canvas_container.configure(yscrollcommand=scrollbar.set)

        # Pack scrollable components
        canvas_container.pack(side="left", fill="both", expand=True, padx=5, pady=5)
        scrollbar.pack(side="right", fill="y")

        # Store references for later use
        self.scrollable_canvas = canvas_container
        self.main_frame = scrollable_frame

        # Configure tab grid weights for proper resizing
        self.main_tab.columnconfigure(0, weight=1)
        self.main_tab.rowconfigure(0, weight=1)

        # Bind mousewheel to canvas for scroll support
        def _on_mousewheel(event):
            canvas_container.yview_scroll(int(-1*(event.delta/120)), "units")

        def _bind_to_mousewheel(event):
            canvas_container.bind_all("<MouseWheel>", _on_mousewheel)

        def _unbind_from_mousewheel(event):
            canvas_container.unbind_all("<MouseWheel>")

        canvas_container.bind('<Enter>', _bind_to_mousewheel)
        canvas_container.bind('<Leave>', _unbind_from_mousewheel)

        # Title bar
        title_frame = ttk.Frame(self.main_frame)
        title_frame.pack(fill=tk.X, pady=(0, 5))

        # Title and status
        title_text = f"📊 {self.spectrum_name}"
        if self.result_data.get('is_reference', False):
            title_text += " [REFERENCE]"

        ttk.Label(title_frame, text=title_text, font=('TkDefaultFont', 14, 'bold')).pack(side=tk.LEFT)

        status_icon = "✅" if self.result_data.get('status') == 'success' else "❌"
        status_text = f"{status_icon} {self.result_data.get('status', 'Unknown').title()}"
        ttk.Label(title_frame, text=status_text, font=('TkDefaultFont', 12)).pack(side=tk.RIGHT)

        # ENHANCED: Main content area with horizontal layout - spectrum on left, peak list on right
        content_frame = ttk.Frame(self.main_frame)
        content_frame.pack(fill=tk.BOTH, expand=True)

        # ENHANCED: Create horizontal paned window for spectrum plot (left) and peak analysis (right)
        main_paned = ttk.PanedWindow(content_frame, orient='horizontal')
        main_paned.pack(fill=tk.BOTH, expand=True)

        # Left panel: Spectrum plot and controls
        left_panel = ttk.Frame(main_paned)
        main_paned.add(left_panel, weight=3)  # 75% of space

        # Right panel: Peak analysis (moved to right side)
        right_panel = ttk.Frame(main_paned)
        main_paned.add(right_panel, weight=1)  # 25% of space

        # Left panel contents
        # ENHANCED: Spectrum plot
        plot_container = ttk.LabelFrame(left_panel, text="📈 Spectrum View", padding=5)
        plot_container.pack(fill=tk.BOTH, expand=True, pady=(0, 10))

        # Setup spectrum plot area
        self.setup_spectrum_plot(plot_container)

        # Setup contour controls below the plot in left panel
        contour_frame = ttk.LabelFrame(left_panel, text="🎨 Contour & Navigation Controls", padding=5)
        contour_frame.pack(fill=tk.X, pady=(0, 0))
        self.setup_contour_controls(contour_frame)

        # Right panel contents - Peak analysis (now on the right side next to spectrum)
        analysis_container = ttk.LabelFrame(right_panel, text="🔬 Peak Analysis", padding=5)
        analysis_container.pack(fill=tk.BOTH, expand=True)

        # Setup peak analysis panel on the right
        self.setup_peak_analysis_panel(analysis_container)

        # Bottom control panel
        control_frame = ttk.Frame(self.main_frame)
        control_frame.pack(fill=tk.X, pady=(5, 0))

        self.setup_controls(control_frame)

    def setup_voigt_analysis_tab(self):
        """Setup the Voigt analysis tab (integrated directly instead of popup)"""
        # ENHANCED: Add navigation controls similar to main GUI
        # Enhanced Peak Navigation Controls (like main GUI)
        nav_frame = ttk.LabelFrame(self.voigt_tab, text="🎯 Enhanced Peak Navigation", padding=10)
        nav_frame.pack(fill=tk.X, pady=(0, 10))

        # Peak selector with enhanced info
        peak_select_frame = ttk.Frame(nav_frame)
        peak_select_frame.pack(fill=tk.X, pady=2)

        # Peak number selector (will be set by peak selection)
        ttk.Label(peak_select_frame, text="Peak Number:").pack(side=tk.LEFT)
        self.voigt_peak_number = tk.IntVar(value=1)
        self.voigt_peak_spin = tk.Spinbox(peak_select_frame, from_=1, to=500, width=8,
                                         textvariable=self.voigt_peak_number,
                                         command=self.voigt_navigate_to_peak)
        self.voigt_peak_spin.pack(side=tk.LEFT, padx=5)

        # Peak info display
        self.voigt_peak_info_label = ttk.Label(peak_select_frame, text="Peak: -/-",
                                              font=('TkDefaultFont', 9, 'bold'))
        self.voigt_peak_info_label.pack(side=tk.LEFT, padx=20)

        # Enhanced navigation buttons
        nav_button_frame = ttk.Frame(nav_frame)
        nav_button_frame.pack(fill=tk.X, pady=5)

        ttk.Button(nav_button_frame, text="◀◀ Prev",
                  command=self.voigt_prev_peak, width=8).pack(side=tk.LEFT, padx=2)
        ttk.Button(nav_button_frame, text="🎯 Center",
                  command=self.voigt_center_on_peak, width=8).pack(side=tk.LEFT, padx=2)
        ttk.Button(nav_button_frame, text="Next ▶▶",
                  command=self.voigt_next_peak, width=8).pack(side=tk.LEFT, padx=2)
        ttk.Button(nav_button_frame, text="🔍 Zoom",
                  command=self.voigt_zoom_to_peak, width=8).pack(side=tk.LEFT, padx=2)
        ttk.Button(nav_button_frame, text="🔬 Analysis",
                  command=self.voigt_show_peak_analysis, width=10).pack(side=tk.LEFT, padx=2)

        # Create 2x2 grid for Voigt analysis plots like main GUI
        plot_container = ttk.Frame(self.voigt_tab)
        plot_container.pack(fill=tk.BOTH, expand=True)

        #self.fig_voigt, self.axes_voigt = plt.subplots(2, 2, figsize=(12, 8))
        self.fig_voigt, self.axes_voigt = plt.subplots(2, 1, figsize=(8, 6))

        self.canvas_voigt = FigureCanvasTkAgg(self.fig_voigt, plot_container)
        self.canvas_voigt.get_tk_widget().pack(fill=tk.BOTH, expand=True)

        # Add toolbar for navigation
        toolbar_voigt = NavigationToolbar2Tk(self.canvas_voigt, self.voigt_tab)
        toolbar_voigt.update()

        # Initialize with placeholder
        for ax in self.axes_voigt.flat:
            ax.text(0.5, 0.5, 'Select a peak to view Voigt analysis',
                   ha='center', va='center', transform=ax.transAxes, fontsize='small')
            ax.set_title('Voigt Analysis - No Peak Selected')

        self.canvas_voigt.draw()

    def populate_voigt_analysis_tab(self):
        """Populate the Voigt analysis tab with current peak data (reused from existing logic)"""
        try:
            # Get peak coordinates
            peak_x = self.selected_peak.get('Position_X', self.selected_peak.get('X_ppm', self.selected_peak.get('ppm_x')))
            peak_y = self.selected_peak.get('Position_Y', self.selected_peak.get('Y_ppm', self.selected_peak.get('ppm_y')))
            assignment = self.selected_peak.get('Assignment', 'Unknown')

            if peak_x is None or peak_y is None:
                # Show error message in plots
                for ax in self.axes_voigt.flat:
                    ax.clear()
                    ax.text(0.5, 0.5, "❌ Peak coordinates not found",
                           ha='center', va='center', transform=ax.transAxes, fontsize='small')
                self.canvas_voigt.draw()
                return

            # Clear previous plots
            for ax in self.axes_voigt.flat:
                ax.clear()

            # Populate with Voigt analysis using existing logic (same as popup window)
            if hasattr(self.integrator, 'nmr_data') and self.integrator.nmr_data is not None:
                data = self.integrator.nmr_data
                ppm_x = getattr(self.integrator, 'ppm_x_axis', None)
                ppm_y = getattr(self.integrator, 'ppm_y_axis', None)

                if ppm_x is not None and ppm_y is not None:
                    # Show zoomed spectrum around peak
                    zoom_x = 1.0  # ±1 ppm
                    zoom_y = 20.0  # ±20 ppm

                    # Calculate data indices for zoom region
                    x_indices = np.where((ppm_x >= peak_x - zoom_x) & (ppm_x <= peak_x + zoom_x))[0]
                    y_indices = np.where((ppm_y >= peak_y - zoom_y) & (ppm_y <= peak_y + zoom_y))[0]

                    if len(x_indices) > 0 and len(y_indices) > 0:
                        # Extract zoom data
                        zoom_data = data[y_indices[0]:y_indices[-1]+1, x_indices[0]:x_indices[-1]+1]
                        zoom_ppm_x = ppm_x[x_indices[0]:x_indices[-1]+1]
                        zoom_ppm_y = ppm_y[y_indices[0]:y_indices[-1]+1]

                        # CRITICAL FIX: Use stored results from series processing
                        assignment_key = self.selected_peak.get('Assignment', self.selected_peak.get('assignment', assignment))
                        stored_result = self.find_stored_voigt_result(assignment_key)
                        print(f"🔍 Using stored result for Voigt analysis tab: {stored_result is not None}")

                        # Use existing 1D analysis logic
                        self.populate_1d_voigt_analysis(self.axes_voigt, zoom_data, zoom_ppm_x, zoom_ppm_y,
                                                       peak_x, peak_y, assignment, stored_voigt_result=stored_result)

                        # Update figure title
                        self.fig_voigt.suptitle(f'Voigt Analysis: {assignment} at ({peak_x:.3f}, {peak_y:.1f}) ppm',
                                               fontsize='small', fontweight='bold')

            self.canvas_voigt.draw()

        except Exception as e:
            print(f"Error populating Voigt analysis tab: {e}")
            import traceback
            traceback.print_exc()

            # Show error in plots
            for ax in self.axes_voigt.flat:
                ax.clear()
                ax.text(0.5, 0.5, f"Error: {str(e)[:50]}",
                       ha='center', va='center', transform=ax.transAxes, fontsize='small')
            self.canvas_voigt.draw()

    def setup_spectrum_plot(self, parent):
        """Setup the spectrum plotting area with constrained height for scrollability"""
        if not MATPLOTLIB_AVAILABLE:
            ttk.Label(parent, text="Matplotlib not available for plotting",
                     foreground='red').pack(expand=True)
            return

        # ENHANCED: Calculate optimal figure size for small screen compatibility
        try:
            screen_height = self.window.winfo_screenheight()
            optimal_height = min(6, max(4, screen_height / 180))  # Adaptive height 4-6 inches
        except:
            optimal_height = 5  # Fallback

        # Create matplotlib figure with constrained size
        self.fig = Figure(figsize=(12, optimal_height), dpi=100, tight_layout=True)
        self.ax = self.fig.add_subplot(111)

        # Create canvas with constrained height for scrollability
        self.canvas = FigureCanvasTkAgg(self.fig, parent)
        canvas_widget = self.canvas.get_tk_widget()

        # ENHANCED: Constrain canvas height for better scrolling
        max_canvas_height = min(480, int(optimal_height * 100))  # Convert inches to pixels
        canvas_widget.configure(height=max_canvas_height)
        canvas_widget.pack(fill=tk.X, expand=False, padx=2, pady=2)  # Don't expand vertically

        # Add toolbar below the plot
        toolbar_frame = ttk.Frame(parent)
        toolbar_frame.pack(fill=tk.X, pady=(2, 0))

        toolbar = NavigationToolbar2Tk(self.canvas, toolbar_frame)
        toolbar.update()

        # Initial plot setup
        self.ax.set_title("Loading spectrum...", fontsize='small')
        self.ax.text(0.5, 0.5, "Spectrum will appear here",
                    transform=self.ax.transAxes, ha='center', va='center',
                    fontsize='small', color='gray')

        # Force initial draw
        self.canvas.draw()

        # Bind events for interaction
        self.canvas.mpl_connect('button_press_event', self.on_plot_click)
        self.canvas.mpl_connect('scroll_event', self.on_zoom)
        self.canvas.mpl_connect('motion_notify_event', self.on_mouse_move)

        # Store original limits for reset
        self.original_xlim = None
        self.original_ylim = None

    def setup_peak_analysis_panel(self, parent):
        """Setup the peak analysis panel"""
        # Peak info frame
        info_frame = ttk.LabelFrame(parent, text="🎯 Peak Information", padding=10)
        info_frame.pack(fill=tk.X, pady=(0, 5))

        self.peak_info_text = tk.Text(info_frame, height=8, width=40, wrap=tk.WORD,
                                     font=('Courier', 9))
        info_scrollbar = ttk.Scrollbar(info_frame, orient=tk.VERTICAL, command=self.peak_info_text.yview)
        self.peak_info_text.configure(yscrollcommand=info_scrollbar.set)

        self.peak_info_text.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        info_scrollbar.pack(side=tk.RIGHT, fill=tk.Y)

        # Peak list frame
        list_frame = ttk.LabelFrame(parent, text="📋 Peak List", padding=10)
        list_frame.pack(fill=tk.BOTH, expand=True, pady=(5, 0))

        # Create peak list treeview
        #peak_columns = ('Peak', 'Assignment', 'Height', 'SNR', 'Quality')
        #peak_columns = ('Peak', 'Assignment', 'Position_X', 'Position_Y', 'Height', 'SNR', 'Quality')
        peak_columns = ('Peak', 'Assignment', 'Position_X', 'Position_Y', 'Height', 'Volume', 'SNR', 'Quality', 'R²')
        self.peak_tree = ttk.Treeview(list_frame, columns=peak_columns, show='headings', height=12)

        column_widths = {
            'Peak': 50, 'Assignment': 100, 'Position_X': 80, 'Position_Y': 80,
            'Height': 80, 'Volume': 80, 'SNR': 60, 'Quality': 80, 'R²': 60
        }

        for col in peak_columns:
            self.peak_tree.heading(col, text=col)
            self.peak_tree.column(col, width=column_widths.get(col, 70))


        #for col in peak_columns:
        #    self.peak_tree.heading(col, text=col)
        #    self.peak_tree.column(col, width=70)

        peak_scrollbar = ttk.Scrollbar(list_frame, orient=tk.VERTICAL, command=self.peak_tree.yview)
        self.peak_tree.configure(yscrollcommand=peak_scrollbar.set)

        self.peak_tree.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        peak_scrollbar.pack(side=tk.RIGHT, fill=tk.Y)

###

        # Enhanced peak selection with coordinate display
        def on_peak_select_enhanced(event):
            """Enhanced peak selection showing coordinates"""
            selection = self.peak_tree.selection()
            if selection:
                item = self.peak_tree.item(selection[0])
                values = item['values']

                if len(values) >= 9:  # Updated for new column count
                    peak_num, assignment, pos_x, pos_y, height, volume, snr, quality, r_squared = values

                    # Update info text with selected peak details
                    info_text = f"Selected Peak Details:\n\n"
                    info_text += f"Peak Number: {peak_num}\n"
                    info_text += f"Assignment: {assignment}\n"
                    info_text += f"Position X: {pos_x} ppm\n"
                    info_text += f"Position Y: {pos_y} ppm\n"
                    info_text += f"Height: {height}\n"
                    info_text += f"SNR: {snr}\n"
                    info_text += f"Quality: {quality}\n\n"
                    info_text += "Click on another peak to view its details."
                    info_text += f"Volume: {volume}\n"
                    info_text += f"R-squared: {r_squared}\n"

                    self.peak_info_text.config(state=tk.NORMAL)
                    self.peak_info_text.delete(1.0, tk.END)
                    self.peak_info_text.insert(tk.END, info_text)
                    self.peak_info_text.config(state=tk.DISABLED)
                    # CRITICAL FIX: Trigger zoom functionality
                    self.on_peak_select(event)

        # Replace the basic handler with enhanced version

###
        # Bind peak selection
        self.peak_tree.bind('<Button-1>', self.on_peak_select)
        self.peak_tree.bind('<Button-1>', on_peak_select_enhanced)
        self.peak_tree.bind('<<TreeviewSelect>>', on_peak_select_enhanced)

        # Initialize with placeholder text
        self.peak_info_text.insert(tk.END, "Click on a peak in the spectrum or select from the list below to view detailed information.")
        self.peak_info_text.config(state=tk.DISABLED)


    def setup_controls(self, parent):
        """Setup the control panel"""
        # Navigation frame
        nav_frame = ttk.LabelFrame(parent, text="🧭 Navigation", padding=5)
        nav_frame.pack(side=tk.LEFT, fill=tk.Y, padx=(0, 10))

        ttk.Button(nav_frame, text="◀ Previous", command=self.previous_spectrum).pack(side=tk.LEFT, padx=(0, 5))
        ttk.Button(nav_frame, text="Next ▶", command=self.next_spectrum).pack(side=tk.LEFT, padx=(0, 5))
        ttk.Button(nav_frame, text="🏠 Center View", command=self.center_view).pack(side=tk.LEFT)

        # Analysis frame
        analysis_frame = ttk.LabelFrame(parent, text="🔬 Analysis", padding=5)
        analysis_frame.pack(side=tk.LEFT, fill=tk.Y, padx=(0, 10))

        ttk.Button(analysis_frame, text="📊 Voigt Analysis", command=self.switch_to_voigt_analysis).pack(side=tk.LEFT, padx=(0, 5))
        ttk.Button(analysis_frame, text="📈 Peak Summary", command=self.show_peak_summary).pack(side=tk.LEFT, padx=(0, 5))
        ttk.Button(analysis_frame, text="💾 Export Data", command=self.export_spectrum_data).pack(side=tk.LEFT)

        # Close button
        ttk.Button(parent, text="Close", command=self.window.destroy).pack(side=tk.RIGHT)

    def setup_contour_controls(self, parent):
        """Setup contour controls (imported from main window)"""
        # Contour controls frame
        contour_frame = ttk.LabelFrame(parent, text="🎨 Contour Settings", padding=5)
        contour_frame.pack(fill=tk.X, pady=(5, 0))

        # Main contour grid
        contour_grid = ttk.Frame(contour_frame)
        contour_grid.pack(fill=tk.X)

        # First row: Levels and Min Level
        ttk.Label(contour_grid, text="Levels:").grid(row=0, column=0, sticky=tk.W, padx=(0, 5))
        contour_levels_spin = tk.Spinbox(contour_grid, from_=3, to=50, increment=1, width=6,
                                        textvariable=self.contour_levels, command=self.update_spectrum_plot)
        contour_levels_spin.grid(row=0, column=1, sticky=tk.W, padx=(0, 15))
        contour_levels_spin.bind('<KeyRelease>', lambda e: self.update_spectrum_plot())
        contour_levels_spin.bind('<FocusOut>', lambda e: self.update_spectrum_plot())

        ttk.Label(contour_grid, text="Min Level:").grid(row=0, column=2, sticky=tk.W, padx=(0, 5))
        contour_min_spin = tk.Spinbox(contour_grid, from_=0.001, to=1.0, increment=0.001, width=8,
                                     textvariable=self.contour_min, format="%.3f", command=self.update_spectrum_plot)
        contour_min_spin.grid(row=0, column=3, sticky=tk.W, padx=(0, 15))
        contour_min_spin.bind('<KeyRelease>', lambda e: self.update_spectrum_plot())
        contour_min_spin.bind('<FocusOut>', lambda e: self.update_spectrum_plot())

        # Second row: Increment and Colormap
        ttk.Label(contour_grid, text="Increment:").grid(row=1, column=0, sticky=tk.W, padx=(0, 5))
        contour_inc_spin = tk.Spinbox(contour_grid, from_=0.01, to=10.0, increment=0.01, width=6,
                                     textvariable=self.contour_increment, format="%.2f", command=self.update_spectrum_plot)
        contour_inc_spin.grid(row=1, column=1, sticky=tk.W, padx=(0, 15))
        contour_inc_spin.bind('<KeyRelease>', lambda e: self.update_spectrum_plot())
        contour_inc_spin.bind('<FocusOut>', lambda e: self.update_spectrum_plot())

        ttk.Label(contour_grid, text="Colormap:").grid(row=1, column=2, sticky=tk.W, padx=(0, 5))
        colormap_combo = ttk.Combobox(contour_grid, textvariable=self.colormap_var, width=10, state='readonly')
        colormap_combo['values'] = ('viridis', 'plasma', 'inferno', 'magma', 'cividis', 'turbo',
                                   'hot', 'cool', 'spring', 'summer', 'autumn', 'winter')
        colormap_combo.grid(row=1, column=3, sticky=tk.W)
        colormap_combo.bind('<<ComboboxSelected>>', lambda e: self.update_spectrum_plot())

        # Quick presets
        preset_frame = ttk.Frame(contour_frame)
        preset_frame.pack(fill=tk.X, pady=(5, 0))

        ttk.Label(preset_frame, text="Quick Presets:").pack(side=tk.LEFT)
        ttk.Button(preset_frame, text="Low Detail", width=8, command=self.preset_low).pack(side=tk.LEFT, padx=2)
        ttk.Button(preset_frame, text="Medium", width=8, command=self.preset_medium).pack(side=tk.LEFT, padx=2)
        ttk.Button(preset_frame, text="High Detail", width=8, command=self.preset_high).pack(side=tk.LEFT, padx=2)
        ttk.Button(preset_frame, text="Auto", width=8, command=self.preset_auto).pack(side=tk.LEFT, padx=2)

    def load_spectrum_data_proper(self):
        """Load spectrum data properly like main window with performance optimization"""
        try:
            if not CORE_MODULES_AVAILABLE:
                raise Exception("Core NMR modules not available")

            # Show loading message
            if hasattr(self, 'ax'):
                self.ax.clear()
                self.ax.text(0.5, 0.5, "⏳ Loading spectrum data...",
                            transform=self.ax.transAxes, ha='center', va='center',
                            fontsize='small', bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.8))
                self.ax.set_title(f"Loading {self.spectrum_name}...", fontsize='small')
                self.canvas.draw_idle()
                self.window.update_idletasks()  # Keep responsive

            print(f"🔍 Loading spectrum data for: {self.spectrum_name}")

            # Create integrator instance with MINIMAL copying
            if (hasattr(self.series_processor, 'main_gui_integrator') and
                self.series_processor.main_gui_integrator is not None):
                IntegratorClass = type(self.series_processor.main_gui_integrator)
                self.integrator = IntegratorClass()

                # Copy essential settings AND enhanced fitter for local quality assessment
                main_integrator = self.series_processor.main_gui_integrator
                essential_attrs = ['processing_mode', 'threshold_multiplier', 'enhanced_fitter', 'fitting_parameters']
                for attr_name in essential_attrs:
                    if hasattr(main_integrator, attr_name):
                        try:
                            setattr(self.integrator, attr_name, getattr(main_integrator, attr_name))
                        except:
                            pass
            else:
                self.integrator = EnhancedVoigtIntegrator()

            # Check file exists
            if not self.spectrum_file_path or not os.path.exists(self.spectrum_file_path):
                raise Exception(f"Spectrum file not found: {self.spectrum_file_path}")

            # Load using nmrglue (same as main window)
            success = self.integrator.load_nmr_file(self.spectrum_file_path)
            if not success:
                raise Exception(f"Failed to load NMR file")

            # Force data loading if needed
            if (hasattr(self.integrator, 'nmr_data') and self.integrator.nmr_data is None and
                hasattr(self.integrator, '_load_nmr_data_only')):
                success = self.integrator._load_nmr_data_only(self.spectrum_file_path)
                if not success:
                    raise Exception(f"Failed to load spectrum data")

            # Extract data from nmr_dict if needed
            if hasattr(self.integrator, 'nmr_dict') and self.integrator.nmr_dict:
                for key, value in self.integrator.nmr_dict.items():
                    if hasattr(value, 'shape') and len(getattr(value, 'shape', [])) == 2:
                        self.integrator.nmr_data = value
                        break

            # Verify data loaded
            if not hasattr(self.integrator, 'nmr_data') or self.integrator.nmr_data is None:
                raise Exception("NMR data is None after loading")

            print(f"✅ NMR data loaded, shape: {getattr(self.integrator.nmr_data, 'shape', 'unknown')}")

            # Set peak data from results - FIXED to use integration_results for both
            integration_results = self.result_data.get('integration_results', [])
            fitted_peaks = self.result_data.get('fitted_peaks', [])

            # Use integration_results as primary source (has proper format)
            if integration_results:
                # Set both fitted_peaks and integration_results for visualization compatibility
                self.integrator.fitted_peaks = integration_results.copy()
                self.integrator.integration_results = integration_results.copy()
                print(f"✅ Set peak data from integration_results: {len(integration_results)} peaks")

                # Debug: Check first peak structure
                if integration_results:
                    first_peak = integration_results[0]
                    print(f"   📊 First peak fields: {list(first_peak.keys())}")
                    print(f"   📊 Sample values: height={first_peak.get('height', 'MISSING')}, snr={first_peak.get('snr', 'MISSING')}")

            elif fitted_peaks:
                self.integrator.fitted_peaks = fitted_peaks.copy()
                self.integrator.integration_results = fitted_peaks.copy()
                print(f"✅ Set peak data from fitted_peaks fallback: {len(fitted_peaks)} peaks")
            else:
                print("⚠️ No peak data found in result_data")
                self.integrator.fitted_peaks = []
                self.integrator.integration_results = []

            # Plot spectrum using SpectrumPlotter (same as main window)
            self.plot_spectrum_like_main_window()

            # Populate peak list properly
            self.populate_peak_list()

            print(f"✅ Spectrum browser loaded successfully!")

            # CRITICAL: Ensure fitted_peaks has proper format for visualization
            if hasattr(self.integrator, 'fitted_peaks') and self.integrator.fitted_peaks:
                # Debug peak data format
                first_peak = self.integrator.fitted_peaks[0]
                print(f"   🔍 Peak data format for visualization:")
                print(f"      assignment: {first_peak.get('assignment', 'MISSING')}")
                print(f"      ppm_x: {first_peak.get('ppm_x', 'MISSING')}")
                print(f"      ppm_y: {first_peak.get('ppm_y', 'MISSING')}")
                print(f"      detected: {first_peak.get('detected', 'MISSING')}")

        except Exception as e:
            print(f"❌ Error loading spectrum: {e}")
            self.show_error(f"Error loading spectrum:\n{str(e)}")
            # Still try to populate peak list with available data
            self.populate_peak_list_simple()

    def plot_spectrum_like_main_window(self):
        """Plot spectrum using same logic as main window"""
        try:
            self.ax.clear()

            if not self.integrator or not hasattr(self.integrator, 'nmr_data') or self.integrator.nmr_data is None:
                self.show_error("No spectrum data available")
                return

            # Use SpectrumPlotter exactly like main window
            plotter = SpectrumPlotter(self.fig, self.ax)
            plotter.plot_spectrum(self.integrator)

            # Plot peaks with overlays
            if hasattr(self.integrator, 'fitted_peaks') and self.integrator.fitted_peaks:
                plotter.plot_peaks(self.integrator, show_detected=True, show_assigned=True)

            # Add title with status
            detection_rate = self.result_data.get('detection_rate', 0)
            status_icon = "✅" if self.result_data.get('status') == 'success' else "❌"
            title = f"{status_icon} {self.spectrum_name} | {self.result_data.get('detected_peaks', 0)}/{self.result_data.get('total_peaks', 0)} peaks ({detection_rate:.1f}%)"
            if self.result_data.get('is_reference', False):
                title += " [REFERENCE]"
            self.ax.set_title(title, fontsize='small', fontweight='bold')

            # Store original limits for reset
            if self.original_xlim is None:
                self.original_xlim = self.ax.get_xlim()
                self.original_ylim = self.ax.get_ylim()

            self.canvas.draw_idle()

        except Exception as e:
            self.show_error(f"Plotting error: {str(e)}")

    def show_error(self, error_msg):
        """Show error message"""
        if hasattr(self, 'ax'):
            self.ax.clear()
            self.ax.text(0.5, 0.5, f"❌ {error_msg}",
                        transform=self.ax.transAxes, ha='center', va='center',
                        fontsize='small', color='red')
            self.ax.set_title("Error", fontsize='small')
            self.canvas.draw_idle()

    def load_spectrum_data_async_old(self):
        """Load spectrum data asynchronously to keep UI responsive"""
        try:
            # Show loading message immediately
            if hasattr(self, 'ax'):
                self.ax.clear()
                self.ax.text(0.5, 0.5, "🔄 Loading spectrum data...\n\nButtons are active while loading",
                            transform=self.ax.transAxes, ha='center', va='center',
                            fontsize='small', color='blue')
                self.ax.set_title(f"Loading {self.spectrum_name}...", fontsize='small')
                self.canvas.draw_idle()

            # Schedule actual loading to run after UI is ready
            self.window.after(50, self.do_simple_load)

        except Exception as e:
            self.show_simple_error(f"Failed to start loading: {str(e)}")

    def do_simple_load(self):
        """Simplified, fast loading without heavy processing"""
        try:
            if not CORE_MODULES_AVAILABLE:
                self.show_simple_error("Core NMR modules not available")
                return

            print(f"🔄 Fast loading: {self.spectrum_name}")

            # Quick integrator setup - copy enhanced fitter for local quality assessment
            if (hasattr(self.series_processor, 'main_gui_integrator') and
                self.series_processor.main_gui_integrator is not None):
                IntegratorClass = type(self.series_processor.main_gui_integrator)
                self.integrator = IntegratorClass()

                # Copy essential settings including enhanced fitter
                main_integrator = self.series_processor.main_gui_integrator
                essential_attrs = ['enhanced_fitter', 'fitting_parameters']
                for attr_name in essential_attrs:
                    if hasattr(main_integrator, attr_name):
                        try:
                            setattr(self.integrator, attr_name, getattr(main_integrator, attr_name))
                        except:
                            pass
            else:
                self.integrator = EnhancedVoigtIntegrator()

            # Simple file existence check
            if not self.spectrum_file_path or not os.path.exists(self.spectrum_file_path):
                self.show_simple_error(f"File not found: {os.path.basename(self.spectrum_file_path)}")
                return

            # Try to load data (this should work since main GUI works)
            success = self.integrator.load_nmr_file(self.spectrum_file_path)
            if not success:
                self.show_simple_error("Failed to load NMR file")
                return

            # Force data loading if needed
            if (hasattr(self.integrator, 'nmr_data') and self.integrator.nmr_data is None and
                hasattr(self.integrator, '_load_nmr_data_only')):
                success = self.integrator._load_nmr_data_only(self.spectrum_file_path)
                if not success:
                    self.show_simple_error("Failed to load spectrum data")
                    return

            # Quick data extraction
            if hasattr(self.integrator, 'nmr_dict') and self.integrator.nmr_dict:
                for key, value in self.integrator.nmr_dict.items():
                    if hasattr(value, 'shape') and len(getattr(value, 'shape', [])) == 2:
                        self.integrator.nmr_data = value
                        break

            # Set results data
            if self.result_data.get('fitted_peaks'):
                self.integrator.fitted_peaks = self.result_data['fitted_peaks']
            if self.result_data.get('integration_results'):
                self.integrator.integration_results = self.result_data['integration_results']

            # Simple initial plot
            self.show_simple_spectrum()

            # Populate peak list without heavy processing
            self.populate_peak_list_simple()

            print(f"✅ Fast loading complete!")

        except Exception as e:
            self.show_simple_error(f"Loading error: {str(e)}")

    def show_simple_spectrum(self):
        """Show spectrum with minimal processing for responsiveness"""
        try:
            if not hasattr(self.integrator, 'nmr_data') or self.integrator.nmr_data is None:
                self.show_simple_error("No spectrum data available")
                return

            self.ax.clear()

            # Fast display using imshow
            data = np.abs(self.integrator.nmr_data)
            self.ax.imshow(data, aspect='auto', cmap='viridis', origin='lower')

            # Add basic info
            detection_rate = self.result_data.get('detection_rate', 0)
            status_icon = "✅" if self.result_data.get('status') == 'success' else "❌"
            title = f"{status_icon} {self.spectrum_name} ({detection_rate:.1f}%)"
            self.ax.set_title(title, fontsize='small')
            self.ax.set_xlabel('F2 (points)')
            self.ax.set_ylabel('F1 (points)')

            # Store original limits
            self.original_xlim = self.ax.get_xlim()
            self.original_ylim = self.ax.get_ylim()

            self.canvas.draw_idle()

        except Exception as e:
            self.show_simple_error(f"Plot error: {str(e)}")

    def populate_peak_list_simple(self):
        """Enhanced fallback peak list population with coordinate extraction"""
        try:
            # Clear existing items
            for item in self.peak_tree.get_children():
                self.peak_tree.delete(item)

            integration_results = self.result_data.get('integration_results', [])

            if integration_results:
                first_peak = integration_results[0]
                print(f"   First peak keys: {list(first_peak.keys())}")
                print(f"   First peak Position_X: {first_peak.get('Position_X', 'MISSING')}")
                print(f"   First peak Position_Y: {first_peak.get('Position_Y', 'MISSING')}")
                print(f"   First peak Assignment: {first_peak.get('Assignment', 'MISSING')}")
            else:
                print(f"   ⚠️ No integration_results found!")
                # Check for alternative data sources
                fitted_results = self.result_data.get('fitted_results', [])
                fitted_peaks = self.result_data.get('fitted_peaks', [])
                print(f"   Alternative sources: fitted_results={len(fitted_results)}, fitted_peaks={len(fitted_peaks)}")

            fitted_peaks = self.result_data.get('fitted_peaks', [])
            fitted_results = self.result_data.get('fitted_results', [])

            # Use whichever has data, prioritizing integration_results
            peak_data = integration_results or fitted_results or fitted_peaks

              # Determine data source for logging
            if integration_results:
                  data_source = 'integration_results'
            elif fitted_results:
                  data_source = 'fitted_results'
            elif fitted_peaks:
                  data_source = 'fitted_peaks'
            else:
                  data_source = 'no_source'

            print(f"📊 Simple fallback using {len(peak_data)} peaks from {data_source}")


            for i, peak in enumerate(peak_data):
                # STANDARDIZED assignment extraction (primary first)
                assignment = str(peak.get('assignment') or
                               peak.get('Assignment') or
                               peak.get('name') or
                               f'Peak_{i+1}')

               # STANDARDIZED coordinate extraction (primary fields first)
                pos_x = (peak.get('ppm_x') or
                         peak.get('Position_X') or
                         peak.get('x') or 0)
                pos_y = (peak.get('ppm_y') or
                         peak.get('Position_Y') or
                         peak.get('y') or 0)

                # Format coordinates with proper error handling
                try:
                    pos_x_val = float(pos_x) if pos_x != 0 else 0
                    pos_y_val = float(pos_y) if pos_y != 0 else 0
                    pos_x_str = f"{pos_x_val:.3f}" if pos_x_val != 0 else "N/A"
                    pos_y_str = f"{pos_y_val:.1f}" if pos_y_val != 0 else "N/A"
                except (ValueError, TypeError):
                    pos_x_str = "Error"
                    pos_y_str = "Error"

                # STANDARDIZED property extraction (primary fields first)
                height = peak.get('height') or peak.get('Height') or 'N/A'
                volume = peak.get('volume') or peak.get('Volume') or 'N/A'
                snr = peak.get('snr') or peak.get('SNR') or 'N/A'
                quality = peak.get('quality') or peak.get('Quality') or peak.get('fitting_quality') or 'Unknown'
                r_squared = peak.get('r_squared') or peak.get('R_Squared') or peak.get('avg_r_squared') or 'N/A'

                self.peak_tree.insert('', 'end',
                                    values=(f"{i+1}", assignment, pos_x_str, pos_y_str,
                                          str(height)[:8], str(snr)[:6], str(quality)[:10]))

            print(f"✅ Enhanced simple peak list populated with {len(peak_data)} peaks")

        except Exception as e:
            print(f"❌ Error in enhanced simple peak list population: {e}")
            # Add diagnostic entry showing what data is available
            available_keys = list(self.result_data.keys()) if self.result_data else []
            self.peak_tree.insert('', 'end',
                                values=("1", "Data Error", "N/A", "N/A", "N/A", "N/A",
                                      f"Keys: {', '.join(available_keys[:3])}" if available_keys else "No data"))

    def populate_peak_list(self):
        """Populate the peak list with integration results (full version)"""
        try:
            # Clear existing items
            for item in self.peak_tree.get_children():
                self.peak_tree.delete(item)

            integration_results = self.result_data.get('integration_results', [])

            for i, result in enumerate(integration_results):
                try:
                   # STANDARDIZED assignment extraction (primary first)
                    assignment = str(result.get('assignment') or
                                   result.get('Assignment') or
                                   f'Peak_{i+1}')

                    # Height extraction (primary → legacy → default)
                    height_raw = result.get('height') or result.get('Height') or 0
                    if hasattr(height_raw, '__iter__') and not isinstance(height_raw, str):
                        height_val = float(height_raw[0]) if len(height_raw) > 0 else 0
                    else:
                        height_val = float(height_raw) if height_raw is not None else 0
                    height = f"{height_val:.2e}" if height_val > 0 else "N/A"

                    # SNR extraction (primary → legacy → default)
                    snr_raw = result.get('snr') or result.get('SNR') or 0
                    if hasattr(snr_raw, '__iter__') and not isinstance(snr_raw, str):
                        snr_val = float(snr_raw[0]) if len(snr_raw) > 0 else 0
                    else:
                        snr_val = float(snr_raw) if snr_raw is not None else 0
                    snr = f"{snr_val:.1f}" if snr_val > 0 else "N/A"

                    # Quality extraction (primary → legacy → default)
                    quality = str(result.get('quality') or
                                result.get('Quality') or
                                result.get('fitting_quality') or
                                'Unknown')[:10]

                    # Volume extraction (add to display - not currently shown)
                    volume_raw = result.get('volume') or result.get('Volume') or 0
                    volume_val = float(volume_raw) if volume_raw else 0

                    # R-squared extraction
                    r_squared_raw = result.get('r_squared') or result.get('R_Squared') or result.get('avg_r_squared') or 0
                    r_squared_val = float(r_squared_raw) if r_squared_raw else 0

                    print(f"   📊 Peak {i+1} display: height={height}, snr={snr}, quality={quality}")

                    # Add color coding based on quality
                    tags = []
                    if snr_val >= 10:
                        tags.append('good')
                    elif snr_val >= 5:
                        tags.append('fair')
                    elif snr_val > 0:
                        tags.append('poor')


                    # Use the SAME logic as the working peak tree population
                    pos_x = (result.get('Position_X') or
                             result.get('ppm_x') or
                             result.get('x') or 0)
                    pos_y = (result.get('Position_Y') or
                             result.get('ppm_y') or
                             result.get('y') or 0)

                    #  With proper float conversion and validation
                    try:
                        pos_x_val = float(pos_x) if pos_x != 0 else 0
                        pos_y_val = float(pos_y) if pos_y != 0 else 0
                    except (ValueError, TypeError):
                        pos_x_val = pos_y_val = 0


                    # Format coordinates for display with validation
                    try:
                        pos_x_val = float(pos_x) if pos_x != 0 else 0
                        pos_y_val = float(pos_y) if pos_y != 0 else 0
                        pos_x_str = f"{pos_x_val:.3f}" if pos_x_val != 0 else "N/A"
                        pos_y_str = f"{pos_y_val:.1f}" if pos_y_val != 0 else "N/A"
                    except (ValueError, TypeError):
                        pos_x_str = "Error"
                        pos_y_str = "Error"

                    # Format additional properties for display
                    volume_str = f"{volume_val:.1e}" if volume_val > 0 else "N/A"
                    r_squared_str = f"{r_squared_val:.3f}" if r_squared_val > 0 else "N/A"

                    item_id = self.peak_tree.insert('', 'end',
                                                  values=(f"{i+1}", assignment[:12], pos_x_str, pos_y_str, height, volume_str, snr, quality, r_squared_str),
                                                  tags=tags)

                    # ENHANCED DEBUG: Log coordinate extraction for first few peaks
                    if i < 3:
                        print(f"   📊 Peak {i+1} tree display (ENHANCED):")
                        print(f"      Assignment: {assignment}")
                        print(f"      Available coord fields: Position_X={result.get('Position_X')}, ppm_x={result.get('ppm_x')}, x={result.get('x')}")
                        print(f"      Selected coordinates: X={pos_x} → {pos_x_str}, Y={pos_y} → {pos_y_str}")
                        print(f"      All result keys: {list(result.keys())}")


                except Exception as e:
                    print(f"Error processing peak {i}: {e}")
                    # Add basic entry without color coding
                    self.peak_tree.insert('', 'end',
                                        values=(f"{i+1}", f"Peak_{i+1}", "0.0e+00", "0.0", "Unknown"))

            # Configure tags for color coding
            self.peak_tree.tag_configure('good', background='#f1f8e9')
            self.peak_tree.tag_configure('fair', background='#fffde7')
            self.peak_tree.tag_configure('poor', background='#fff3e0')


            # ENHANCED INFO TEXT with comprehensive coordinate analysis
            info_text = "Click on peaks in the plot or select from the list below for detailed information.\n\n"
            info_text += f"Spectrum: {self.spectrum_name}\n"
            info_text += f"Total Peaks: {len(integration_results)}\n"
            info_text += f"Detection Rate: {self.result_data.get('detection_rate', 0):.1f}%\n\n"

            # COMPREHENSIVE coordinate analysis with multiple field checking
            if integration_results:
                # Check multiple coordinate field possibilities
                x_coords = []
                y_coords = []
                coord_sources = []

                for r in integration_results:
                    # Try multiple coordinate sources (using working logic from peak tree)
                    pos_x = (r.get('Position_X') or
                            r.get('ppm_x') or
                            r.get('x') or 0)
                    pos_y = (r.get('Position_Y') or
                            r.get('ppm_y') or
                            r.get('y') or 0)

                    # Validate and convert to float (same as peak tree logic)
                    try:
                         x_val = float(pos_x) if pos_x != 0 else 0
                         y_val = float(pos_y) if pos_y != 0 else 0
                    except (ValueError, TypeError):
                         x_val = y_val = 0

                    if x_val != 0 and y_val != 0:
                    #if x_val != 0 and y_val != 0:
                        x_coords.append(float(x_val))
                        y_coords.append(float(y_val))

                        # Identify coordinate source for debugging
                        if 'Position_X' in r:
                            coord_sources.append('Position_X/Y')
                        elif 'ppm_x' in r:
                            coord_sources.append('ppm_x/y')
                        elif 'x' in r:
                            coord_sources.append('x/y')

                if x_coords and y_coords:
                    info_text += f"📍 Coordinate Ranges ({len(x_coords)} peaks with coordinates):\n"
                    info_text += f"  X: {min(x_coords):.3f} to {max(x_coords):.3f} ppm\n"
                    info_text += f"  Y: {min(y_coords):.1f} to {max(y_coords):.1f} ppm\n"

                    # Show coordinate sources for debugging
                    unique_sources = list(set(coord_sources))
                    info_text += f"  Sources: {', '.join(unique_sources)}\n"
                else:
                    info_text += "⚠️ Warning: No coordinate data found in peak results\n"

                    # DEBUG: Show what fields ARE available
                    if integration_results:
                        sample_fields = list(integration_results[0].keys())
                        info_text += f"  Available fields: {', '.join(sample_fields[:8])}\n"
                        if len(sample_fields) > 8:
                            info_text += f"  ... and {len(sample_fields)-8} more\n"
            else:
                info_text += "⚠️ No peak data available\n"

            self.peak_info_text.config(state=tk.NORMAL)
            self.peak_info_text.delete(1.0, tk.END)
            self.peak_info_text.insert(tk.END, info_text)
            self.peak_info_text.config(state=tk.DISABLED)


        except Exception as e:
            print(f"Error populating peak list: {e}")
            # Fall back to simple version
            try:
                self.populate_peak_list_simple()
            except:
                pass

    def show_simple_error(self, error_msg):
        """Show error with minimal overhead"""
        print(f"❌ {error_msg}")

        if hasattr(self, 'ax'):
            self.ax.clear()
            self.ax.text(0.5, 0.5, f"❌ {error_msg}",
                        transform=self.ax.transAxes, ha='center', va='center',
                        fontsize='small', color='red')
            self.ax.set_title("Error", fontsize='small')
            self.canvas.draw_idle()

    def plot_spectrum(self):
        """Optimized plotting with responsive UI"""
        if not self.integrator or not hasattr(self.integrator, 'nmr_data') or self.integrator.nmr_data is None:
            self.show_simple_error("No data to plot")
            return

        try:
            # Show progress during plotting
            self.ax.clear()
            self.ax.text(0.5, 0.5, "📊 Rendering spectrum...",
                        transform=self.ax.transAxes, ha='center', va='center',
                        fontsize='small', bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8))
            self.canvas.draw_idle()
            self.window.update_idletasks()  # Keep UI responsive

            # Get settings (with sanity checks) - expanded range for user-requested defaults
            levels = max(3, min(self.contour_levels.get(), 50))  # Support up to 50 levels
            min_level = max(0.01, min(self.contour_min.get(), 0.8))  # Support down to 0.01
            colormap = self.colormap_var.get()

            data = self.integrator.nmr_data
            ppm_x = getattr(self.integrator, 'ppm_x_axis', None)
            ppm_y = getattr(self.integrator, 'ppm_y_axis', None)

            if ppm_x is not None and ppm_y is not None:
                # Clear again for actual plot
                self.ax.clear()

                # Optimized contour calculation
                data_abs = np.abs(data)
                max_intensity = np.max(data_abs)
                min_intensity = min_level * max_intensity

                # Use requested levels but with reasonable performance cap
                effective_levels = min(levels, 40)  # Cap at 40 to support user-requested default
                level_values = np.geomspace(min_intensity, max_intensity * 0.8, effective_levels)

                # Thinner lines for faster rendering
                self.ax.contour(ppm_x, ppm_y, data, levels=level_values, cmap=colormap, linewidths=0.3)
                self.ax.set_xlabel('¹H (ppm)', fontsize='small')
                self.ax.set_ylabel('¹⁵N (ppm)', fontsize='small')
                self.ax.invert_xaxis()
                self.ax.invert_yaxis()  # Fix: Also invert Y-axis for proper NMR orientation

                # Fix: Add peak markers after contour plot
                self.plot_peak_markers()
            else:
                # Fallback
                self.ax.imshow(np.abs(data), aspect='auto', cmap=colormap, origin='lower')
                self.ax.set_xlabel('F2 (points)')
                self.ax.set_ylabel('F1 (points)')

            # Add title
            detection_rate = self.result_data.get('detection_rate', 0)
            status_icon = "✅" if self.result_data.get('status') == 'success' else "❌"
            title = f"{status_icon} {self.spectrum_name} ({detection_rate:.1f}%)"
            self.ax.set_title(title, fontsize='small')

            self.canvas.draw_idle()

        except Exception as e:
            self.show_simple_error(f"Plot error: {str(e)}")

    def plot_peak_markers(self):
        """Plot peak markers on the spectrum - called after contour plotting"""
        try:
            print(f"🔍 DEBUG: Attempting to plot peak markers...")

            # ENHANCED: First try to use integration_results (always available from our fixes)
            if self.result_data.get('integration_results'):
                integration_results = self.result_data['integration_results']
                x_coords = []
                y_coords = []
                labels = []
                colors = []

                print(f"🔍 DEBUG: Found {len(integration_results)} integration results")

                for result in integration_results:
                    # Extract coordinates using the standard integration format
                    x_ppm = result.get('Position_X')
                    y_ppm = result.get('Position_Y')
                    assignment = result.get('Assignment', '')
                    detected = result.get('Detected', True)
                    quality = result.get('Quality', 'Unknown')

                    print(f"🔍 DEBUG: Result {assignment}: X={x_ppm}, Y={y_ppm}, detected={detected}")

                    if x_ppm is not None and y_ppm is not None and x_ppm != 0 and y_ppm != 0:
                        x_coords.append(float(x_ppm))
                        y_coords.append(float(y_ppm))
                        labels.append(assignment)

                        # Color code by quality
                        if quality == 'Excellent':
                            colors.append('lime')
                        elif quality == 'Good':
                            colors.append('green')
                        elif quality == 'Fair':
                            colors.append('orange')
                        else:
                            colors.append('red')

                # Plot integration results as colored circles
                if x_coords and y_coords:
                    self.ax.scatter(x_coords, y_coords, c=colors, s=40, marker='o',
                                  alpha=0.8, edgecolors='black', linewidths=1,
                                  label='Peak Positions', zorder=10)

                    # Add peak labels
                    for x, y, label in zip(x_coords, y_coords, labels):
                        if label and label != 'Unknown':
                            self.ax.annotate(label, (x, y), xytext=(5, 5),
                                           textcoords='offset points', fontsize='small',
                                           color='black', weight='bold',
                                           bbox=dict(boxstyle='round,pad=0.2', facecolor='white', alpha=0.8))

                print(f"✅ Plotted {len(x_coords)} peak markers from integration results")

            # FALLBACK: Try fitted_peaks with multiple field name variations
            elif hasattr(self.integrator, 'fitted_peaks') and self.integrator.fitted_peaks:
                x_coords = []
                y_coords = []
                labels = []

                print(f"🔍 DEBUG: Fallback to fitted_peaks: {len(self.integrator.fitted_peaks)} peaks")

                for peak in self.integrator.fitted_peaks:
                    if isinstance(peak, dict):
                        # Try multiple field name variations
                        x_ppm = (peak.get('ppm_x') or peak.get('x_ppm') or
                                peak.get('Position_X') or peak.get('center_x', 0))
                        y_ppm = (peak.get('ppm_y') or peak.get('y_ppm') or
                                peak.get('Position_Y') or peak.get('center_y', 0))
                        assignment = peak.get('assignment', peak.get('Assignment', ''))
                        detected = peak.get('detected', True)

                        print(f"🔍 DEBUG: Peak {assignment}: X={x_ppm}, Y={y_ppm}")

                        if x_ppm and y_ppm and x_ppm != 0 and y_ppm != 0:
                            x_coords.append(float(x_ppm))
                            y_coords.append(float(y_ppm))
                            labels.append(assignment)

                # Plot fitted peaks as blue circles
                if x_coords and y_coords:
                    self.ax.scatter(x_coords, y_coords, c='cyan', s=30, marker='s',
                                  alpha=0.7, edgecolors='blue', linewidths=1,
                                  label='Fitted Peaks', zorder=10)

                    # Add peak labels
                    for x, y, label in zip(x_coords, y_coords, labels):
                        if label:
                            self.ax.annotate(label, (x, y), xytext=(3, 3),
                                           textcoords='offset points', fontsize='small',
                                           color='blue', weight='bold')

                print(f"✅ Plotted {len(x_coords)} fitted peak markers")
            else:
                print("⚠️ No peak data available for plotting")

        except Exception as e:
            print(f"❌ Error plotting peak markers: {e}")
            import traceback
            traceback.print_exc()

    def on_plot_click(self, event):
        """Handle clicks on the spectrum plot to select peaks"""
        if event.inaxes != self.ax:
            return

        # Find the closest peak to the click
        click_x, click_y = event.xdata, event.ydata
        if click_x is None or click_y is None:
            return

        integration_results = self.result_data.get('integration_results', [])
        if not integration_results:
            return

        # Find closest peak
        min_distance = float('inf')
        closest_peak = None

        for peak in integration_results:
            peak_x = peak.get('Position_X', 0)
            peak_y = peak.get('Position_Y', 0)

            # Calculate distance (weighted more heavily on X-axis for NMR)
            distance = ((click_x - peak_x) * 2) ** 2 + (click_y - peak_y) ** 2

            if distance < min_distance:
                min_distance = distance
                closest_peak = peak

        # Select the closest peak if it's reasonably close
        if closest_peak and min_distance < 1.0:  # Adjust threshold as needed
            self.select_peak(closest_peak)

    def on_peak_select(self, event):
        """Handle peak selection from the peak list"""
        selection = self.peak_tree.selection()
        if not selection:
            return

        item = selection[0]
        peak_index = int(self.peak_tree.item(item, 'values')[0]) - 1  # Convert to 0-based index

        integration_results = self.result_data.get('integration_results', [])
        if 0 <= peak_index < len(integration_results):
            self.select_peak(integration_results[peak_index])




    def select_peak(self, peak_data):
        """Select and display information for a specific peak"""
        self.selected_peak = peak_data

        # Update peak information display
        self.update_peak_info(peak_data)

        # Highlight peak in the plot if possible
        self.highlight_peak_in_plot(peak_data)

        # Zoom to peak location
        self.zoom_to_peak(peak_data)

        # TABBED INTERFACE: Auto-update Voigt analysis tab if currently selected
        try:
            current_tab = self.notebook.select()
            voigt_tab_id = str(self.voigt_tab)
            if current_tab == voigt_tab_id:
                # If Voigt analysis tab is currently active, update it with the new peak
                self.populate_voigt_analysis_tab()
        except Exception:
            # Gracefully handle any tab selection errors
            pass

    def find_stored_voigt_result(self, assignment):
        """
        Find stored Voigt fitting results for a peak by assignment
        Returns the detailed fitted_results data from series processing
        """
        if not assignment:
            return None

        # First check if we have fitted_results from series processing
        fitted_results = self.result_data.get('fitted_results', [])
        if fitted_results:
            for result in fitted_results:
                # STANDARDIZED assignment and R-squared access (primary first)
                result_assignment = result.get('assignment') or result.get('Assignment')
                if result_assignment == assignment:
                    r_squared = result.get('r_squared') or result.get('R_Squared') or result.get('avg_r_squared') or 0
                    print(f"✅ Found stored Voigt result for {assignment}: R²={r_squared:.3f}")

                #if result.get('assignment') == assignment:
                #    print(f"✅ Found stored Voigt result for {assignment}: R²={result.get('avg_r_squared', 0):.3f}")
                    return result

        # Fallback: check if integrator has fitted_peaks
        if hasattr(self.integrator, 'fitted_peaks') and self.integrator.fitted_peaks:
            for result in self.integrator.fitted_peaks:
                if result.get('assignment') == assignment:
                    print(f"✅ Found integrator Voigt result for {assignment}")
                    return result

        print(f"⚠️ No stored Voigt result found for {assignment}")
        return None

    def zoom_to_peak(self, peak_data):
        """Zoom the main spectrum plot to focus on the selected peak"""
        try:
            # Extract peak coordinates using same field names as click detection
            peak_x = peak_data.get('Position_X', peak_data.get('X_ppm', peak_data.get('ppm_x')))
            peak_y = peak_data.get('Position_Y', peak_data.get('Y_ppm', peak_data.get('ppm_y')))

            if peak_x is None or peak_y is None:
                print(f"Warning: Could not find peak coordinates for zoom")
                return

            # Define zoom window (adjust these values for desired zoom level)
            zoom_window_x = 0.3  # ±0.5 ppm in X dimension
            zoom_window_y = 2.5  # ±10 ppm in Y dimension

            # Calculate zoom bounds
            x_min = peak_x - zoom_window_x
            x_max = peak_x + zoom_window_x
            y_min = peak_y - zoom_window_y
            y_max = peak_y + zoom_window_y

            # Apply zoom to the plot
            self.ax.set_xlim(x_max, x_min)  # Inverted for NMR convention
            self.ax.set_ylim(y_max, y_min)  # Inverted for NMR convention

            # Redraw the plot
            self.canvas.draw_idle()

            print(f"✅ Zoomed to peak at ({peak_x:.3f}, {peak_y:.1f}) ppm")

        except Exception as e:
            print(f"Error zooming to peak: {e}")

    def update_peak_info(self, peak_data):
        """Update the peak information panel"""
        self.peak_info_text.config(state=tk.NORMAL)
        self.peak_info_text.delete(1.0, tk.END)

        # Generate peak information text
        info_text = f"""PEAK INFORMATION
{'='*30}

Assignment: {peak_data.get('Assignment', 'Unknown')}
Position X: {peak_data.get('Position_X', 0):.3f} ppm
Position Y: {peak_data.get('Position_Y', 0):.3f} ppm

INTEGRATION RESULTS:
Height: {peak_data.get('Height', 0):.2e}
Volume: {peak_data.get('Volume', 0):.2e}
SNR: {peak_data.get('SNR', 0):.2f}
Quality: {peak_data.get('Quality', 'Unknown')}

PEAK PROPERTIES:
Detected: {peak_data.get('Detected', False)}
Line Width: {peak_data.get('Line_Width', 'N/A')}
"""

        if 'Voigt_Amplitude' in peak_data:
            info_text += f"\nVOIGT FIT PARAMETERS:\n"
            info_text += f"Amplitude: {peak_data.get('Voigt_Amplitude', 'N/A')}\n"
            info_text += f"Center: {peak_data.get('Voigt_Center', 'N/A')}\n"
            info_text += f"Gamma: {peak_data.get('Voigt_Gamma', 'N/A')}\n"
            info_text += f"Sigma: {peak_data.get('Voigt_Sigma', 'N/A')}\n"
            info_text += f"R-squared: {peak_data.get('R_squared', 'N/A')}\n"

        self.peak_info_text.insert(tk.END, info_text)
        self.peak_info_text.config(state=tk.DISABLED)

    def highlight_peak_in_plot(self, peak_data):
        """Highlight the selected peak in the spectrum plot"""
        # This would add a highlight or marker to the selected peak
        # Implementation depends on the plotting framework used
        pass

    def previous_spectrum(self):
        """Navigate to previous spectrum"""
        messagebox.showinfo("Navigation", "Previous spectrum navigation not yet implemented")

    def next_spectrum(self):
        """Navigate to next spectrum"""
        messagebox.showinfo("Navigation", "Next spectrum navigation not yet implemented")

    def center_view(self):
        """Center the view on the spectrum"""
        if hasattr(self, 'ax'):
            self.ax.autoscale()
            self.canvas.draw()

    def switch_to_voigt_analysis(self):
        """Switch to Voigt analysis tab and populate with current peak data"""
        if not self.selected_peak:
            messagebox.showwarning("No Peak Selected", "Please select a peak first by clicking on the spectrum or selecting from the peak list.")
            return

        # Switch to Voigt analysis tab
        self.notebook.select(self.voigt_tab)

        # Populate the Voigt analysis tab with current peak data
        self.populate_voigt_analysis_tab()

    def show_voigt_analysis_old(self):
        """Show detailed Voigt analysis with proper tab system like main window"""
        if not self.selected_peak:
            messagebox.showwarning("No Peak Selected", "Please select a peak first by clicking on the spectrum or selecting from the peak list.")
            return

        try:
            # Create Voigt analysis window with tab system
            voigt_window = tk.Toplevel(self.window)
            voigt_window.title(f"Voigt Analysis - {self.selected_peak.get('Assignment', 'Unknown Peak')}")
            voigt_window.geometry("1200x800")
            voigt_window.minsize(1000, 600)
            # Allow Voigt window to be independent for better interaction
            # voigt_window.transient(self.window)

            # ENHANCEMENT: Ensure enhanced fitter has current GUI parameters for consistent display
            if (hasattr(self.main_gui, 'integrator') and 
                hasattr(self.main_gui.integrator, 'enhanced_fitter') and
                hasattr(self.main_gui.integrator.enhanced_fitter, 'set_gui_parameters')):
                
                if hasattr(self.main_gui.integrator, 'fitting_parameters'):
                    self.main_gui.integrator.enhanced_fitter.set_gui_parameters(
                        self.main_gui.integrator.fitting_parameters
                    )
                    print(f"   📊 Spectrum browser: GUI parameters synchronized for Voigt analysis display")

            # Create notebook with tabs (same as main window)
            notebook = ttk.Notebook(voigt_window)
            notebook.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)

            # Tab 1: Main spectrum view
            main_tab = ttk.Frame(notebook)
            notebook.add(main_tab, text="🔬 Main Spectrum")

            # Create main spectrum plot
            fig_main, ax_main = plt.subplots(1, 1, figsize=(8, 6))
            canvas_main = FigureCanvasTkAgg(fig_main, main_tab)
            canvas_main.get_tk_widget().pack(fill=tk.BOTH, expand=True)
            toolbar_main = NavigationToolbar2Tk(canvas_main, main_tab)
            toolbar_main.update()

            # Tab 2: Voigt analysis (same layout as main window)
            voigt_tab = ttk.Frame(notebook)
            notebook.add(voigt_tab, text="📈 Voigt Analysis")

            # Create 2x2 subplot layout like main window
            fig_voigt, axes_voigt = plt.subplots(2, 2, figsize=(3, 6))
            canvas_voigt = FigureCanvasTkAgg(fig_voigt, voigt_tab)
            canvas_voigt.get_tk_widget().pack(fill=tk.BOTH, expand=True)
            toolbar_voigt = NavigationToolbar2Tk(canvas_voigt, voigt_tab)
            toolbar_voigt.update()

            # Perform Voigt analysis and populate tabs
            self.populate_voigt_analysis_tabs(ax_main, axes_voigt, canvas_main, canvas_voigt)

        except Exception as e:
            messagebox.showerror("Error", f"Could not show Voigt analysis: {str(e)}")
            import traceback
            traceback.print_exc()

    def populate_voigt_analysis_tabs(self, ax_main, axes_voigt, canvas_main, canvas_voigt):
        """Populate the Voigt analysis tabs with actual data and fits"""
        try:
            # Get peak coordinates
            # STANDARDIZED coordinate and assignment extraction (primary first)
            peak_x = (self.selected_peak.get('ppm_x') or
                     self.selected_peak.get('Position_X') or
                     self.selected_peak.get('X_ppm') or 0)
            peak_y = (self.selected_peak.get('ppm_y') or
                     self.selected_peak.get('Position_Y') or
                     self.selected_peak.get('Y_ppm') or 0)
            assignment = (self.selected_peak.get('assignment') or
                         self.selected_peak.get('Assignment') or
                         'Unknown')


            if peak_x is None or peak_y is None:
                # Show error message in plots
                ax_main.text(0.5, 0.5, "❌ Peak coordinates not found",
                           ha='center', va='center', transform=ax_main.transAxes, fontsize='small')
                for ax in axes_voigt.flat:
                    ax.text(0.5, 0.5, "❌ No data available",
                           ha='center', va='center', transform=ax.transAxes, fontsize='small')
                canvas_main.draw()
                canvas_voigt.draw()
                return

            # Tab 1: Main spectrum with zoomed view around peak
            if hasattr(self.integrator, 'nmr_data') and self.integrator.nmr_data is not None:
                data = self.integrator.nmr_data
                ppm_x = getattr(self.integrator, 'ppm_x_axis', None)
                ppm_y = getattr(self.integrator, 'ppm_y_axis', None)

                if ppm_x is not None and ppm_y is not None:
                    # Show zoomed spectrum around peak
                    zoom_x = 1.0  # ±1 ppm
                    zoom_y = 20.0  # ±20 ppm

                    # Calculate data indices for zoom region
                    x_center_idx = np.argmin(np.abs(ppm_x - peak_x))
                    y_center_idx = np.argmin(np.abs(ppm_y - peak_y))

                    x_range = int(zoom_x * len(ppm_x) / (ppm_x[0] - ppm_x[-1]))
                    y_range = int(zoom_y * len(ppm_y) / (ppm_y[0] - ppm_y[-1]))

                    x_min = max(0, x_center_idx - x_range)
                    x_max = min(len(ppm_x), x_center_idx + x_range)
                    y_min = max(0, y_center_idx - y_range)
                    y_max = min(len(ppm_y), y_center_idx + y_range)

                    # Extract zoom region
                    zoom_data = data[y_min:y_max, x_min:x_max]
                    zoom_ppm_x = ppm_x[x_min:x_max]
                    zoom_ppm_y = ppm_y[y_min:y_max]

                    # Plot contour
                    levels = np.geomspace(np.max(np.abs(zoom_data)) * 0.1,
                                        np.max(np.abs(zoom_data)) * 0.8, 10)
                    ax_main.contour(zoom_ppm_x, zoom_ppm_y, zoom_data, levels=levels,
                                  cmap='viridis', linewidths=0.5)

                    # Mark the peak
                    ax_main.scatter([peak_x], [peak_y], c='red', s=50, marker='x',
                                  linewidths=3, label=f'Peak {assignment}')

                    ax_main.set_xlabel('¹H (ppm)')
                    ax_main.set_ylabel('¹⁵N (ppm)')
                    ax_main.set_title(f"Zoomed Spectrum - Peak {assignment}")
                    ax_main.invert_xaxis()
                    ax_main.invert_yaxis()
                    ax_main.legend()

                    # Tab 2: Voigt analysis with 1D fits (2x2 layout)
                    # CRITICAL FIX: Use stored results from series processing
                    assignment_key = self.selected_peak.get('Assignment', self.selected_peak.get('assignment', assignment))
                    stored_result = self.find_stored_voigt_result(assignment_key)
                    print(f"🔍 Using stored result for Voigt analysis: {stored_result is not None}")

                    self.populate_1d_voigt_analysis(axes_voigt, zoom_data, zoom_ppm_x, zoom_ppm_y,
                                                   peak_x, peak_y, assignment, stored_voigt_result=stored_result)

            canvas_main.draw()
            canvas_voigt.draw()

        except Exception as e:
            print(f"Error populating Voigt analysis tabs: {e}")
            import traceback
            traceback.print_exc()

## here goes populate_1d_voigt_analysis

    def populate_1d_voigt_analysis(self, axes, data, ppm_x, ppm_y, peak_x, peak_y, assignment, stored_voigt_result=None):
        """
        Populates the 1D Voigt analysis plots, including the CORRECT residual
        on a shared Y-axis.
        """
        try:
            # --- Data Preparation ---
            # Find center indices and define a data window to display
            x_idx = np.argmin(np.abs(ppm_x - peak_x))
            y_idx = np.argmin(np.abs(ppm_y - peak_y))
            window_size = 20
            x_min, x_max = max(0, x_idx - window_size), min(len(ppm_x), x_idx + window_size + 1)
            y_min, y_max = max(0, y_idx - window_size), min(len(ppm_y), y_idx + window_size + 1)

            # This is the raw data from the currently browsed spectrum
            current_x_ppm = ppm_x[x_min:x_max]
            current_x_slice = data[y_idx, x_min:x_max]
            current_y_ppm = ppm_y[y_min:y_max]
            current_y_slice = data[y_min:y_max, x_idx]

            # --- Plot 1: X-Dimension (Top Plot) ---
            ax_x = axes[0]
            ax_x.clear()
            ax_x.plot(current_x_ppm, current_x_slice, 'b-', label='Current Raw Data', alpha=0.5)

            if stored_voigt_result and stored_voigt_result.get('x_fit'):
                x_fit = stored_voigt_result['x_fit']
                # These are the results from the original fit
                stored_curve = x_fit.get('fitted_curve')
                stored_ppm = x_fit.get('ppm_scale')
                stored_data = x_fit.get('cross_section') # The data used for the fit

                if stored_curve is not None and stored_ppm is not None:
                    # Plot the stored fitted curve
                    r_squared = x_fit.get('r_squared', 0)
                    ax_x.plot(stored_ppm, stored_curve, 'r--', linewidth=2, label=f"Stored Fit (R²={r_squared:.3f})")

                    # --- INTERPOLATED RESIDUAL CALCULATION FOR ALIGNMENT ---
                    # Check if the original data slice was stored with the fit
                    if stored_data is not None and len(stored_data) == len(stored_curve):
                        from scipy.interpolate import interp1d

                        # Create interpolation functions for proper x-axis alignment
                        stored_data_interp = interp1d(stored_ppm, stored_data,
                                                    bounds_error=False, fill_value='extrapolate')
                        stored_curve_interp = interp1d(stored_ppm, stored_curve,
                                                     bounds_error=False, fill_value='extrapolate')

                        # Interpolate to current ppm scale for alignment
                        stored_data_aligned = stored_data_interp(current_x_ppm)
                        stored_curve_aligned = stored_curve_interp(current_x_ppm)

                        # Calculate residual on aligned scales
                        residual = stored_data_aligned - stored_curve_aligned
                        ax_x.plot(current_x_ppm, residual, 'g-', linewidth=1, alpha=0.8,
                                 label='Aligned Residual')

                        # Replace stored fit plot with aligned version
                        ax_x.lines[-2].remove()  # Remove the previous stored fit line
                        ax_x.plot(current_x_ppm, stored_curve_aligned, 'r--', linewidth=2,
                                 label=f"Aligned Fit (R²={r_squared:.3f})")

                    # Set plot limits based on the fit's data range for clarity
                    ax_x.set_xlim(current_x_ppm.max(), current_x_ppm.min())

            ax_x.axhline(y=0, color='k', linestyle='--', alpha=0.5)
            ax_x.set_title(f'X-Dimension Slice (¹H) at {peak_y:.1f} ppm', fontsize='small')
            #ax_x.set_ylabel('Intensity / Residual', fontsize='small')
            ax_x.set_ylabel('Intensity / Residual', fontsize='small', rotation=90, labelpad=10)
            ax_x.grid(True, alpha=0.4)
            ax_x.legend(fontsize='small')

            # --- Plot 2: Y-Dimension (Bottom Plot) ---
            ax_y = axes[1]
            ax_y.clear()
            ax_y.plot(current_y_ppm, current_y_slice, 'b-', label='Current Raw Data', alpha=0.5)

            if stored_voigt_result and stored_voigt_result.get('y_fit'):
                y_fit = stored_voigt_result['y_fit']
                stored_curve = y_fit.get('fitted_curve')
                stored_ppm = y_fit.get('ppm_scale')
                stored_data = y_fit.get('cross_section')

                if stored_curve is not None and stored_ppm is not None:
                    ax_y.plot(stored_ppm, stored_curve, 'r--', linewidth=2, label=f"Stored Fit (R²={y_fit.get('r_squared', 0):.3f})")

                    #if stored_data is not None and len(stored_data) == len(stored_curve):
                    #    residual = stored_data - stored_curve
                    #    ax_y.plot(stored_ppm, residual, 'g-', linewidth=1, alpha=0.8, label='Original Residual')
## gm added

                    if stored_data is not None and len(stored_data) == len(stored_curve):
                        from scipy.interpolate import interp1d

                        # Create interpolation functions for proper x-axis alignment
                        stored_data_interp = interp1d(stored_ppm, stored_data,
                                                    bounds_error=False, fill_value='extrapolate')
                        stored_curve_interp = interp1d(stored_ppm, stored_curve,
                                                     bounds_error=False, fill_value='extrapolate')

                        # Interpolate to current ppm scale for alignment
                        stored_data_aligned = stored_data_interp(current_y_ppm)
                        stored_curve_aligned = stored_curve_interp(current_y_ppm)

                        # Calculate residual on aligned scales
                        residual = stored_data_aligned - stored_curve_aligned
                        ax_y.plot(current_y_ppm, residual, 'g-', linewidth=1, alpha=0.8,
                                 label='Aligned Residual')

                        # Replace stored fit plot with aligned version
                        ax_y.lines[-2].remove()  # Remove the previous stored fit line
                        ax_y.plot(current_y_ppm, stored_curve_aligned, 'r--', linewidth=2,
                                 label=f"Aligned Fit (R²={y_fit.get('r_squared', 0):.3f})")

## gm added
                    #ax_y.set_xlim(stored_ppm.max(), stored_ppm.min())
                    ax_y.set_xlim(current_y_ppm.max(), current_y_ppm.min())

            ax_y.axhline(y=0, color='k', linestyle='--', alpha=0.5)
            ax_y.set_title(f'Y-Dimension Slice (¹⁵N/¹³C) at {peak_x:.3f} ppm', fontsize='small')
            ax_y.set_xlabel('¹⁵N/¹³C (ppm)', fontsize='small')
            #ax_y.set_ylabel('Intensity / Residual', fontsize='small')
            ax_y.set_ylabel('Intensity / Residual', fontsize='small', rotation=90, labelpad=10)
            ax_y.grid(True, alpha=0.4)
            ax_y.legend(fontsize='small')

            #self.fig_voigt.tight_layout(pad=3.0)
            # Fix overlapping Y-axis labels and cramped spacing
            self.fig_voigt.subplots_adjust(
                  left=0.15,      # More space for Y-axis labels
                  right=0.95,     # Reasonable right margin
                  top=0.92,       # Space for top title
                  bottom=0.12,    # Space for bottom X-axis label
                  hspace=0.4      # More vertical space between subplots
            )


        except Exception as e:
            print(f"Error in 1D Voigt analysis: {e}")
            import traceback
            traceback.print_exc()
            for ax in axes.flat:
                ax.clear()
                ax.text(0.5, 0.5, f"Analysis Error:\n{e}", ha='center', va='center', transform=ax.transAxes, color='red')
##

    def show_peak_summary(self):
        """Show summary of all peaks in this spectrum"""
        # Create summary window
        summary_window = tk.Toplevel(self.window)
        summary_window.title(f"Peak Summary - {self.spectrum_name}")
        summary_window.geometry("800x600")
        # Allow summary window to be independent for better interaction
        # summary_window.transient(self.window)

        # Summary content
        text_widget = tk.Text(summary_window, wrap=tk.WORD, font=('Courier', 10))
        scrollbar = ttk.Scrollbar(summary_window, orient=tk.VERTICAL, command=text_widget.yview)
        text_widget.configure(yscrollcommand=scrollbar.set)

        # Generate summary
        summary_text = self.generate_peak_summary()
        text_widget.insert(tk.END, summary_text)
        text_widget.config(state=tk.DISABLED)

        text_widget.pack(side=tk.LEFT, fill=tk.BOTH, expand=True, padx=10, pady=10)
        scrollbar.pack(side=tk.RIGHT, fill=tk.Y, pady=10)

    def generate_peak_summary(self):
        """Generate peak summary text"""
        integration_results = self.result_data.get('integration_results', [])

        summary = f"""PEAK SUMMARY REPORT
{'='*60}

Spectrum: {self.spectrum_name}
Status: {self.result_data.get('status', 'Unknown').title()}
Detection: {self.result_data.get('detected_peaks', 0)}/{self.result_data.get('total_peaks', 0)} peaks ({self.result_data.get('detection_rate', 0):.1f}%)

"""

        if integration_results:
            summary += f"INTEGRATION STATISTICS:\n"
            summary += f"{'='*40}\n"

            # Calculate statistics
            heights = [p.get('Height', 0) for p in integration_results]
            snrs = [p.get('SNR', 0) for p in integration_results]

            summary += f"Total Peaks Integrated: {len(integration_results)}\n"
            if heights:
                summary += f"Height Range: {min(heights):.2e} - {max(heights):.2e}\n"
                summary += f"Average Height: {np.mean(heights):.2e}\n"
            if snrs:
                summary += f"SNR Range: {min(snrs):.1f} - {max(snrs):.1f}\n"
                summary += f"Average SNR: {np.mean(snrs):.1f}\n"
                good_quality = sum(1 for snr in snrs if snr >= 5)
                summary += f"Good Quality (SNR ≥ 5): {good_quality}/{len(snrs)} ({good_quality/len(snrs)*100:.1f}%)\n"

            summary += f"\nPEAK DETAILS:\n"
            summary += f"{'='*40}\n"
            summary += f"{'#':<3} {'Assignment':<15} {'Height':<12} {'SNR':<8} {'Quality':<10}\n"
            summary += f"{'-'*60}\n"

            for i, peak in enumerate(integration_results, 1):
                assignment = peak.get('Assignment', 'Unknown')[:14]
                height = f"{peak.get('Height', 0):.2e}"
                snr = f"{peak.get('SNR', 0):.1f}"
                quality = peak.get('Quality', 'Unknown')[:9]

                summary += f"{i:<3} {assignment:<15} {height:<12} {snr:<8} {quality:<10}\n"

        return summary

    def export_spectrum_data(self):
        """Export spectrum data and analysis"""
        try:
            from tkinter.filedialog import asksaveasfilename

            # Ask user where to save
            file_path = asksaveasfilename(
                title="Export Spectrum Data",
                filetypes=[("CSV files", "*.csv"), ("Text files", "*.txt"), ("All files", "*.*")],
                defaultextension=".csv",
                initialname=f"{self.spectrum_name}_export.csv"
            )

            if file_path:
                # Export integration results
                integration_results = self.result_data.get('integration_results', [])

                if integration_results:
                    import pandas as pd
                    df = pd.DataFrame(integration_results)
                    df.to_csv(file_path, index=False, float_format='%.6f')
                    messagebox.showinfo("Export Complete", f"Spectrum data exported to:\n{file_path}")
                else:
                    messagebox.showwarning("No Data", "No integration results to export")

        except Exception as e:
            messagebox.showerror("Export Error", f"Failed to export data:\n{str(e)}")

    # =================== INTERACTIVE FUNCTIONALITY ===================

    def update_spectrum_plot(self):
        """Update spectrum plot with current contour settings"""
        if self.integrator and hasattr(self.integrator, 'nmr_data') and self.integrator.nmr_data is not None:
            self.plot_spectrum()

    def preset_low(self):
        """Low detail contour preset"""
        self.contour_levels.set(10)
        self.contour_min.set(0.2)
        self.contour_increment.set(2.0)
        self.update_spectrum_plot()

    def preset_medium(self):
        """Medium detail contour preset"""
        self.contour_levels.set(20)
        self.contour_min.set(0.1)
        self.contour_increment.set(1.0)
        self.update_spectrum_plot()

    def preset_high(self):
        """High detail contour preset"""
        self.contour_levels.set(40)
        self.contour_min.set(0.05)
        self.contour_increment.set(0.5)
        self.update_spectrum_plot()

    def preset_auto(self):
        """Auto contour preset based on data"""
        try:
            if self.integrator and hasattr(self.integrator, 'nmr_data') and self.integrator.nmr_data is not None:
                data = np.abs(self.integrator.nmr_data)
                noise_level = np.mean(data[data < np.percentile(data, 20)])
                signal_level = np.max(data)

                # Calculate optimal settings
                snr = signal_level / noise_level if noise_level > 0 else 100

                # Optimized adaptive levels for responsiveness
                if snr > 50:
                    levels = 10  # Reduced from 30
                    min_level = 0.15  # Higher for less complexity
                elif snr > 20:
                    levels = 8   # Reduced from 20
                    min_level = 0.25  # Higher for less complexity
                else:
                    levels = 6   # Reduced from 15
                    min_level = 0.4   # Higher for less complexity

                self.contour_levels.set(levels)
                self.contour_min.set(min_level)
                self.contour_increment.set(1.0)
                self.update_spectrum_plot()

        except Exception as e:
            print(f"Auto preset failed: {e}")
            self.preset_medium()  # Fallback

    def plot_peak_overlays(self):
        """Plot peak markers on spectrum"""
        if not hasattr(self.integrator, 'fitted_peaks'):
            return

        try:
            for peak in self.integrator.fitted_peaks:
                if peak.get('detected', False):
                    x = peak.get('ppm_x', 0)
                    y = peak.get('ppm_y', 0)
                    assignment = peak.get('assignment', '')

                    # Plot marker
                    self.ax.plot(x, y, 'ro', markersize=6, alpha=0.8)

                    # Add label
                    if assignment:
                        self.ax.annotate(assignment, (x, y), xytext=(5, 5),
                                       textcoords='offset points', fontsize='small',
                                       bbox=dict(boxstyle='round,pad=0.2', facecolor='yellow', alpha=0.7))
        except Exception as e:
            print(f"Peak overlay error: {e}")

    def on_plot_click(self, event):
        """Handle plot click events for peak selection"""
        if not event.inaxes or not hasattr(self.integrator, 'fitted_peaks'):
            return

        try:
            x_click, y_click = event.xdata, event.ydata
            min_dist = float('inf')
            selected_peak = None

            # Find closest peak
            for peak in self.integrator.fitted_peaks:
                if peak.get('detected', False):
                    x_peak = peak.get('ppm_x', 0)
                    y_peak = peak.get('ppm_y', 0)

                    dist = np.sqrt((x_click - x_peak)**2 + (y_click - y_peak)**2)
                    if dist < min_dist and dist < 1.0:  # Within 1 ppm
                        min_dist = dist
                        selected_peak = peak

            if selected_peak:
                self.select_peak_by_data(selected_peak)

        except Exception as e:
            print(f"Plot click error: {e}")

    def on_zoom(self, event):
        """Handle mouse wheel zoom"""
        if not event.inaxes:
            return

        try:
            scale = 1.1 if event.button == 'down' else 1/1.1

            xlim = self.ax.get_xlim()
            ylim = self.ax.get_ylim()

            x_center, y_center = event.xdata, event.ydata

            # Calculate new limits
            x_range = (xlim[1] - xlim[0]) * scale
            y_range = (ylim[1] - ylim[0]) * scale

            new_xlim = (x_center - x_range/2, x_center + x_range/2)
            new_ylim = (y_center - y_range/2, y_center + y_range/2)

            self.ax.set_xlim(new_xlim)
            self.ax.set_ylim(new_ylim)
            self.canvas.draw_idle()

        except Exception as e:
            print(f"Zoom error: {e}")

    def on_mouse_move(self, event):
        """Handle mouse movement for cursor info"""
        # This could be extended to show cursor coordinates
        pass

    def center_view(self):
        """Reset view to center on full spectrum"""
        if self.original_xlim and self.original_ylim:
            self.ax.set_xlim(self.original_xlim)
            self.ax.set_ylim(self.original_ylim)
            self.canvas.draw_idle()

    def select_peak_by_data(self, peak_data):
        """Select a peak and update displays"""
        self.selected_peak = peak_data

        # Update peak info display
        peak_info = f"""SELECTED PEAK INFORMATION
{'='*40}

Assignment: {peak_data.get('assignment', 'Unknown')}
Position: {peak_data.get('ppm_x', 0):.3f} x {peak_data.get('ppm_y', 0):.1f} ppm
Detected: {'Yes' if peak_data.get('detected', False) else 'No'}
SNR: {peak_data.get('snr', 0):.1f}

"""

        # Find corresponding integration result
        if hasattr(self.integrator, 'integration_results'):
            assignment = peak_data.get('assignment', '')
            for result in self.integrator.integration_results:
                if result.get('Assignment') == assignment:
                    peak_info += f"INTEGRATION RESULTS:\n{'-'*20}\n"
                    peak_info += f"Height: {result.get('Height', 0):.2e}\n"
                    peak_info += f"Volume: {result.get('Volume', 0):.2e}\n"
                    peak_info += f"Quality: {result.get('Quality', 'Unknown')}\n"
                    break

        # Update text display
        self.peak_info_text.config(state=tk.NORMAL)
        self.peak_info_text.delete(1.0, tk.END)
        self.peak_info_text.insert(tk.END, peak_info)
        self.peak_info_text.config(state=tk.DISABLED)

    def previous_spectrum(self):
        """Navigate to previous spectrum in series"""
        messagebox.showinfo("Navigation", "Previous spectrum navigation will be implemented")

    def next_spectrum(self):
        """Navigate to next spectrum in series"""
        messagebox.showinfo("Navigation", "Next spectrum navigation will be implemented")

    def generate_voigt_analysis(self, peak_data):
        """Generate detailed Voigt analysis for a peak"""
        analysis = f"""VOIGT FIT ANALYSIS REPORT
{'='*60}

Peak Information:
• Assignment: {peak_data.get('assignment', 'Unknown')}
• Position: {peak_data.get('ppm_x', 0):.3f} × {peak_data.get('ppm_y', 0):.1f} ppm
• Detection Status: {'✅ Detected' if peak_data.get('detected', False) else '❌ Not Detected'}
• Signal-to-Noise Ratio: {peak_data.get('snr', 0):.1f}

"""

        # Find corresponding integration result
        integration_found = False
        if hasattr(self.integrator, 'integration_results') and self.integrator.integration_results:
            assignment = peak_data.get('assignment', '')
            for result in self.integrator.integration_results:
                if result.get('Assignment') == assignment:
                    integration_found = True
                    analysis += f"""Integration Results:
• Height: {result.get('Height', 0):.2e}
• Volume: {result.get('Volume', 0):.2e}
• Quality: {result.get('Quality', 'Unknown')}
• SNR: {result.get('SNR', 0):.1f}

"""
                    break

        if not integration_found:
            analysis += """Integration Results:
• No integration data available for this peak

"""

        # Check for Voigt fit data
        voigt_found = False
        if hasattr(self.integrator, 'voigt_fits') and self.integrator.voigt_fits:
            assignment = peak_data.get('assignment', '')
            for fit in self.integrator.voigt_fits:
                if fit.get('assignment') == assignment or fit.get('peak_id') == assignment:
                    voigt_found = True
                    # Use local R-squared for quality assessment if available
                    r_squared_local = fit.get('avg_r_squared_local', fit.get('r_squared', 0))
                    r_squared_global = fit.get('avg_r_squared_global', fit.get('r_squared', 0))
                    quality = 'Excellent' if r_squared_local > 0.9 else 'Good' if r_squared_local > 0.8 else 'Fair' if r_squared_local > 0.6 else 'Poor'

                    analysis += f"""Voigt Fit Parameters:
• R-squared (Local): {r_squared_local:.3f} ({quality})
• R-squared (Global): {r_squared_global:.3f}
• Amplitude: {fit.get('amplitude', 0):.2e}
• Center X: {fit.get('center_x', 0):.3f} ppm
• Center Y: {fit.get('center_y', 0):.1f} ppm
• Sigma X: {fit.get('sigma_x', 0):.3f} ppm
• Sigma Y: {fit.get('sigma_y', 0):.1f} ppm
• Gamma X: {fit.get('gamma_x', 0):.3f} ppm
• Gamma Y: {fit.get('gamma_y', 0):.1f} ppm

Fit Quality Assessment:
• Overall Quality: {quality}
• Convergence: {'✅ Converged' if fit.get('converged', False) else '❌ Failed to converge'}
• Iterations: {fit.get('iterations', 'Unknown')}
• Chi-square: {fit.get('chi_squared', 'Not available')}

"""
                    break

        if not voigt_found:
            analysis += """Voigt Fit Parameters:
• No Voigt fit data available for this peak
• Peak may not have been fitted or fitting may have failed

"""

        # Add spectrum context
        analysis += f"""Spectrum Context:
• File: {self.spectrum_name}
• Total Peaks in Spectrum: {self.result_data.get('total_peaks', 0)}
• Detected Peaks: {self.result_data.get('detected_peaks', 0)}
• Detection Rate: {self.result_data.get('detection_rate', 0):.1f}%
• Noise Level: {self.result_data.get('noise_level', 0):.2e}
• Processing Status: {self.result_data.get('status', 'Unknown').title()}

Analysis Notes:
• This analysis is based on the series integration results
• Voigt fitting combines Gaussian and Lorentzian line shapes
• R-squared values > 0.8 indicate reliable fits
• SNR values > 5 typically indicate good peak quality
• Peak assignments are from the reference peak list

"""

        return analysis

    # =================== VOIGT ANALYSIS NAVIGATION METHODS ===================
    def voigt_navigate_to_peak(self):
        """Navigate to the specified peak in Voigt analysis"""
        try:
            peak_number = self.voigt_peak_number.get()
            integration_results = self.result_data.get('integration_results', [])

            if 1 <= peak_number <= len(integration_results):
                # Get peak data and select it
                peak_data = integration_results[peak_number - 1]
                self.select_peak(peak_data)

                # Update the peak info display
                assignment = peak_data.get('Assignment', f'Peak_{peak_number}')
                self.voigt_peak_info_label.config(text=f"Peak: {peak_number}/{len(integration_results)} - {assignment}")

                # Update the Voigt analysis display
                self.populate_voigt_analysis_tab()

                print(f"🎯 Navigated to peak {peak_number}: {assignment}")
            else:
                print(f"⚠️ Invalid peak number: {peak_number}")

        except Exception as e:
            print(f"❌ Error navigating to peak: {e}")

    def voigt_prev_peak(self):
        """Navigate to previous peak in Voigt analysis"""
        try:
            current_peak = self.voigt_peak_number.get()
            if current_peak > 1:
                self.voigt_peak_number.set(current_peak - 1)
                self.voigt_navigate_to_peak()
        except Exception as e:
            print(f"❌ Error going to previous peak: {e}")

    def voigt_next_peak(self):
        """Navigate to next peak in Voigt analysis"""
        try:
            current_peak = self.voigt_peak_number.get()
            integration_results = self.result_data.get('integration_results', [])
            if current_peak < len(integration_results):
                self.voigt_peak_number.set(current_peak + 1)
                self.voigt_navigate_to_peak()
        except Exception as e:
            print(f"❌ Error going to next peak: {e}")

    def voigt_center_on_peak(self):
        """Center the main spectrum view on the current peak"""
        try:
            if self.selected_peak:
                # STANDARDIZED coordinate extraction (primary first)
                peak_x = (self.selected_peak.get('pmp_x') or
                         self.selected_peak.get('Position_X') or 0)
                peak_y = (self.selected_peak.get('pmp_y') or
                         self.selected_peak.get('Position_Y') or 0)

                # Switch to Main Spectrum tab first
                self.notebook.select(0)  # Select first tab (Main Spectrum)

                # Center the plot on the peak
                if hasattr(self, 'ax') and self.ax:
                    # Set zoom window around peak (±0.5 ppm in both dimensions)
                    x_margin = 0.3
                    y_margin = 2.5

                    self.ax.set_xlim(peak_x + x_margin, peak_x - x_margin)  # Inverted for 1H
                    self.ax.set_ylim(peak_y - y_margin, peak_y + y_margin)

                    # Mark the peak
                    self.ax.plot(peak_x, peak_y, 'ro', markersize=8, markeredgewidth=2,
                               markeredgecolor='white', alpha=0.8)

                    self.canvas.draw()
                    print(f"🎯 Centered on peak at {peak_x:.3f}, {peak_y:.1f} ppm")

        except Exception as e:
            print(f"❌ Error centering on peak: {e}")

    def voigt_zoom_to_peak(self):
        """Zoom to the current peak in the main spectrum view"""
        try:
            if self.selected_peak:
                peak_x = self.selected_peak.get('Position_X', 0)
                peak_y = self.selected_peak.get('Position_Y', 0)

                # Switch to Main Spectrum tab first
                self.notebook.select(0)  # Select first tab (Main Spectrum)

                # Tight zoom around peak (±0.2 ppm in both dimensions)
                if hasattr(self, 'ax') and self.ax:
                    x_margin = 0.3
                    y_margin = 2.5

                    self.ax.set_xlim(peak_x + x_margin, peak_x - x_margin)  # Inverted for 1H
                    self.ax.set_ylim(peak_y - y_margin, peak_y + y_margin)

                    # Mark the peak prominently
                    self.ax.plot(peak_x, peak_y, 'ro', markersize=12, markeredgewidth=3,
                               markeredgecolor='yellow', alpha=0.9)

                    self.canvas.draw()
                    print(f"🔍 Zoomed to peak at {peak_x:.3f}, {peak_y:.1f} ppm")

        except Exception as e:
            print(f"❌ Error zooming to peak: {e}")

    def voigt_show_peak_analysis(self):
        """Show detailed analysis of the current peak"""
        try:
            if self.selected_peak:
                assignment = self.selected_peak.get('Assignment', 'Unknown')

                # Generate detailed analysis
                analysis_text = self.generate_detailed_peak_analysis()

                # Show in a dialog window
                analysis_window = tk.Toplevel(self.window)
                analysis_window.title(f"Peak Analysis - {assignment}")
                analysis_window.geometry("600x500")

                # Create scrollable text area
                text_frame = ttk.Frame(analysis_window)
                text_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

                analysis_display = tk.Text(text_frame, wrap=tk.WORD, font=('Courier', 9))
                scrollbar_analysis = ttk.Scrollbar(text_frame, orient=tk.VERTICAL, command=analysis_display.yview)
                analysis_display.configure(yscrollcommand=scrollbar_analysis.set)

                analysis_display.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
                scrollbar_analysis.pack(side=tk.RIGHT, fill=tk.Y)

                # Insert analysis text
                analysis_display.insert(tk.END, analysis_text)
                analysis_display.config(state=tk.DISABLED)

                # Close button
                ttk.Button(analysis_window, text="Close",
                          command=analysis_window.destroy).pack(pady=10)

                print(f"🔬 Showing detailed analysis for {assignment}")
            else:
                print("⚠️ No peak selected for analysis")

        except Exception as e:
            print(f"❌ Error showing peak analysis: {e}")

#!/usr/bin/env python3
"""
GUI Components Module

This module contains all the reusable GUI components for the NMR Peak Series application.
These components are separated from the main GUI logic for better modularity and reusability.

Components:
- ScrollableFrame: Enhanced scrollable frame for controls
- EnhancedFileListFrame: File selection with metadata and preview
- AdvancedProgressDialog: Progress tracking with detailed logging
- StatisticsPanel: Statistics display and analysis

Author: Guillaume Mas
Date: 2025
"""

import tkinter as tk
from tkinter import ttk, filedialog, messagebox
import os
import glob
import json
from datetime import datetime, timedelta
from pathlib import Path

class ScrollableFrame(ttk.Frame):
    """Enhanced scrollable frame widget for control panels"""
    def __init__(self, container, *args, **kwargs):
        super().__init__(container, *args, **kwargs)

        self.canvas = tk.Canvas(self, highlightthickness=0, bg='white')
        self.scrollbar = ttk.Scrollbar(self, orient="vertical", command=self.canvas.yview)
        self.scrollable_frame = ttk.Frame(self.canvas)

        self.scrollable_frame.bind(
            "<Configure>",
            lambda e: self.canvas.configure(scrollregion=self.canvas.bbox("all"))
        )

        self.canvas.create_window((0, 0), window=self.scrollable_frame, anchor="nw")
        self.canvas.configure(yscrollcommand=self.scrollbar.set)

        self.canvas.pack(side="left", fill="both", expand=True)
        self.scrollbar.pack(side="right", fill="y")

        # Enhanced mouse wheel scrolling
        self.bind_mousewheel()

    def bind_mousewheel(self):
        def _on_mousewheel(event):
            self.canvas.yview_scroll(int(-1*(event.delta/120)), "units")

        def _bind_to_mousewheel(event):
            self.canvas.bind_all("<MouseWheel>", _on_mousewheel)

        def _unbind_from_mousewheel(event):
            self.canvas.unbind_all("<MouseWheel>")

        self.canvas.bind('<Enter>', _bind_to_mousewheel)
        self.canvas.bind('<Leave>', _unbind_from_mousewheel)

class EnhancedFileListFrame(ttk.Frame):
    """Enhanced file list frame with preview and metadata"""
    def __init__(self, parent, title, file_types=None, height=5, width=20):
        super().__init__(parent)
        self.file_types = file_types or []
        self.current_folder = None
        self.callback = None
        self.file_metadata = {}

        # Title with icon
        self.title_label = ttk.Label(self, text=f"📁 {title}", font=('TkDefaultFont', 10, 'bold'))
        self.title_label.pack(anchor=tk.W, pady=(0, 5))

        # Enhanced folder selection
        folder_frame = ttk.Frame(self)
        folder_frame.pack(fill=tk.X, pady=(0, 5))

        self.folder_button = ttk.Button(folder_frame, text=f"Select {title} Folder",
                                       command=self.select_folder, width=20)
        self.folder_button.pack(side=tk.LEFT, padx=(0, 5))

        # Refresh button
        self.refresh_button = ttk.Button(folder_frame, text="🔄", width=3,
                                       command=self.refresh_file_list, state='disabled')
        self.refresh_button.pack(side=tk.LEFT)

        # Current folder display
        self.folder_label = ttk.Label(self, text="📂 No folder selected",
                                     foreground='gray', font=('TkDefaultFont', 9))
        self.folder_label.pack(anchor=tk.W, pady=(0, 5))

        # Enhanced file list with metadata
        list_frame = ttk.Frame(self)
        list_frame.pack(fill=tk.BOTH, expand=True)

        # File listbox with custom styling
        self.file_listbox = tk.Listbox(list_frame, height=height,width=width, font=('TkDefaultFont', 8),
                                      selectmode=tk.SINGLE, activestyle='dotbox')

        # Scrollbars
        scrollbar_y = ttk.Scrollbar(list_frame, orient="vertical", command=self.file_listbox.yview)
        scrollbar_x = ttk.Scrollbar(list_frame, orient="horizontal", command=self.file_listbox.xview)

        self.file_listbox.config(yscrollcommand=scrollbar_y.set, xscrollcommand=scrollbar_x.set)

        # Grid layout for scrollbars
        self.file_listbox.grid(row=0, column=0, sticky='nsew')
        scrollbar_y.grid(row=0, column=1, sticky='ns')
        scrollbar_x.grid(row=1, column=0, sticky='ew')

        list_frame.grid_rowconfigure(0, weight=1)
        list_frame.grid_columnconfigure(0, weight=1)

        # Bind selection events
        self.file_listbox.bind('<<ListboxSelect>>', self.on_file_select)
        self.file_listbox.bind('<Double-Button-1>', self.on_file_double_click)

        # Status and metadata display
        status_frame = ttk.Frame(self)
        status_frame.pack(fill=tk.X, pady=(5, 0))

        self.status_label = ttk.Label(status_frame, text="", foreground='blue', font=('TkDefaultFont', 8))
        self.status_label.pack(side=tk.LEFT)

        self.metadata_label = ttk.Label(status_frame, text="", foreground='gray', font=('TkDefaultFont', 8))
        self.metadata_label.pack(side=tk.RIGHT)

        #File info GM added comment
        #File preview area (optional)
        self.preview_frame = ttk.LabelFrame(self, text="📋 File Info", padding=5)
        self.preview_frame.pack(fill=tk.X, pady=(5, 0))

        self.preview_text = tk.Text(self.preview_frame, height=3, width=40, font=('Courier', 8),
                                   wrap=tk.WORD, state='disabled')
        self.preview_text.pack(fill=tk.X)

    def select_folder(self):
        """Enhanced folder selection with validation"""
        initial_dir = self.current_folder if self.current_folder else os.getcwd()
        folder = filedialog.askdirectory(
            title=f"Select {self.title_label['text'].replace('📁 ', '')} Folder",
            initialdir=initial_dir
        )

        if folder:
            self.current_folder = folder
            folder_name = os.path.basename(folder)
            try:
                item_count = len(os.listdir(folder))
                self.folder_label.config(
                    text=f"📂 {folder_name} ({item_count} items)",
                    foreground='black'
                )
            except PermissionError:
                self.folder_label.config(
                    text=f"📂 {folder_name} (access denied)",
                    foreground='red'
                )

            self.refresh_button.config(state='normal')
            self.refresh_file_list()

    def refresh_file_list(self):
        """Enhanced file list refresh with metadata collection"""
        if not self.current_folder:
            return

        # Clear current list
        self.file_listbox.delete(0, tk.END)
        self.file_metadata.clear()

        try:
            # Find matching files
            files = []
            for file_type in self.file_types:
                pattern = os.path.join(self.current_folder, f"*.{file_type}")
                files.extend(glob.glob(pattern))

            # Sort files naturally
            files.sort(key=lambda x: os.path.basename(x).lower())

            # Add files to listbox with metadata
            for file_path in files:
                filename = os.path.basename(file_path)
                try:
                    # Collect metadata
                    stat_info = os.stat(file_path)
                    size_mb = stat_info.st_size / (1024 * 1024)
                    mod_time = datetime.fromtimestamp(stat_info.st_mtime)

                    self.file_metadata[filename] = {
                        'path': file_path,
                        'size_mb': size_mb,
                        'modified': mod_time,
                        'extension': os.path.splitext(filename)[1]
                    }

                    # Display format: filename (size)
                    display_name = f"{filename} ({size_mb:.1f} MB)"
                    self.file_listbox.insert(tk.END, display_name)

                except (OSError, PermissionError):
                    self.file_listbox.insert(tk.END, f"{filename} (error)")

            # Update status
            self.status_label.config(text=f"{len(files)} file(s) found")

        except Exception as e:
            self.status_label.config(text=f"Error: {str(e)}")

    def on_file_select(self, event=None):
        """Enhanced file selection with preview"""
        selection = self.file_listbox.curselection()
        if not selection or not self.current_folder:
            return

        # Get selected filename (remove size info)
        display_name = self.file_listbox.get(selection[0])
        filename = display_name.split(' (')[0]  # Remove size info

        if filename in self.file_metadata:
            metadata = self.file_metadata[filename]

            # Update metadata display
            mod_time = metadata['modified'].strftime("%Y-%m-%d %H:%M")
            self.metadata_label.config(text=f"Modified: {mod_time}")

            # Update preview
            self.update_preview(metadata)

            # Call callback if set
            if self.callback:
                self.callback(metadata['path'], filename)

    def on_file_double_click(self, event=None):
        """Handle double-click events"""
        self.on_file_select(event)

    def update_preview(self, metadata):
        """Update file preview area"""
        self.preview_text.config(state='normal')
        self.preview_text.delete(1.0, tk.END)

        preview_info = (
            f"📁 Path: {metadata['path']}\n"
            f"📏 Size: {metadata['size_mb']:.2f} MB\n"
            f"📅 Modified: {metadata['modified'].strftime('%Y-%m-%d %H:%M:%S')}"
        )

        self.preview_text.insert(1.0, preview_info)
        self.preview_text.config(state='disabled')

    def set_callback(self, callback):
        """Set callback function for file selection"""
        self.callback = callback

    def get_current_folder(self):
        """Get currently selected folder"""
        return self.current_folder

    def get_all_files(self):
        """Get all files in current folder"""
        if not self.current_folder:
            return []

        files = []
        for file_type in self.file_types:
            pattern = os.path.join(self.current_folder, f"*.{file_type}")
            files.extend(glob.glob(pattern))

        return sorted(files)

class AdvancedProgressDialog:
    """Enhanced progress dialog with detailed logging and statistics"""
    def __init__(self, parent, title="Processing", show_details=True):
        self.top = tk.Toplevel(parent)
        self.top.title(title)
        self.top.geometry("600x400")
        # Make dialog independent (not modal)
        # self.top.transient(parent)  # Commented out to prevent window following
        # self.top.grab_set()  # Commented out to make non-modal

        # Center the dialog
        self.center_dialog()

        # Progress variables
        self.progress_var = tk.DoubleVar()
        self.current_task = tk.StringVar()
        self.total_tasks = 0
        self.completed_tasks = 0
        self.failed_tasks = 0
        self.start_time = datetime.now()

        # Create UI
        self.create_widgets(show_details)

        # Control flags
        self.cancelled = False
        self.completed = False
        self.paused = False

    def center_dialog(self):
        """Center dialog on screen"""
        self.top.update_idletasks()
        x = (self.top.winfo_screenwidth() // 2) - (600 // 2)
        y = (self.top.winfo_screenheight() // 2) - (400 // 2)
        self.top.geometry(f"600x400+{x}+{y}")

    def create_widgets(self, show_details):
        """Create dialog widgets"""
        # Main progress section
        progress_frame = ttk.Frame(self.top)
        progress_frame.pack(fill=tk.X, padx=20, pady=10)

        # Progress bar
        self.progress_bar = ttk.Progressbar(progress_frame, variable=self.progress_var, maximum=100)
        self.progress_bar.pack(fill=tk.X, pady=(0, 10))

        # Current task display
        task_frame = ttk.Frame(progress_frame)
        task_frame.pack(fill=tk.X, pady=(0, 10))

        ttk.Label(task_frame, text="Current Task:", font=('TkDefaultFont', 9, 'bold')).pack(anchor=tk.W)
        self.task_label = ttk.Label(task_frame, textvariable=self.current_task, font=('TkDefaultFont', 9))
        self.task_label.pack(anchor=tk.W, padx=(20, 0))

        # Statistics frame
        stats_frame = ttk.LabelFrame(progress_frame, text="📊 Statistics", padding=5)
        stats_frame.pack(fill=tk.X, pady=(0, 10))

        self.stats_label = ttk.Label(stats_frame, text="Starting...", font=('TkDefaultFont', 9))
        self.stats_label.pack(anchor=tk.W)

        if show_details:
            # Detailed log section
            log_frame = ttk.LabelFrame(self.top, text="📋 Detailed Log", padding=5)
            log_frame.pack(fill=tk.BOTH, expand=True, padx=20, pady=(0, 10))

            # Log text with scrollbar
            log_container = ttk.Frame(log_frame)
            log_container.pack(fill=tk.BOTH, expand=True)

            self.log_text = tk.Text(log_container, height=10, font=('Courier', 8), wrap=tk.WORD)
            log_scroll = ttk.Scrollbar(log_container, orient="vertical", command=self.log_text.yview)

            self.log_text.config(yscrollcommand=log_scroll.set)
            self.log_text.grid(row=0, column=0, sticky='nsew')
            log_scroll.grid(row=0, column=1, sticky='ns')

            log_container.grid_rowconfigure(0, weight=1)
            log_container.grid_columnconfigure(0, weight=1)
        else:
            self.log_text = None

        # Control buttons
        button_frame = ttk.Frame(self.top)
        button_frame.pack(fill=tk.X, padx=20, pady=10)

        self.pause_button = ttk.Button(button_frame, text="⏸️ Pause", command=self.toggle_pause)
        self.pause_button.pack(side=tk.LEFT, padx=(0, 5))

        self.cancel_button = ttk.Button(button_frame, text="❌ Cancel", command=self.cancel)
        self.cancel_button.pack(side=tk.LEFT)

        self.close_button = ttk.Button(button_frame, text="✅ Close", command=self.close, state='disabled')
        self.close_button.pack(side=tk.RIGHT)

        self.details_button = ttk.Button(button_frame, text="💾 Save Log", command=self.save_log, state='disabled')
        self.details_button.pack(side=tk.RIGHT, padx=(0, 5))

    def update_progress(self, value, task="", log_message="", failed=False):
        """Update progress dialog with enhanced statistics - THREAD SAFE"""
        def _do_gui_update():
            self.progress_var.set(value)

            if task:
                self.current_task.set(task)

            # Update statistics
            elapsed = datetime.now() - self.start_time
            elapsed_str = str(elapsed).split('.')[0]  # Remove microseconds

            if value > 0:
                eta_total = elapsed.total_seconds() * (100 / value)
                eta_remaining = eta_total - elapsed.total_seconds()
                eta_str = str(timedelta(seconds=int(eta_remaining)))
            else:
                eta_str = "Calculating..."

            stats_text = (
                f"⏱️  Elapsed: {elapsed_str} | 🔮 ETA: {eta_str}\n"
                f"✅ Completed: {self.completed_tasks} | ❌ Failed: {self.failed_tasks}"
            )
            self.stats_label.config(text=stats_text)

            # Update log
            if log_message and self.log_text:
                timestamp = datetime.now().strftime('%H:%M:%S')
                status_icon = "❌" if failed else "✅"
                log_entry = f"[{timestamp}] {status_icon} {log_message}\n"

                self.log_text.insert(tk.END, log_entry)
                self.log_text.see(tk.END)

            self.top.update_idletasks()

        # Thread-safe GUI update
        self.top.after(0, _do_gui_update)

    def complete(self, message="Processing completed"):
        """Mark processing as completed - THREAD SAFE"""
        def _do_completion():
            self.completed = True
            self.progress_var.set(100)
            self.current_task.set(message)

            # Final statistics
            total_time = datetime.now() - self.start_time
            total_time_str = str(total_time).split('.')[0]

            final_stats = (
                f"🏁 Completed in: {total_time_str}\n"
                f"✅ Successful: {self.completed_tasks} | ❌ Failed: {self.failed_tasks}"
            )
            self.stats_label.config(text=final_stats)

            # Enable controls
            self.cancel_button.config(state='disabled')
            self.pause_button.config(state='disabled')
            self.close_button.config(state='normal')

            if self.log_text:
                self.details_button.config(state='normal')

        # Thread-safe GUI update
        self.top.after(0, _do_completion)

    def toggle_pause(self):
        """Toggle pause state"""
        self.paused = not self.paused
        if self.paused:
            self.pause_button.config(text="▶️ Resume")
            self.current_task.set("⏸️ Processing paused...")
        else:
            self.pause_button.config(text="⏸️ Pause")

    def cancel(self):
        """Cancel processing"""
        self.cancelled = True
        self.close()

    def close(self):
        """Close dialog"""
        self.top.destroy()

    def save_log(self):
        """Save detailed log to file"""
        if not self.log_text:
            return

        filename = filedialog.asksaveasfilename(
            title="Save Processing Log",
            defaultextension=".txt",
            filetypes=[("Text files", "*.txt"), ("All files", "*.*")]
        )

        if filename:
            try:
                with open(filename, 'w', encoding='utf-8') as f:
                    f.write(self.log_text.get(1.0, tk.END))
                messagebox.showinfo("Success", f"Log saved to:\n{filename}")
            except Exception as e:
                messagebox.showerror("Error", f"Failed to save log:\n{str(e)}")

class StatisticsPanel(ttk.LabelFrame):
    """Statistics display and analysis panel"""
    def __init__(self, parent, title="📊 Statistics", **kwargs):
        super().__init__(parent, text=title, padding=10, **kwargs)

        self.stats = {}
        self.create_widgets()

    def create_widgets(self):
        """Create statistics display widgets"""
        # Main statistics display
        self.stats_text = tk.Text(self, height=8, width=50, font=('Courier', 9),
                                 state='disabled', wrap=tk.WORD)
        stats_scroll = ttk.Scrollbar(self, orient="vertical", command=self.stats_text.yview)

        self.stats_text.config(yscrollcommand=stats_scroll.set)
        self.stats_text.grid(row=0, column=0, sticky='nsew', padx=(0, 5))
        stats_scroll.grid(row=0, column=1, sticky='ns')

        # Control buttons
        button_frame = ttk.Frame(self)
        button_frame.grid(row=1, column=0, columnspan=2, sticky='ew', pady=(10, 0))

        ttk.Button(button_frame, text="🔄 Refresh", command=self.refresh_display).pack(side=tk.LEFT, padx=(0, 5))
        ttk.Button(button_frame, text="💾 Export", command=self.export_stats).pack(side=tk.LEFT, padx=(0, 5))
        ttk.Button(button_frame, text="🗑️ Clear", command=self.clear_stats).pack(side=tk.LEFT)

        self.grid_rowconfigure(0, weight=1)
        self.grid_columnconfigure(0, weight=1)

    def update_stats(self, new_stats):
        """Update statistics data"""
        self.stats.update(new_stats)
        self.refresh_display()

    def refresh_display(self):
        """Refresh statistics display"""
        self.stats_text.config(state='normal')
        self.stats_text.delete(1.0, tk.END)

        if self.stats:
            display_text = self.format_stats()
            self.stats_text.insert(1.0, display_text)
        else:
            self.stats_text.insert(1.0, "No statistics available")

        self.stats_text.config(state='disabled')

    def format_stats(self):
        """Format statistics for display"""
        formatted = "📈 PROCESSING STATISTICS\n" + "="*40 + "\n\n"

        for category, data in self.stats.items():
            formatted += f"📊 {category.upper()}:\n"
            if isinstance(data, dict):
                for key, value in data.items():
                    if isinstance(value, float):
                        formatted += f"  • {key}: {value:.3f}\n"
                    else:
                        formatted += f"  • {key}: {value}\n"
            else:
                formatted += f"  • Value: {data}\n"
            formatted += "\n"

        return formatted

    def export_stats(self):
        """Export statistics to file"""
        if not self.stats:
            messagebox.showwarning("No Data", "No statistics to export")
            return

        filename = filedialog.asksaveasfilename(
            title="Export Statistics",
            defaultextension=".json",
            filetypes=[("JSON files", "*.json"), ("Text files", "*.txt")]
        )

        if filename:
            try:
                if filename.endswith('.json'):
                    with open(filename, 'w', encoding='utf-8') as f:
                        json.dump(self.stats, f, indent=2, default=str)
                else:
                    with open(filename, 'w', encoding='utf-8') as f:
                        f.write(self.format_stats())

                messagebox.showinfo("Success", f"Statistics exported to:\n{filename}")
            except Exception as e:
                messagebox.showerror("Error", f"Export failed:\n{str(e)}")

    def clear_stats(self):
        """Clear all statistics"""
        if messagebox.askyesno("Clear Statistics", "Clear all statistics data?"):
            self.stats.clear()
            self.refresh_display()

class ModeSelectionFrame(ttk.LabelFrame):
    """Processing mode selection frame with radio buttons"""
    def __init__(self, parent, title="🔧 Processing Mode", **kwargs):
        super().__init__(parent, text=title, padding=10, **kwargs)

        self.mode_var = tk.StringVar(value='full_detection')
        self.callback = None

        self.create_widgets()

    def create_widgets(self):
        """Create mode selection widgets"""
        # Mode description
        desc_frame = ttk.Frame(self)
        desc_frame.pack(fill=tk.X, pady=(0, 10))

        ttk.Label(desc_frame, text="Select the peak processing approach:",
                 font=('TkDefaultFont', 9, 'bold')).pack(anchor=tk.W)

        # Radio buttons
        radio_frame = ttk.Frame(self)
        radio_frame.pack(fill=tk.X, pady=(0, 10))

        # Full Detection mode
        self.full_radio = ttk.Radiobutton(
            radio_frame,
            text="🔍 Full Peak Detection",
            variable=self.mode_var,
            value='full_detection',
            command=self.on_mode_change
        )
        self.full_radio.pack(anchor=tk.W, pady=2)

        # Full mode description
        full_desc = ttk.Label(
            radio_frame,
            text="  Complete peak detection across entire spectrum with validation",
            font=('TkDefaultFont', 8),
            foreground='gray'
        )
        full_desc.pack(anchor=tk.W, padx=(20, 0))

        # In-Place mode
        self.inplace_radio = ttk.Radiobutton(
            radio_frame,
            text="🎯 In-Place Fitting",
            variable=self.mode_var,
            value='in_place',
            command=self.on_mode_change
        )
        self.inplace_radio.pack(anchor=tk.W, pady=(10, 2))

        # In-place mode description
        inplace_desc = ttk.Label(
            radio_frame,
            text="  Direct fitting to reference peak positions with advanced analysis",
            font=('TkDefaultFont', 8),
            foreground='gray'
        )
        inplace_desc.pack(anchor=tk.W, padx=(20, 0))

        # Status display
        self.status_frame = ttk.Frame(self)
        self.status_frame.pack(fill=tk.X, pady=(10, 0))

        self.status_label = ttk.Label(
            self.status_frame,
            text="📋 Mode: Full Peak Detection",
            font=('TkDefaultFont', 9, 'bold'),
            foreground='blue'
        )
        self.status_label.pack(anchor=tk.W)

    def on_mode_change(self):
        """Handle mode change"""
        mode = self.mode_var.get()

        if mode == 'full_detection':
            self.status_label.config(text="📋 Mode: Full Peak Detection")
        else:
            self.status_label.config(text="📋 Mode: In-Place Fitting")

        if self.callback:
            self.callback(mode)

    def set_callback(self, callback):
        """Set callback for mode changes"""
        self.callback = callback

    def get_mode(self):
        """Get current mode"""
        return self.mode_var.get()

    def set_mode(self, mode):
        """Set processing mode programmatically"""
        if mode in ['full_detection', 'in_place']:
            self.mode_var.set(mode)
            self.on_mode_change()


class PeakNavigator(ttk.Frame):
    """Peak Navigator panel for interactive peak list display and selection"""

    def __init__(self, parent, *args, **kwargs):
        super().__init__(parent, *args, **kwargs)

        # State variables
        self.spectrum_controller = None
        self.selected_peak_index = None
        self.selected_peak_type = "reference"
        self.reference_peaks = []
        self.detected_peaks = []

        # Setup UI
        self.setup_ui()

    def setup_ui(self):
        """Create the peak navigator UI"""

        # Header frame
        header_frame = ttk.Frame(self)
        header_frame.pack(fill=tk.X, pady=(0, 10))

        # Title
        title_label = ttk.Label(header_frame, text="🧭 Peak Navigator",
                               font=("Arial", 12, "bold"))
        title_label.pack(anchor="w")

        # Peak type selector
        selector_frame = ttk.Frame(header_frame)
        selector_frame.pack(fill=tk.X, pady=(5, 0))

        ttk.Label(selector_frame, text="Peak List:").pack(side=tk.LEFT)

        self.peak_type_var = tk.StringVar(value="Reference Peaks")
        self.type_combo = ttk.Combobox(selector_frame, textvariable=self.peak_type_var,
                                      values=["Reference Peaks", "Detected Peaks"],
                                      state="disabled", width=15) #readonly to "disabled"
        self.type_combo.pack(side=tk.LEFT, padx=(5, 0))
        self.type_combo.bind("<<ComboboxSelected>>", self.on_peak_type_changed)

        # Peak table frame
        table_frame = ttk.Frame(self)
        table_frame.pack(fill=tk.BOTH, expand=True, pady=(0, 10))

        # Create treeview with scrollbar
        columns = ("assignment", "x_coord", "y_coord")
        self.tree = ttk.Treeview(table_frame, columns=columns, show="headings", height=20)

        # Configure columns
        self.tree.heading("assignment", text="Assignment", anchor="center")
        self.tree.heading("x_coord", text="X (1H)", anchor="center")
        self.tree.heading("y_coord", text="Y (15N/13C)", anchor="center")

        self.tree.column("assignment", width=90, anchor="center", minwidth=70)
        self.tree.column("x_coord", width=90, anchor="center", minwidth=70)
        self.tree.column("y_coord", width=90, anchor="center", minwidth=70)

        # Scrollbar for table
        tree_scrollbar = ttk.Scrollbar(table_frame, orient="vertical", command=self.tree.yview)
        self.tree.configure(yscrollcommand=tree_scrollbar.set)

        # Pack table and scrollbar
        self.tree.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        tree_scrollbar.pack(side=tk.RIGHT, fill=tk.Y)

        # Bind selection event
        self.tree.bind("<<TreeviewSelect>>", self.on_peak_selected)

        # Status frame
        status_frame = ttk.Frame(self)
        status_frame.pack(fill=tk.X, pady=(0, 5))

        self.status_label = ttk.Label(status_frame, text="No peaks loaded",
                                     font=("Arial", 9), foreground="gray")
        self.status_label.pack(anchor="w")

        # Button frame for refresh and analysis
        button_frame = ttk.Frame(self)
        button_frame.pack(fill=tk.X)

        self.refresh_btn = ttk.Button(button_frame, text="🔄 Refresh",
                                     command=self.refresh_peak_list, width=12)
        self.refresh_btn.pack(side=tk.LEFT, padx=(0, 5))

        self.analysis_btn = ttk.Button(button_frame, text="🔬 Analyze",
                                      command=self.analyze_selected_peak, width=12)
        self.analysis_btn.pack(side=tk.LEFT)

    def set_spectrum_controller(self, controller):
        """Set reference to main GUI for spectrum control"""
        self.spectrum_controller = controller

    def on_peak_type_changed(self, event=None):
        """Handle peak type dropdown selection change"""
        peak_type = self.peak_type_var.get()
        print(f"Peak Navigator: Dropdown changed to: {peak_type}")
        self.selected_peak_type = "reference" if "Reference" in peak_type else "detected"
        print(f"Peak Navigator: Set selected_peak_type to: {self.selected_peak_type}")
        self.refresh_peak_list()

    def on_peak_selected(self, event):
        """Handle peak selection from table"""
        selection = self.tree.selection()
        if not selection:
            return

        # Get selected item data
        item = self.tree.item(selection[0])
        values = item['values']

        if len(values) >= 3:
            try:
                # Extract coordinates
                peak_x = float(values[1])
                peak_y = float(values[2])

                # Get peak index from table position (simpler and more reliable)
                children = self.tree.get_children()
                peak_index = children.index(selection[0])  # Get actual row index

                print(f"Peak Navigator: Selected peak {peak_index}: {values[0]} ({peak_x:.3f}, {peak_y:.1f})")

                # Update selected peak
                self.selected_peak_index = peak_index

                # Center spectrum on peak
                self.center_spectrum_on_peak(peak_x, peak_y)

                # Update status
                assignment = values[0] if values[0] else f"Peak {peak_index + 1}"
                self.status_label.config(text=f"Selected: {assignment} ({peak_x:.3f}, {peak_y:.1f})")

                # Notify spectrum controller for coordination
                if hasattr(self.spectrum_controller, 'set_selected_peak'):
                    self.spectrum_controller.set_selected_peak(peak_index, self.selected_peak_type, source="navigator")

            except (ValueError, IndexError) as e:
                print(f"Error processing peak selection: {e}")

    def center_spectrum_on_peak(self, x, y):
        """Center main spectrum plot on specified coordinates"""
        if self.spectrum_controller and hasattr(self.spectrum_controller, 'center_on_coordinates'):
            self.spectrum_controller.center_on_coordinates(x, y)
        elif self.spectrum_controller and hasattr(self.spectrum_controller, 'center_on_selected_peak'):
            # Fallback to existing center method
            self.spectrum_controller.center_on_selected_peak()

    def load_reference_peaks(self, peak_data):
        """Load reference peaks data"""
        self.reference_peaks = []

        if peak_data is None:
            print("Peak Navigator: No reference peak data provided")
            return

        try:
            # Handle pandas DataFrame (most common format from integrator.peak_list)
            if hasattr(peak_data, 'iloc') and hasattr(peak_data, 'columns'):
                for i, row in peak_data.iterrows():
                    assignment = row.get('Assignment', f"Ref{i+1}")
                    x_coord = float(row.get('Position_X', 0))
                    y_coord = float(row.get('Position_Y', 0))
                    self.reference_peaks.append([assignment, x_coord, y_coord])

            # Handle numpy array format
            elif hasattr(peak_data, 'shape') and len(peak_data.shape) == 2:
                for i, peak in enumerate(peak_data):
                    if len(peak) >= 2:
                        assignment = f"Ref{i+1}"  # Default assignment
                        x_coord = float(peak[0])
                        y_coord = float(peak[1])
                        self.reference_peaks.append([assignment, x_coord, y_coord])

            # Handle list/tuple format
            elif isinstance(peak_data, (list, tuple)):
                for i, peak in enumerate(peak_data):
                    if len(peak) >= 2:
                        assignment = f"Ref{i+1}"
                        x_coord = float(peak[0])
                        y_coord = float(peak[1])
                        self.reference_peaks.append([assignment, x_coord, y_coord])
            else:
                print(f"Peak Navigator Warning: Unrecognized reference peak data format: {type(peak_data)}")
                return

        except Exception as e:
            print(f"Peak Navigator: Error loading reference peaks: {e}")
            import traceback
            traceback.print_exc()

        # Refresh display if currently showing reference peaks
        if self.selected_peak_type == "reference":
            self.refresh_peak_list()

    def load_detected_peaks(self, fitted_peaks):
        """Load detected peaks data"""
        self.detected_peaks = []

        if fitted_peaks is None:
            print("Peak Navigator: No detected peak data provided")
            return

        try:
            # Handle list of dictionaries (most common format from integrator.fitted_peaks)
            if isinstance(fitted_peaks, list) and fitted_peaks and isinstance(fitted_peaks[0], dict):
                for i, peak in enumerate(fitted_peaks):
                    # Try different possible coordinate keys
                    x_coord = peak.get('ppm_x', peak.get('Position_X', peak.get('position_x', 0)))
                    y_coord = peak.get('ppm_y', peak.get('Position_Y', peak.get('position_y', 0)))
                    assignment = peak.get('assignment', peak.get('Assignment', f"Det{i+1}"))

                    if x_coord != 0 and y_coord != 0:  # Skip empty coordinates
                        self.detected_peaks.append([assignment, float(x_coord), float(y_coord)])

            # Handle pandas DataFrame
            elif hasattr(fitted_peaks, 'iloc') and hasattr(fitted_peaks, 'columns'):
                for i, row in fitted_peaks.iterrows():
                    assignment = row.get('Assignment', row.get('assignment', f"Det{i+1}"))
                    x_coord = float(row.get('Position_X', row.get('ppm_x', 0)))
                    y_coord = float(row.get('Position_Y', row.get('ppm_y', 0)))
                    if x_coord != 0 and y_coord != 0:
                        self.detected_peaks.append([assignment, x_coord, y_coord])

            # Handle numpy array or list of lists
            elif isinstance(fitted_peaks, (list, tuple)) and fitted_peaks:
                for i, peak in enumerate(fitted_peaks):
                    if hasattr(peak, '__len__') and len(peak) >= 2:
                        assignment = f"Det{i+1}"
                        x_coord = float(peak[0])
                        y_coord = float(peak[1])
                        self.detected_peaks.append([assignment, x_coord, y_coord])
            else:
                print(f"Peak Navigator Warning: Unrecognized detected peak data format: {type(fitted_peaks)}")
                if fitted_peaks:
                    print(f"Sample data: {fitted_peaks[0] if len(fitted_peaks) > 0 else 'Empty'}")
                return

        except Exception as e:
            print(f"Peak Navigator: Error loading detected peaks: {e}")
            import traceback
            traceback.print_exc()

        # Refresh display if currently showing detected peaks
        if self.selected_peak_type == "detected":
            self.refresh_peak_list()

    def refresh_peak_list(self):
        """Refresh the peak list display"""
        print(f"Peak Navigator: refresh_peak_list() called")
        print(f"Peak Navigator: selected_peak_type = {self.selected_peak_type}")

        # Clear existing items
        for item in self.tree.get_children():
            self.tree.delete(item)

        # Get current peak list
        if self.selected_peak_type == "reference":
            peaks = self.reference_peaks
            peak_type_name = "Reference Peaks"
            print(f"Peak Navigator: Using reference peaks: {len(peaks)} peaks")
        else:
            peaks = self.detected_peaks
            peak_type_name = "Detected Peaks"
            print(f"Peak Navigator: Using detected peaks: {len(peaks)} peaks")

        # Populate table
        for i, peak in enumerate(peaks):
            if len(peak) >= 3:
                assignment, x_coord, y_coord = peak[0], peak[1], peak[2]
                self.tree.insert("", "end", values=(assignment, f"{x_coord:.3f}", f"{y_coord:.1f}"))
                print(f"Peak Navigator: Added peak {i+1}: {assignment} ({x_coord:.3f}, {y_coord:.1f})")

        # Update status
        count = len(peaks)
        self.status_label.config(text=f"{peak_type_name}: {count} peak{'s' if count != 1 else ''}")
        print(f"Peak Navigator: Updated status: {count} peaks")

    def update_selection(self, peak_index, peak_type):
        """Update navigator selection from external source"""
        if peak_type != self.selected_peak_type:
            # Switch peak type if needed
            type_name = "Reference Peaks" if peak_type == "reference" else "Detected Peaks"
            self.peak_type_var.set(type_name)
            self.on_peak_type_changed()

        # Select the peak in the table
        try:
            children = self.tree.get_children()
            if 0 <= peak_index < len(children):
                # Clear previous selection
                self.tree.selection_remove(self.tree.selection())
                # Select new peak
                self.tree.selection_add(children[peak_index])
                self.tree.focus(children[peak_index])
                self.tree.see(children[peak_index])

        except Exception as e:
            print(f"Error updating navigator selection: {e}")

    def analyze_selected_peak(self):
        """Show analysis for selected peak - same function as navigation panel button"""
        from tkinter import messagebox

        selection = self.tree.selection()
        if not selection:
            messagebox.showinfo("No Selection", "Please select a peak from the table first.")
            return

        try:
            children = self.tree.get_children()
            peak_index = children.index(selection[0])

            if hasattr(self, 'spectrum_controller') and self.spectrum_controller:
                self.spectrum_controller.navigator_show_peak_analysis(
                    self.selected_peak_type, peak_index
                )
            else:
                messagebox.showerror("Error", "No connection to main GUI controller.")

        except Exception as e:
            messagebox.showerror("Error", f"Failed to analyze peak: {str(e)}")
            print(f"Peak Navigator: Analysis error: {e}")
            import traceback
            traceback.print_exc()

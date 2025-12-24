#!/usr/bin/env python3
"""
Standalone script for plotting dual-field NMR relaxation data - J(0.87ωH) Version
from ZZ_multi_2fields_density_087.py output files.

J(0.87ωH) VERSION: Updated to plot J(0.87ωH) instead of J(ωH) for more accurate
spectral density visualization in NMR relaxation analysis.

User can select which columns to plot and which field to display.
Multiple datasets can be plotted side by side for comparison.

MISSING RESIDUE HANDLING (RES VERSION):
- Missing residues are shown as grey bars from 0 to max value (preserving sign)
- X-axis spans complete residue range from residue 1 to max residue numbers
- Ensures visual continuity and highlights data gaps

Y-AXIS LIMITS CONFIGURATION (SS VERSION):
- Supports JSON configuration files for custom y-axis limits
- Backward compatible: uses matplotlib auto-scaling if no config provided
- Supports parameter-specific and field-specific limits
- Example config: {"ylimits": {"R1": {"min": 0, "max": 5}}}

SECONDARY STRUCTURE VISUALIZATION (SS VERSION):
- Displays protein secondary structure at the top of plots
- Requires text file with H (helix), B (beta sheet), L (loop) symbols
- One symbol per residue, matching the length of the protein sequence
- Example: HHHHHLLLLBBBBBLLnWTLLHHHH (for 22 residues)
"""

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import argparse
import sys
import os
import json
import matplotlib.patches as mpatches
import matplotlib.path as mpath
import colorsys

# Configure matplotlib for Adobe Illustrator compatibility
# Keep text as editable text objects rather than converting to paths
matplotlib.rcParams['pdf.fonttype'] = 42  # TrueType fonts
matplotlib.rcParams['ps.fonttype'] = 42   # PostScript compatibility
matplotlib.rcParams['font.family'] = 'sans-serif'
# Use fonts that support Unicode subscripts/superscripts
matplotlib.rcParams['font.sans-serif'] = ['DejaVu Sans', 'Liberation Sans', 'Bitstream Vera Sans', 'Arial Unicode MS', 'Arial']

def adjust_color(hex_color, lighten_amount=0.8, desaturate_amount=0.3):
    """
    Adjust the given hex color by lightening and/or desaturating.

    Parameters:
    - lighten_amount: float (0 to 1), amount to lighten the color
    - desaturate_amount: float (0 to 1), amount to desaturate the color

    Returns:
    - Adjusted hex color string.
    """
    hex_color = hex_color.lstrip('#')
    r, g, b = tuple(int(hex_color[i:i+2], 16) / 255.0 for i in (0, 2, 4))

    # Convert RGB to HLS
    h, l, s = colorsys.rgb_to_hls(r, g, b)

    # Adjust lightness
    l = min(1.0, l + lighten_amount * (1.0 - l))

    # Adjust saturation (reduce toward 0)
    s = max(0.0, s * (1.0 - desaturate_amount))

    # Convert back to RGB
    r, g, b = colorsys.hls_to_rgb(h, l, s)
    return '#{0:02x}{1:02x}{2:02x}'.format(int(r * 255), int(g * 255), int(b * 255))

def load_secondary_structure(ss_file):
    """
    Load secondary structure from text file
    
    Parameters:
    -----------
    ss_file : str
        Path to secondary structure file with H/B/L symbols
        
    Returns:
    --------
    list : List of secondary structure elements with spans
    """
    if not ss_file or not os.path.exists(ss_file):
        return None
    
    try:
        with open(ss_file, 'r') as f:
            ss_string = f.read().strip()
        
        print(f"Loaded secondary structure from: {ss_file}")
        print(f"Secondary structure length: {len(ss_string)}")
        
        # Convert string to spans
        ss_map = []
        current_type = None
        current_start = None
        
        for i, ss_char in enumerate(ss_string):
            ss_type = None
            if ss_char.upper() == 'H':
                ss_type = 'helix'
            elif ss_char.upper() == 'B':
                ss_type = 'sheet'
            elif ss_char.upper() == 'L':
                ss_type = 'loop'
            
            if ss_type != current_type:
                # Close previous span if any
                if current_type is not None:
                    ss_map.append({
                        "type": current_type,
                        "span": (current_start + 1, i)  # Convert to 1-based indexing
                    })
                # Start new span
                current_type = ss_type
                current_start = i
        
        # Append last span
        if current_type is not None:
            ss_map.append({
                "type": current_type,
                "span": (current_start + 1, len(ss_string))  # Convert to 1-based indexing
            })
        
        return ss_map
        
    except IOError as e:
        print(f"Warning: Could not load secondary structure file {ss_file}: {e}")
        return None

def draw_secondary_structure(target_ax, range, secondary_structure_map: list) -> None:
    """
    Draw secondary structure representation (adapted from ss_plot.py)
    """
    def draw_rect(ax, girth:float, range: tuple, color: str=None) -> None:
        range = (range[0]-0.5, range[1]+0.5)
        span = range[1] - range[0] # span
        bottom_left_coords = (range[0], -girth/2) # start x, y
        ax.add_artist(mpatches.Rectangle(bottom_left_coords, span, girth, color=color, clip_on=False))
        return

    def draw_helix(ax, r: tuple, h:float, annot:str ="", strand_width:float = 1.0, angle:float = 10, c: list=["red", "white"], ec: str="black", lw: float=0.5) -> None:
        r = (r[0]-0.5, r[1]+0.5)
        w = r[1] - r[0]
        coords = (r[0], -h/2)

        w_para = strand_width
        l = h / np.sin(np.deg2rad(angle))
        d = l * np.cos(np.deg2rad(angle))
        max_iter = int(w)
        offset = 0.1

        start_x = r[0] + d + offset
        i = 0
        while start_x + w_para + d < r[1] and i < max_iter:
            codes, verts = zip(*[
                (mpath.Path.MOVETO, [start_x, h/2]),
                (mpath.Path.LINETO, [start_x + d, -h/2]),
                (mpath.Path.LINETO, [start_x + d + w_para, -h/2]),
                (mpath.Path.LINETO, [start_x + w_para, h/2]),
                (mpath.Path.CLOSEPOLY, [start_x, h/2])
            ])
            ax.add_artist(mpatches.PathPatch(mpath.Path(verts, codes), facecolor=c[1], lw=lw, ec=ec, clip_on=True))
            start_x += d*2 + offset*2
            i += 1

        start_x = r[0]
        i = 0
        while start_x + w_para + d < r[1] and i < max_iter:
            codes, verts = zip(*[
                (mpath.Path.MOVETO, [start_x, -h/2]),
                (mpath.Path.LINETO, [start_x + d, h/2]),
                (mpath.Path.LINETO, [start_x + d + w_para, h/2]),
                (mpath.Path.LINETO, [start_x + w_para, -h/2]),
                (mpath.Path.CLOSEPOLY, [start_x, -h/2])
            ])
            ax.add_artist(mpatches.PathPatch(mpath.Path(verts, codes), facecolor=c[0], lw=lw, ec=ec, clip_on=True))
            start_x += d*2 + offset*2
            i += 1

        return

    def draw_sheet(ax, r: tuple, h:float, annot:str="", a_w: float=2, a_h:float=1, c: str="dodgerblue", ec: str= "black", lw: float = 0.5) -> None:
        r = (r[0]-0.5, r[1]+0.5)
        w = r[1] - r[0]
        w_arrow_head = max(a_w, w*0.1)
        h_arrow_head = a_h
        w_rec = r[1] - r[0] - w_arrow_head

        codes, verts = zip(*[
            (mpath.Path.MOVETO, [r[0], h/2]),
            (mpath.Path.LINETO, [r[0] + w_rec, h/2]),
            (mpath.Path.LINETO, [r[0] + w_rec,  h_arrow_head/2]),
            (mpath.Path.LINETO, [r[1],  0]),
            (mpath.Path.LINETO, [r[0] + w_rec,  -h_arrow_head/2]),
            (mpath.Path.LINETO, [r[0] + w_rec,  -h/2]),
            (mpath.Path.LINETO, [r[0],  -h/2]),
            (mpath.Path.CLOSEPOLY, [r[0], h/2])
        ])
        ax.add_artist(mpatches.PathPatch(mpath.Path(verts, codes), facecolor=c, lw=lw, ec=ec, clip_on=True))

        return

    if not secondary_structure_map:
        return

    target_ax.set_xlim(range[0] - 0.5, range[1] - 0.5)
    target_ax.set_ylim([-0.6, 0.6])
    target_ax.set_yticks([])
    target_ax.set_yticklabels([])
    target_ax.set_xticks([])
    target_ax.set_xticklabels([])
    target_ax.set_axis_off()

    lw = 1
    helix_width = 1.5
    helix_height = 1
    helix_angle = 35
    helix_color_primary = "#dc143c"
    helix_color_secondary = "#ffc0cb"
    helix_edge_color = "black"
    sheet_height = 0.65
    sheet_arrow_width = 2
    sheet_arrow_height = 1
    sheet_color = "#1e90ff"
    sheet_edge_color = "black"
    
    # Draw backbone
    draw_rect(target_ax, girth=0.05, range=target_ax.get_xlim(), color="gray")

    for secondary_structure in secondary_structure_map:
        if secondary_structure["type"] == "sheet":
            draw_sheet(
                target_ax,
                secondary_structure["span"],
                sheet_height,
                annot="",
                a_w=sheet_arrow_width,
                a_h=sheet_arrow_height,
                c=sheet_color,
                ec=sheet_edge_color,
                lw=lw
            )
        elif secondary_structure["type"] == "helix":
            helix_colors = [helix_color_primary, helix_color_secondary]
            draw_helix(
                target_ax,
                secondary_structure["span"],
                helix_height,
                annot="",
                strand_width=helix_width,
                angle=helix_angle,
                c=helix_colors,
                lw=lw,
                ec=helix_edge_color
            )
    return

def load_ylimit_config(config_file="ylimits_config.json"):
    """
    Load y-limit configuration with fallback to auto-scaling
    
    Parameters:
    -----------
    config_file : str
        Path to JSON configuration file
        
    Returns:
    --------
    dict or None : Configuration dictionary or None if not found/invalid
    """
    if not config_file or not os.path.exists(config_file):
        return None
    
    try:
        with open(config_file, 'r') as f:
            config = json.load(f)
            print(f"Loaded y-limits configuration from: {config_file}")
            return config
    except (json.JSONDecodeError, IOError) as e:
        print(f"Warning: Could not load {config_file}: {e}")
        return None

def extract_base_column_name(column):
    """
    Extract base column name (remove _f1, _f2, _err suffixes)
    
    Parameters:
    -----------
    column : str
        Column name potentially with suffixes
        
    Returns:
    --------
    str : Base column name
    """
    base = column.replace('_f1', '').replace('_f2', '').replace('_err', '')
    return base

def apply_ylimits(ax, column, ylimit_config):
    """
    Apply y-limits from config with fallback to auto-scaling
    
    Parameters:
    -----------
    ax : matplotlib axis
        Plot axis to modify
    column : str
        Column name being plotted
    ylimit_config : dict or None
        Configuration dictionary
    """
    if not ylimit_config or "ylimits" not in ylimit_config:
        return  # Use matplotlib auto-scaling
    
    # Extract base column name for lookup
    base_column = extract_base_column_name(column)
    
    # Try exact match first, then base column
    limits = ylimit_config["ylimits"].get(column)  # Exact match (e.g., R1_f1)
    if not limits:
        limits = ylimit_config["ylimits"].get(base_column)  # Base match (e.g., R1)
    
    # Also check field_specific
    if not limits and "field_specific" in ylimit_config:
        limits = ylimit_config["field_specific"].get(column)
    
    if limits and isinstance(limits, dict) and "min" in limits and "max" in limits:
        try:
            ax.set_ylim(limits["min"], limits["max"])
            print(f"  Applied y-limits for {column}: [{limits['min']}, {limits['max']}]")
        except (ValueError, TypeError) as e:
            print(f"  Warning: Invalid y-limits for {column}: {e}")

def create_example_config(filename="example_ylimits_config.json"):
    """
    Create an example configuration file
    
    Parameters:
    -----------
    filename : str
        Output filename for example config
    """
    example_config = {
        "_comment": "Y-axis limits configuration for NMR data plots",
        "ylimits": {
            "R1": {"min": 0, "max": 5},
            "R2": {"min": 0, "max": 50},
            "hetNOE": {"min": 0, "max": 1.2},
            "J0": {"min": 0, "max": 5},
            "JwN": {"min": 0, "max": 2},
            "JwH_087": {"min": 0, "max": 0.5},
            "S2": {"min": 0, "max": 1},
            "te": {"min": 0, "max": 100},
            "Rex": {"min": 0, "max": 15},
            "tc": {"min": 0, "max": 50}
        },
        "field_specific": {
            "_comment": "Field-specific limits override general limits",
            "R1_f1": {"min": 0, "max": 4},
            "R1_f2": {"min": 0, "max": 6},
            "R2_f1": {"min": 0, "max": 40},
            "R2_f2": {"min": 0, "max": 60}
        }
    }
    
    try:
        with open(filename, 'w') as f:
            json.dump(example_config, f, indent=2)
        print(f"Example y-limits configuration created: {filename}")
        print("Edit this file to customize y-axis limits for your data.")
    except IOError as e:
        print(f"Error creating example config: {e}")

def create_complete_residue_data(df, data_column, error_column=None):
    """
    Create complete residue data with missing residues filled with NaN
    
    Parameters:
    -----------
    df : pd.DataFrame
        Input dataframe
    data_column : str
        Name of the data column
    error_column : str
        Name of the error column (optional)
        
    Returns:
    --------
    dict : Contains complete residue data
    """
    if 'Residue' not in df.columns or data_column not in df.columns:
        # Fallback to index if no Residue column
        residues = df.index.values
        data = df[data_column].values
        errors = df[error_column].values if error_column and error_column in df.columns else None
        return {'residues': residues, 'data': data, 'errors': errors}
    
    # Get actual residues and data
    actual_residues = df['Residue'].values
    actual_data = df[data_column].values
    actual_errors = df[error_column].values if error_column and error_column in df.columns else None
    
    # Create complete residue range starting from 1
    min_res = 1  # Always start from residue 1
    max_res = int(max(actual_residues))
    complete_residues = np.arange(min_res, max_res + 1)
    
    # Initialize complete data arrays with NaN
    complete_data = np.full(len(complete_residues), np.nan)
    complete_errors = np.full(len(complete_residues), np.nan) if actual_errors is not None else None
    
    # Fill in existing data
    for i, res in enumerate(actual_residues):
        idx = int(res) - min_res  # min_res is now always 1
        if 0 <= idx < len(complete_residues):
            complete_data[idx] = actual_data[i]
            if complete_errors is not None:
                complete_errors[idx] = actual_errors[i]
    
    return {
        'residues': complete_residues,
        'data': complete_data, 
        'errors': complete_errors,
        'max_abs_value': np.nanmax(np.abs(actual_data)) if len(actual_data) > 0 else 1.0
    }

def plot_missing_residues_bars(ax, residue_data, color='lightgrey', alpha=0.3):
    """
    Plot grey bars for missing residues from 0 to max value (preserving sign)
    
    Parameters:
    -----------
    ax : matplotlib axis
        Plot axis
    residue_data : dict
        Complete residue data from create_complete_residue_data
    color : str
        Color for missing residue bars
    alpha : float
        Transparency for missing residue bars
    """
    residues = residue_data['residues']
    data = residue_data['data']
    
    # Find missing residues (NaN values)
    missing_mask = np.isnan(data)
    missing_residues = residues[missing_mask]
    
    if len(missing_residues) > 0:
        # Find the actual max and min values (preserving signs)
        valid_data = data[~missing_mask]
        if len(valid_data) > 0:
            max_positive = np.nanmax(valid_data) if np.any(valid_data > 0) else 0
            min_negative = np.nanmin(valid_data) if np.any(valid_data < 0) else 0
            
            # Plot grey bars from 0 to max positive value if exists
            if max_positive > 0:
                ax.bar(missing_residues, [max_positive] * len(missing_residues), 
                       bottom=0, color=color, alpha=alpha, width=0.8, label='Missing residues')
            
            # Plot grey bars from 0 to min negative value if exists  
            if min_negative < 0:
                ax.bar(missing_residues, [min_negative] * len(missing_residues), 
                       bottom=0, color=color, alpha=alpha, width=0.8)

def get_available_columns(df):
    """Get available data columns for plotting"""
    # Define standard column groups
    relaxation_cols = ['R1_f1', 'R2_f1', 'hetNOE_f1', 'R1_f2', 'R2_f2', 'hetNOE_f2']
    spectral_density_cols = ['J0_f1', 'JwN_f1', 'JwH_087_f1', 'J0_f2', 'JwN_f2', 'JwH_087_f2']
    model_free_cols = ['S2', 'te', 'Rex_f1', 'Rex_f2']
    
    available = {}
    available['relaxation'] = [col for col in relaxation_cols if col in df.columns]
    available['spectral_density'] = [col for col in spectral_density_cols if col in df.columns]
    available['model_free'] = [col for col in model_free_cols if col in df.columns]
    
    return available

def plot_selected_data(datasets, column, field='f1', output_file=None, ylimit_config=None, ss_map=None):
    """
    Plot selected column data from one or more datasets
    
    Parameters:
    -----------
    datasets : list of tuples
        List of (dataframe, label) pairs
    column : str
        Column name to plot
    field : str
        Field identifier ('f1' or 'f2')
    output_file : str
        Output PDF filename
    ylimit_config : dict
        Y-axis limits configuration
    ss_map : list
        Secondary structure mapping
    """
    n_datasets = len(datasets)
    
    # Add space for secondary structure if provided
    if ss_map:
        # Create layout with SS panel at top
        height_ratios = [1, 4]  # SS panel + data panel
        fig, axes = plt.subplots(2, n_datasets, figsize=(8*n_datasets, 7),
                                gridspec_kw={'height_ratios': height_ratios})
        
        # Ensure axes is always 2D array
        if n_datasets == 1:
            axes = axes.reshape(-1, 1)
        
        ss_axes = axes[0, :]
        data_axes = axes[1, :]
        
        # Make ss_axes and data_axes lists if single dataset
        if n_datasets == 1:
            ss_axes = [ss_axes]
            data_axes = [data_axes]
    else:
        # Original layout without SS
        if n_datasets == 1:
            fig, ax = plt.subplots(1, 1, figsize=(8, 6))
            data_axes = [ax]
        else:
            fig, axes = plt.subplots(1, n_datasets, figsize=(8*n_datasets, 6))
            data_axes = list(axes) if n_datasets == 2 else axes
    
    # Plot each dataset
    for i, (df, label) in enumerate(datasets):
        ax = data_axes[i]
        
        # Get column name with field suffix if needed
        if field == 'none':
            # Field-independent columns (S2, tau_e, Rex, tc)
            plot_column = column
            error_column = f"{column}_err"
        elif field in ['f1', 'f2'] and not column.endswith(f'_{field}'):
            plot_column = f"{column}_{field}"
            error_column = f"{plot_column}_err"
        else:
            plot_column = column
            error_column = f"{column}_err"
        
        # Check if column exists
        if plot_column not in df.columns:
            print(f"Warning: Column '{plot_column}' not found in dataset '{label}'")
            continue
        
        # Get complete residue data (fills missing residues with NaN)
        residue_data = create_complete_residue_data(df, plot_column, error_column)
        residues = residue_data['residues']
        y_data = residue_data['data']
        y_err = residue_data['errors']
        
        # Draw secondary structure if provided
        if ss_map:
            residue_range = (residues.min(), residues.max())
            draw_secondary_structure(ss_axes[i], residue_range, ss_map)
            ss_axes[i].set_title(f'Secondary Structure - {label}')
        
        # Plot grey bars for missing residues first (so they appear behind data)
        plot_missing_residues_bars(ax, residue_data)
        
        # Plot data (only non-NaN values)
        valid_mask = ~np.isnan(y_data)
        if np.any(valid_mask):
            valid_residues = residues[valid_mask]
            valid_y_data = y_data[valid_mask]
            valid_y_err = y_err[valid_mask] if y_err is not None else None
            
            if valid_y_err is not None:
                ax.errorbar(valid_residues, valid_y_data, yerr=valid_y_err, fmt='o', capsize=3, 
                           alpha=0.7, label=label)
            else:
                ax.plot(valid_residues, valid_y_data, 'o', alpha=0.7, label=label)
        
        # Set labels and title
        ax.set_xlabel('Residue')
        ax.set_ylabel(get_ylabel(plot_column))
        ax.set_title(f'{get_plot_title(plot_column)} - {label}')
        ax.grid(True, alpha=0.3)
        
        # Apply y-limits from configuration
        apply_ylimits(ax, plot_column, ylimit_config)
    
    plt.tight_layout()
    
    if output_file:
        # Save with Adobe Illustrator compatible settings
        try:
            plt.savefig(output_file, dpi=300, bbox_inches='tight', fonttype=42)
        except TypeError:
            # Fallback if fonttype parameter is not supported
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"Plot saved to: {output_file}")
    else:
        plt.show()

def plot_multiple_columns(datasets, columns, field='f1', output_file=None, ylimit_config=None, ss_map=None):
    """
    Plot multiple columns from datasets in a single column layout (n rows x 1 column)
    
    Parameters:
    -----------
    datasets : list of tuples
        List of (dataframe, label) pairs
    columns : list
        List of column names to plot
    field : str
        Field identifier ('f1', 'f2', or 'auto')
    output_file : str
        Output PDF filename
    ylimit_config : dict
        Y-axis limits configuration
    ss_map : list
        Secondary structure mapping
    """
    n_plots = len(columns)
    n_datasets = len(datasets)
    
    # Add space for secondary structure if provided
    if ss_map:
        # Create layout with SS panel at top
        height_ratios = [1] + [4] * n_plots  # SS panel + data panels
        fig, axes = plt.subplots(n_plots + 1, 1, figsize=(8, 4*n_plots + 1),
                                gridspec_kw={'height_ratios': height_ratios})
        
        # Ensure axes is always a list
        if n_plots == 0:
            axes = [axes]
        
        # Draw secondary structure in top panel
        if len(datasets) > 0:
            # Use first dataset to determine residue range
            df_sample = datasets[0][0]
            if 'Residue' in df_sample.columns:
                residue_range = (df_sample['Residue'].min(), df_sample['Residue'].max())
            else:
                residue_range = (1, len(df_sample))
            
            draw_secondary_structure(axes[0], residue_range, ss_map)
            axes[0].set_title('Secondary Structure')
        
        # Adjust data panel indices
        data_axes = axes[1:]
    else:
        # Original layout without SS
        fig, axes = plt.subplots(n_plots, 1, figsize=(8, 4*n_plots))
        
        # Ensure axes is always a list
        if n_plots == 1:
            axes = [axes]
        
        data_axes = axes
    
    # Plot each column
    for i, column in enumerate(columns):
        ax = data_axes[i]
        
        # Determine field for this specific column
        field_independent_cols = ['S2', 'te', 'tc'] #'tau_e', 
        rex_cols = ['Rex']  # Rex has field-dependent versions
        
        if column in field_independent_cols:
            current_field = 'none'
        elif column in rex_cols:
            current_field = field  # Rex needs field specification
        elif any(f'_{f}' in column for f in ['f1', 'f2']):
            current_field = 'auto'
        else:
            current_field = field
        
        # Handle both fields plotting
        if current_field == 'both' and column not in field_independent_cols:
            # Plot both f1 and f2 for field-dependent columns
            field_colors = {'f1': 'blue', 'f2': 'red'}
            field_markers = {'f1': 'o', 'f2': 's'}
            
            for field_suffix in ['f1', 'f2']:
                for j, (df, label) in enumerate(datasets):
                    plot_column = f"{column}_{field_suffix}"
                    error_column = f"{plot_column}_err"
                    
                    # Check if column exists
                    if plot_column not in df.columns:
                        print(f"Warning: Column '{plot_column}' not found in dataset '{label}'")
                        continue
                    
                    # Get complete residue data (fills missing residues with NaN)
                    residue_data = create_complete_residue_data(df, plot_column, error_column)
                    residues = residue_data['residues']
                    y_data = residue_data['data']
                    y_err = residue_data['errors']
                    
                    # Create label for legend
                    if n_datasets > 1:
                        plot_label = f"{label} {field_suffix}"
                    else:
                        plot_label = field_suffix
                    
                    # Plot grey bars for missing residues first (only once per axis)
                    if field_suffix == 'f1':  # Only plot grey bars once per axis
                        plot_missing_residues_bars(ax, residue_data)
                    
                    # Plot data (only non-NaN values)
                    valid_mask = ~np.isnan(y_data)
                    if np.any(valid_mask):
                        valid_residues = residues[valid_mask]
                        valid_y_data = y_data[valid_mask]
                        valid_y_err = y_err[valid_mask] if y_err is not None else None
                        
                        if valid_y_err is not None:
                            ax.errorbar(valid_residues, valid_y_data, yerr=valid_y_err, 
                                       fmt=field_markers[field_suffix], capsize=3, 
                                       alpha=0.7, color=field_colors[field_suffix], 
                                       label=plot_label)
                        else:
                            ax.plot(valid_residues, valid_y_data, field_markers[field_suffix], 
                                   alpha=0.7, color=field_colors[field_suffix],
                                   label=plot_label)
        else:
            # Original single field plotting
            for j, (df, label) in enumerate(datasets):
                # Get column name with field suffix if needed
                if current_field == 'none':
                    # Field-independent columns
                    plot_column = column
                    error_column = f"{column}_err"
                elif current_field in ['f1', 'f2'] and not column.endswith(f'_{current_field}'):
                    plot_column = f"{column}_{current_field}"
                    error_column = f"{plot_column}_err"
                else:
                    plot_column = column
                    error_column = f"{column}_err"
                
                # Check if column exists
                if plot_column not in df.columns:
                    print(f"Warning: Column '{plot_column}' not found in dataset '{label}'")
                    continue
                
                # Get complete residue data (fills missing residues with NaN)
                residue_data = create_complete_residue_data(df, plot_column, error_column)
                residues = residue_data['residues']
                y_data = residue_data['data']
                y_err = residue_data['errors']
                
                # Plot grey bars for missing residues first (only once per axis)
                if j == 0:  # Only plot grey bars once per axis
                    plot_missing_residues_bars(ax, residue_data)
                
                # Use different markers for multiple datasets
                markers = ['o', 's', '^', 'v', 'D']
                colors = ['blue', 'red', 'green', 'orange', 'purple']
                marker = markers[j % len(markers)]
                color = colors[j % len(colors)]
                
                # Plot data (only non-NaN values)
                valid_mask = ~np.isnan(y_data)
                if np.any(valid_mask):
                    valid_residues = residues[valid_mask]
                    valid_y_data = y_data[valid_mask]
                    valid_y_err = y_err[valid_mask] if y_err is not None else None
                    
                    if valid_y_err is not None:
                        ax.errorbar(valid_residues, valid_y_data, yerr=valid_y_err, fmt=marker, capsize=3, 
                                   alpha=0.7, color=color, label=label if n_datasets > 1 else None)
                    else:
                        ax.plot(valid_residues, valid_y_data, marker, alpha=0.7, color=color,
                               label=label if n_datasets > 1 else None)
        
        # Set labels and title
        ax.set_xlabel('Residue' if i == n_plots-1 else '')  # Only label bottom plot
        
        # Get ylabel - use the column name without field suffix for 'both' mode
        if current_field == 'both':
            ax.set_ylabel(get_ylabel(column))
            ax.set_title(f'{get_plot_title(column)} - Both Fields')
        else:
            # Use the actual plot_column for single field
            ylabel_col = plot_column if 'plot_column' in locals() else column
            ax.set_ylabel(get_ylabel(ylabel_col))
            ax.set_title(get_plot_title(ylabel_col))
        
        ax.grid(True, alpha=0.3)
        
        # Add legend when needed
        if (current_field == 'both') or (n_datasets > 1 and i == 0):
            ax.legend()
        
        # Apply y-limits from configuration
        if len(datasets) > 0:
            # Use the actual plot_column for y-limits
            plot_column_for_limits = plot_column if 'plot_column' in locals() else column
            apply_ylimits(ax, plot_column_for_limits, ylimit_config)
    
    plt.tight_layout()
    
    if output_file:
        # Save with Adobe Illustrator compatible settings
        try:
            plt.savefig(output_file, dpi=300, bbox_inches='tight', fonttype=42)
        except TypeError:
            # Fallback if fonttype parameter is not supported
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"Plot saved to: {output_file}")
    else:
        plt.show()


def get_ylabel(column):
    """Get appropriate y-label for column"""
    ylabel_map = {
        'R1': 'R₁ (s⁻¹)',
        'R2': 'R₂ (s⁻¹)', 
        'hetNOE': 'hetNOE',
        'J0': 'J(0) (ns/rad²)',
        'JwN': 'J(ωₙ) (ns/rad²)',
        'JwH_087': 'J(0.87ωₕ) (ns/rad²)',
        'S2': 'S²',
        'te': 'τₑ (ps)', #'tau_e' 
        'Rex': 'Rₑₓ (s⁻¹)'
    }
    
    for key, label in ylabel_map.items():
        if key in column:
            return label
    return column

def get_plot_title(column):
    """Get appropriate plot title for column"""
    title_map = {
        'R1': 'Longitudinal Relaxation Rate',
        'R2': 'Transverse Relaxation Rate',
        'hetNOE': 'Heteronuclear NOE',
        'J0': 'Spectral Density J(0)',
        'JwN': 'Spectral Density J(ωₙ)',
        'JwH_087': 'Spectral Density J(0.87ωₕ)',
        'S2': 'Order Parameter',
        'te': 'Internal Correlation Time', #tau_e
        'Rex': 'Chemical Exchange' 
    }
    
    for key, title in title_map.items():
        if key in column:
            return title
    return column

def set_ylimits(ax, column, data):
    """Set reasonable y-limits based on data type (legacy function for compatibility)"""
    if 'R1' in column:
        ax.set_ylim(0, max(3, data.max() * 1.1))
    elif 'R2' in column:
        ax.set_ylim(0, max(30, data.max() * 1.1))
    elif 'hetNOE' in column:
        ax.set_ylim(0, 1.2)
    elif 'J0' in column:
        ax.set_ylim(0, data.max() * 1.2)
    elif 'JwN' in column or 'JwH_087' in column:
        ax.set_ylim(0, data.max() * 1.2)
    elif 'S2' in column:
        ax.set_ylim(0, 1)
    elif 'te' in column: #tau_e
        ax.set_ylim(0, data.max() * 1.2)
    elif 'Rex' in column:
        ax.set_ylim(0, data.max() * 1.2)

def interactive_mode():
    """Run script in interactive mode"""
    print("=== Dual-Field NMR Data Plotting Tool (with Y-limits Configuration and Secondary Structure) ===\n")
    
    # Get input files
    datasets = []
    while True:
        file_path = input(f"Enter CSV file path for dataset {len(datasets)+1} (or 'done' to finish): ").strip()
        if file_path.lower() == 'done':
            break
        
        if not os.path.exists(file_path):
            print(f"File not found: {file_path}")
            continue
        
        try:
            df = pd.read_csv(file_path)
            label = input(f"Enter label for this dataset (default: Dataset {len(datasets)+1}): ").strip()
            if not label:
                label = f"Dataset {len(datasets)+1}"
            datasets.append((df, label))
            print(f"Loaded dataset: {label} ({len(df)} rows)\n")
        except Exception as e:
            print(f"Error loading file: {e}")
            continue
    
    if not datasets:
        print("No datasets loaded. Exiting.")
        return
    
    # Get y-limits configuration
    config_file = input("Enter y-limits config file (optional, press Enter to skip): ").strip()
    ylimit_config = load_ylimit_config(config_file) if config_file else None
    
    # Get secondary structure file
    ss_file = input("Enter secondary structure file (optional, press Enter to skip): ").strip()
    ss_map = load_secondary_structure(ss_file) if ss_file else None
    
    # Show available columns
    print("Available columns in first dataset:")
    available = get_available_columns(datasets[0][0])
    
    print("\nRelaxation parameters:")
    for col in available['relaxation']:
        print(f"  {col}")
    
    print("\nSpectral densities:")
    for col in available['spectral_density']:
        print(f"  {col}")
    
    print("\nModel-free parameters:")
    for col in available['model_free']:
        print(f"  {col}")
    
    # Get column selection - allow multiple columns
    print("\nSelect columns to plot (comma-separated):")
    print("Examples: R1,R2,hetNOE,J0,JwN,JwH_087,Rex,S2,te")
    column_input = input("Enter column names: ").strip()
    columns = [col.strip() for col in column_input.split(',')]
    
    # Get field selection for field-dependent columns
    field_independent_cols = ['S2', 'te', 'tc'] #'tau_e'
    rex_cols = ['Rex']  # Rex is field-dependent
    has_field_dependent = any(col not in field_independent_cols 
                             and not any(f'_{f}' in col for f in ['f1', 'f2']) 
                             for col in columns)
    
    if has_field_dependent:
        field = input("Enter field for field-dependent columns (f1/f2/both, default: f1): ").strip() or 'f1'
        if field.lower() == 'both':
            field = 'both'
    else:
        field = 'auto'
    
    # Get output file
    output_file = input("Enter output filename (optional, press Enter for screen display): ").strip()
    if not output_file:
        output_file = None
    elif not output_file.endswith('.pdf'):
        output_file += '.pdf'
    
    # Create multi-column plot
    try:
        plot_multiple_columns(datasets, columns, field, output_file, ylimit_config, ss_map)
    except Exception as e:
        print(f"Error creating plot: {e}")

def main():
    parser = argparse.ArgumentParser(description='Plot dual-field NMR relaxation data with configurable y-axis limits and secondary structure visualization')
    parser.add_argument('files', nargs='*', help='CSV files to plot')
    parser.add_argument('--column', '-c', help='Column(s) to plot (comma-separated for multiple)')
    parser.add_argument('--field', '-f', default='f1', choices=['f1', 'f2', 'both', 'none'], 
                       help='Field to plot (default: f1, use "both" for both fields, "none" for field-independent columns)')
    parser.add_argument('--labels', '-l', nargs='*', help='Labels for datasets')
    parser.add_argument('--output', '-o', help='Output PDF filename')
    parser.add_argument('--ylim-config', help='JSON config file for y-axis limits')
    parser.add_argument('--ss-file', help='Secondary structure file (H/B/L symbols, one per residue)')
    parser.add_argument('--create-example-config', action='store_true',
                       help='Create example y-limits configuration file and exit')
    parser.add_argument('--interactive', '-i', action='store_true', 
                       help='Run in interactive mode')
    
    args = parser.parse_args()
    
    # Handle example config creation
    if args.create_example_config:
        create_example_config()
        return
    
    if args.interactive or not args.files:
        interactive_mode()
        return
    
    # Load y-limits configuration
    ylimit_config = load_ylimit_config(args.ylim_config) if args.ylim_config else None
    
    # Load secondary structure
    ss_map = load_secondary_structure(args.ss_file) if args.ss_file else None
    
    # Load datasets
    datasets = []
    for i, file_path in enumerate(args.files):
        if not os.path.exists(file_path):
            print(f"File not found: {file_path}")
            continue
        
        try:
            df = pd.read_csv(file_path)
            label = args.labels[i] if args.labels and i < len(args.labels) else f"Dataset {i+1}"
            datasets.append((df, label))
        except Exception as e:
            print(f"Error loading {file_path}: {e}")
            continue
    
    if not datasets:
        print("No valid datasets loaded.")
        return
    
    if not args.column:
        print("Available columns:")
        available = get_available_columns(datasets[0][0])
        for category, cols in available.items():
            if cols:
                print(f"\n{category.replace('_', ' ').title()}:")
                for col in cols:
                    print(f"  {col}")
        return
    
    # Parse multiple columns
    columns = [col.strip() for col in args.column.split(',')]
    
    # Auto-detect field handling
    field_independent_cols = ['S2', 'te', 'tc'] #'tau_e',
    has_field_dependent = any(col not in field_independent_cols 
                             and not any(f'_{f}' in col for f in ['f1', 'f2']) 
                             for col in columns)
    
    field = args.field if has_field_dependent else 'auto'
    
    # Create plot
    if len(columns) == 1:
        plot_selected_data(datasets, columns[0], field, args.output, ylimit_config, ss_map)
    else:
        plot_multiple_columns(datasets, columns, field, args.output, ylimit_config, ss_map)

if __name__ == '__main__':
    main()
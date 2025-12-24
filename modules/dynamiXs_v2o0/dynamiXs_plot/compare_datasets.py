#!/usr/bin/env python3
"""
Standalone script for comparing two NMR datasets by subtraction.
Allows comparison of different protein variants (WT, T5D, T6D) at the same field.

User can select which columns to subtract (e.g., R1-R1, hetNOE-hetNOE) and
the results are plotted showing the differences.

MISSING RESIDUE HANDLING (RES VERSION):
- Missing residues are shown as grey bars in both original and difference plots
- Grey bars span from 0 to max value (positive direction) and 0 to min value (negative direction)
- X-axis spans complete residue range from residue 1 to max residue numbers
- Ensures visual continuity and highlights data gaps in comparisons

DUAL FIELD COMPARISON:
- When "both" field is selected, creates TWO separate PDF files
- filename_f1.pdf: Contains comparison for all f1 parameters + field-independent
- filename_f2.pdf: Contains comparison for all f2 parameters + field-independent
- Each PDF shows 2-column layout: original data overlay (left) vs differences (right)

Y-AXIS LIMITS CONFIGURATION (LIMITS VERSION):
- Supports JSON configuration files for custom y-axis limits
- Backward compatible: uses matplotlib auto-scaling if no config provided
- Supports parameter-specific, field-specific, and difference plot limits
- Example config: {"ylimits": {"R1": {"min": 0, "max": 5}}}

SECONDARY STRUCTURE VISUALIZATION (SS VERSION):
- Displays protein secondary structure at the top of plots
- Requires text file with H (helix), B (beta sheet), L (loop) symbols
- One symbol per residue, matching the length of the protein sequence
- Example: HHHHHLLLLBBBBBLLLLHHHH (for 22 residues)
"""

import pandas as pd
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend for thread safety
import matplotlib.pyplot as plt
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

def apply_ylimits(ax, column, plot_type, ylimit_config):
    """
    Apply y-limits from config with fallback to auto-scaling
    
    Parameters:
    -----------
    ax : matplotlib axis
        Plot axis to modify
    column : str
        Column name being plotted
    plot_type : str
        Either "original" or "difference"
    ylimit_config : dict or None
        Configuration dictionary
    """
    if not ylimit_config:
        return  # Use matplotlib auto-scaling
    
    # Determine the config key
    if plot_type == "difference":
        config_key = "difference_ylimits"
        base_column = extract_base_column_name(column)
    else:
        config_key = "ylimits"
        base_column = extract_base_column_name(column)
    
    # Try exact match first, then base column
    limits = None
    if config_key in ylimit_config:
        limits = ylimit_config[config_key].get(column)  # Exact match (e.g., R1_f1)
        if not limits:
            limits = ylimit_config[config_key].get(base_column)  # Base match (e.g., R1)
    
    # Also check field_specific for original plots
    if not limits and plot_type == "original" and "field_specific" in ylimit_config:
        limits = ylimit_config["field_specific"].get(column)
    
    if limits and isinstance(limits, dict) and "min" in limits and "max" in limits:
        try:
            ax.set_ylim(limits["min"], limits["max"])
            print(f"  Applied y-limits for {column} ({plot_type}): [{limits['min']}, {limits['max']}]")
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
        "_comment": "Y-axis limits configuration for NMR data comparison plots",
        "ylimits": {
            "R1": {"min": 0, "max": 5},
            "R2": {"min": 0, "max": 50},
            "hetNOE": {"min": 0, "max": 1.2},
            "J0": {"min": 0, "max": 5}, #e-09
            "JwN": {"min": 0, "max": 2}, #e-10
            "JwH": {"min": 0, "max": 0.5}, #e-11
            "S2": {"min": 0, "max": 1},
            "te": {"min": 0, "max": 100},
            "Rex": {"min": 0, "max": 15},
            "tc": {"min": 0, "max": 50}
        },
        "difference_ylimits": {
            "R1": {"min": -2, "max": 2},
            "R2": {"min": -10, "max": 10},
            "hetNOE": {"min": -0.5, "max": 0.5},
            "J0": {"min": -1, "max": 1},
            "JwN": {"min": -0.5, "max": 0.5},
            "JwH": {"min": -0.2, "max": 0.2},
            "S2": {"min": -0.3, "max": 0.3},
            "te": {"min": -20, "max": 20},
            "Rex": {"min": -5, "max": 5}
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

def load_and_validate_data(file1, file2, label1, label2):
    """Load and validate two datasets for comparison"""
    try:
        df1 = pd.read_csv(file1)
        df2 = pd.read_csv(file2)
    except Exception as e:
        print(f"Error loading files: {e}")
        return None, None
    
    print(f"Dataset 1 ({label1}): {len(df1)} rows")
    print(f"Dataset 2 ({label2}): {len(df2)} rows")
    
    # Check for common residues
    if 'Residue' in df1.columns and 'Residue' in df2.columns:
        common_residues = set(df1['Residue']) & set(df2['Residue'])
        print(f"Common residues: {len(common_residues)}")
        
        if len(common_residues) == 0:
            print("Warning: No common residues found between datasets")
            return df1, df2
        
        # Filter to common residues
        df1_filtered = df1[df1['Residue'].isin(common_residues)].copy()
        df2_filtered = df2[df2['Residue'].isin(common_residues)].copy()
        
        # Sort by residue for proper alignment
        df1_filtered = df1_filtered.sort_values('Residue').reset_index(drop=True)
        df2_filtered = df2_filtered.sort_values('Residue').reset_index(drop=True)
        
        return df1_filtered, df2_filtered
    else:
        print("Warning: No 'Residue' column found, assuming row-by-row alignment")
        return df1, df2

def get_available_columns_for_comparison(df1, df2):
    """Get columns available in both datasets for comparison"""
    common_cols = set(df1.columns) & set(df2.columns)
    
    # Filter out non-data columns
    exclude_cols = {'Residue', 'Residue_Number', 'index'}
    data_cols = [col for col in common_cols if col not in exclude_cols]
    
    # Group by type
    relaxation_cols = [col for col in data_cols if any(x in col for x in ['R1', 'R2', 'hetNOE'])]
    spectral_density_cols = [col for col in data_cols if any(x in col for x in ['J0', 'JwN', 'JwH'])]
    model_free_cols = [col for col in data_cols if any(x in col for x in ['S2', 'te', 'Rex', 'tc'])]
    other_cols = [col for col in data_cols if col not in relaxation_cols + spectral_density_cols + model_free_cols]
    
    return {
        'relaxation': relaxation_cols,
        'spectral_density': spectral_density_cols,
        'model_free': model_free_cols,
        'other': other_cols
    }

def calculate_difference(df1, df2, column, error_propagation=True):
    """
    Calculate difference between two datasets for a given column with complete residue handling
    
    Parameters:
    -----------
    df1, df2 : pd.DataFrame
        Input datasets
    column : str
        Column name to subtract
    error_propagation : bool
        Whether to propagate errors
        
    Returns:
    --------
    dict : Contains complete residue difference data and errors
    """
    if column not in df1.columns or column not in df2.columns:
        raise ValueError(f"Column '{column}' not found in both datasets")
    
    # Get complete residue data for both datasets
    error_col = f"{column}_err"
    residue_data1 = create_complete_residue_data(df1, column, error_col)
    residue_data2 = create_complete_residue_data(df2, column, error_col)
    
    # Find common residue range starting from 1
    min_res = 1  # Always start from residue 1
    max_res = max(residue_data1['residues'].max(), residue_data2['residues'].max())
    complete_residues = np.arange(min_res, max_res + 1)
    
    # Initialize complete arrays
    data1_complete = np.full(len(complete_residues), np.nan)
    data2_complete = np.full(len(complete_residues), np.nan)
    err1_complete = np.full(len(complete_residues), np.nan)
    err2_complete = np.full(len(complete_residues), np.nan)
    
    # Fill in data1
    for i, res in enumerate(residue_data1['residues']):
        idx = int(res) - min_res
        if 0 <= idx < len(complete_residues):
            data1_complete[idx] = residue_data1['data'][i]
            if residue_data1['errors'] is not None:
                err1_complete[idx] = residue_data1['errors'][i]
    
    # Fill in data2
    for i, res in enumerate(residue_data2['residues']):
        idx = int(res) - min_res
        if 0 <= idx < len(complete_residues):
            data2_complete[idx] = residue_data2['data'][i]
            if residue_data2['errors'] is not None:
                err2_complete[idx] = residue_data2['errors'][i]
    
    # Calculate difference (df1 - df2) for valid data points
    diff = data1_complete - data2_complete
    
    # Error propagation if error columns exist
    diff_err = None
    if (error_propagation and 
        residue_data1['errors'] is not None and 
        residue_data2['errors'] is not None):
        # Propagate errors: σ(A-B) = √(σA² + σB²)
        diff_err = np.sqrt(err1_complete**2 + err2_complete**2)
    
    # Calculate max absolute value for missing residue bars
    valid_data = np.concatenate([data1_complete[~np.isnan(data1_complete)], 
                                data2_complete[~np.isnan(data2_complete)]])
    max_abs_value = np.max(np.abs(valid_data)) if len(valid_data) > 0 else 1.0
    
    return {
        'residues': complete_residues,
        'difference': diff,
        'difference_err': diff_err,
        'data1': data1_complete,
        'data2': data2_complete,
        'data1_err': err1_complete,
        'data2_err': err2_complete,
        'max_abs_value': max_abs_value
    }

def plot_multiple_comparisons(diff_data_list, columns, label1, label2, output_file=None, ylimit_config=None, ss_map=None):
    """
    Plot comparison results for multiple parameters in 2-column layout
    Left column: Original data overlay, Right column: Differences
    
    Parameters:
    -----------
    diff_data_list : list of dicts
        List of difference data from calculate_difference for each column
    columns : list of str
        Column names being compared
    label1, label2 : str
        Dataset labels
    output_file : str
        Output filename
    ylimit_config : dict
        Y-axis limits configuration
    ss_map : list
        Secondary structure mapping
    """
    n_params = len(columns)
    
    # Add extra space for secondary structure if provided
    if ss_map:
        # Create layout with SS panel at top
        height_ratios = [1] + [4] * n_params  # SS panel + data panels
        fig, axes = plt.subplots(n_params + 1, 2, figsize=(16, 4*n_params + 1),
                                gridspec_kw={'height_ratios': height_ratios})
        
        # Ensure axes is always 2D array
        if n_params == 0:
            axes = axes.reshape(-1, 2)
        
        # Draw secondary structure in top panels
        if len(diff_data_list) > 0:
            residue_range = (diff_data_list[0]['residues'].min(), diff_data_list[0]['residues'].max())
            
            # Left panel - secondary structure
            draw_secondary_structure(axes[0, 0], residue_range, ss_map)
            axes[0, 0].set_title('Secondary Structure')
            
            # Right panel - secondary structure (same as left)
            draw_secondary_structure(axes[0, 1], residue_range, ss_map)
            axes[0, 1].set_title('Secondary Structure')
        
        # Adjust data panel indices
        data_axes = axes[1:, :]
    else:
        # Original layout without SS
        fig, axes = plt.subplots(n_params, 2, figsize=(16, 4*n_params))
        
        # Ensure axes is always 2D array
        if n_params == 1:
            axes = axes.reshape(1, -1)
        
        data_axes = axes
    
    print(f"Applying y-axis limits from configuration...")
    
    for i, (diff_data, column) in enumerate(zip(diff_data_list, columns)):
        residues = diff_data['residues']
        
        # Left column: Original data overlay
        ax_orig = data_axes[i, 0]
        
        # Plot grey bars for missing residues first
        residue_data_for_bars = {
            'residues': residues,
            'data': diff_data['data1'],  # Use data1 for missing residue reference
            'max_abs_value': diff_data['max_abs_value']
        }
        plot_missing_residues_bars(ax_orig, residue_data_for_bars)
        
        # Plot original data (only non-NaN values)
        valid_mask1 = ~np.isnan(diff_data['data1'])
        valid_mask2 = ~np.isnan(diff_data['data2'])
        
        if np.any(valid_mask1):
            valid_residues1 = residues[valid_mask1]
            valid_data1 = diff_data['data1'][valid_mask1]
            valid_err1 = diff_data['data1_err'][valid_mask1] if diff_data['data1_err'] is not None else None
            
            if valid_err1 is not None and not np.all(np.isnan(valid_err1)):
                ax_orig.errorbar(valid_residues1, valid_data1, yerr=valid_err1, 
                               fmt='o', capsize=3, label=label1, alpha=0.7, color='blue')
            else:
                ax_orig.plot(valid_residues1, valid_data1, 'o', label=label1, alpha=0.7, color='blue')
        
        if np.any(valid_mask2):
            valid_residues2 = residues[valid_mask2]
            valid_data2 = diff_data['data2'][valid_mask2]
            valid_err2 = diff_data['data2_err'][valid_mask2] if diff_data['data2_err'] is not None else None
            
            if valid_err2 is not None and not np.all(np.isnan(valid_err2)):
                ax_orig.errorbar(valid_residues2, valid_data2, yerr=valid_err2, 
                               fmt='s', capsize=3, label=label2, alpha=0.7, color='red')
            else:
                ax_orig.plot(valid_residues2, valid_data2, 's', label=label2, alpha=0.7, color='red')
        
        ax_orig.set_ylabel(get_ylabel(column))
        ax_orig.set_title(f'Original Data: {get_plot_title(column)}')
        if i == 0:  # Legend only on first row
            ax_orig.legend()
        ax_orig.grid(True, alpha=0.3)
        
        # Apply y-limits for original data
        apply_ylimits(ax_orig, column, "original", ylimit_config)
        
        # Only label x-axis on bottom row
        if i == n_params - 1:
            ax_orig.set_xlabel('Residue')
        
        # Right column: Difference
        ax_diff = data_axes[i, 1]
        
        # Plot grey bars for missing residues first
        residue_data_for_diff_bars = {
            'residues': residues,
            'data': diff_data['difference'],
            'max_abs_value': diff_data['max_abs_value']
        }
        plot_missing_residues_bars(ax_diff, residue_data_for_diff_bars)
        
        # Plot difference data (only non-NaN values)
        valid_diff_mask = ~np.isnan(diff_data['difference'])
        if np.any(valid_diff_mask):
            valid_residues_diff = residues[valid_diff_mask]
            valid_difference = diff_data['difference'][valid_diff_mask]
            valid_diff_err = (diff_data['difference_err'][valid_diff_mask] 
                             if diff_data['difference_err'] is not None else None)
            
            if valid_diff_err is not None and not np.all(np.isnan(valid_diff_err)):
                ax_diff.errorbar(valid_residues_diff, valid_difference, yerr=valid_diff_err, 
                               fmt='o', capsize=3, color='red', alpha=0.7)
            else:
                ax_diff.plot(valid_residues_diff, valid_difference, 'o', color='red', alpha=0.7)
        
        # Add zero line
        ax_diff.axhline(y=0, color='black', linestyle='--', alpha=0.5)
        
        ax_diff.set_ylabel(f'Δ{get_ylabel(column)}')
        ax_diff.set_title(f'Difference: {label1} - {label2}')
        ax_diff.grid(True, alpha=0.3)
        
        # Apply y-limits for difference data
        apply_ylimits(ax_diff, column, "difference", ylimit_config)
        
        # Only label x-axis on bottom row
        if i == n_params - 1:
            ax_diff.set_xlabel('Residue')
        
        # Calculate and display statistics
        mean_diff = np.nanmean(diff_data['difference'])
        std_diff = np.nanstd(diff_data['difference'])
        ax_diff.text(0.02, 0.98, f'Mean Δ: {mean_diff:.4f}\nStd Δ: {std_diff:.4f}', 
                   transform=ax_diff.transAxes, verticalalignment='top',
                   bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
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

def plot_comparison(diff_data, column, label1, label2, output_file=None, show_original=True, ylimit_config=None, ss_map=None):
    """
    Plot comparison results for single parameter (backward compatibility)
    
    Parameters:
    -----------
    diff_data : dict
        Difference data from calculate_difference
    column : str
        Column name being compared
    label1, label2 : str
        Dataset labels
    output_file : str
        Output filename
    show_original : bool
        Whether to show original data alongside difference
    ylimit_config : dict
        Y-axis limits configuration
    ss_map : list
        Secondary structure mapping
    """
    # Use the new multi-parameter function with single parameter
    plot_multiple_comparisons([diff_data], [column], label1, label2, output_file, ylimit_config, ss_map)

def get_ylabel(column):
    """Get appropriate y-label for column"""
    ylabel_map = {
        'R1': 'R₁ (s⁻¹)',
        'R2': 'R₂ (s⁻¹)', 
        'hetNOE': 'hetNOE',
        'J0': 'J(0) (ns/rad²)',
        'JwN': 'J(ωₙ) (ns/rad²)',
        'JwH': 'J(ωₕ) (ns/rad²)',
        'S2': 'S²',
        'te': 'τₑ (ps)',
        'Rex': 'Rₑₓ (s⁻¹)',
        'tc': 'τc (ns)'
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
        'JwH': 'Spectral Density J(ωₕ)',
        'S2': 'Order Parameter',
        'te': 'Internal Correlation Time',
        'Rex': 'Chemical Exchange',
        'tc': 'Overall Correlation Time'
    }
    
    for key, title in title_map.items():
        if key in column:
            return title
    return column

def interactive_mode():
    """Run script in interactive mode"""
    print("=== Dataset Comparison Tool (with Y-limits Configuration and Secondary Structure) ===\n")
    
    # Get input files
    file1 = input("Enter path to first dataset (reference): ").strip()
    if not os.path.exists(file1):
        print(f"File not found: {file1}")
        return
    
    file2 = input("Enter path to second dataset (comparison): ").strip()
    if not os.path.exists(file2):
        print(f"File not found: {file2}")
        return
    
    # Get labels
    label1 = input("Enter label for first dataset (default: Dataset 1): ").strip() or "Dataset 1"
    label2 = input("Enter label for second dataset (default: Dataset 2): ").strip() or "Dataset 2"
    
    # Get y-limits configuration
    config_file = input("Enter y-limits config file (optional, press Enter to skip): ").strip()
    ylimit_config = load_ylimit_config(config_file) if config_file else None
    
    # Get secondary structure file
    ss_file = input("Enter secondary structure file (optional, press Enter to skip): ").strip()
    ss_map = load_secondary_structure(ss_file) if ss_file else None
    
    # Load and validate data
    df1, df2 = load_and_validate_data(file1, file2, label1, label2)
    if df1 is None or df2 is None:
        return
    
    # Show available columns
    available = get_available_columns_for_comparison(df1, df2)
    
    print(f"\nAvailable columns for comparison:")
    print("\nRelaxation parameters:")
    for col in available['relaxation']:
        print(f"  {col}")
    
    print("\nSpectral densities:")
    for col in available['spectral_density']:
        print(f"  {col}")
    
    print("\nModel-free parameters:")
    for col in available['model_free']:
        print(f"  {col}")
    
    if available['other']:
        print("\nOther parameters:")
        for col in available['other']:
            print(f"  {col}")
    
    # Get column selection - allow multiple columns
    print("\nSelect columns to compare (comma-separated):")
    print("Examples: R1,R2,hetNOE,J0,JwN,JwH,S2,te,Rex")
    column_input = input("Enter column names: ").strip()
    columns = [col.strip() for col in column_input.split(',')]
    
    # Get field selection for field-dependent columns
    field_independent_cols = ['S2', 'te', 'tc']
    rex_cols = ['Rex']  # Rex has field-dependent versions
    has_field_dependent = any(col not in field_independent_cols 
                             and not any(f'_{f}' in col for f in ['f1', 'f2']) 
                             for col in columns)
    
    if has_field_dependent:
        field = input("Enter field for field-dependent columns (f1/f2/both, default: f1): ").strip() or 'f1'
        if field.lower() == 'both':
            field = 'both'
    else:
        field = 'auto'
    
    # Convert column names to actual column names in dataset
    actual_columns = []
    for column in columns:
        if column in field_independent_cols:
            # Field-independent columns (S2, te, tc)
            if column in df1.columns and column in df2.columns:
                actual_columns.append(column)
            else:
                print(f"Warning: Column '{column}' not found in both datasets - skipping")
        elif field == 'both' and column not in field_independent_cols:
            # Both fields for field-dependent columns
            for field_suffix in ['f1', 'f2']:
                col_name = f"{column}_{field_suffix}"
                if col_name in df1.columns and col_name in df2.columns:
                    actual_columns.append(col_name)
                else:
                    print(f"Warning: Column '{col_name}' not found in both datasets - skipping")
        elif column in rex_cols or not column.endswith(('_f1', '_f2')):
            # Field-dependent columns (R1, R2, hetNOE, J0, JwN, JwH, Rex, etc.)
            if field in ['f1', 'f2']:
                col_name = f"{column}_{field}"
            else:
                col_name = column
                
            if col_name in df1.columns and col_name in df2.columns:
                actual_columns.append(col_name)
            else:
                print(f"Warning: Column '{col_name}' not found in both datasets - skipping")
        else:
            # Column already has field suffix
            if column in df1.columns and column in df2.columns:
                actual_columns.append(column)
            else:
                print(f"Warning: Column '{column}' not found in both datasets - skipping")
    
    valid_columns = actual_columns
    
    if not valid_columns:
        print("No valid columns found in both datasets")
        return
    
    # Get output file
    output_file = input("Enter output filename (optional, press Enter for screen): ").strip()
    if not output_file:
        output_file = None
    elif not output_file.endswith('.pdf'):
        output_file += '.pdf'
    
    # Calculate differences for all columns
    try:
        if field == 'both':
            # Create separate comparisons for f1 and f2
            for field_suffix in ['f1', 'f2']:
                print(f"\nProcessing field {field_suffix}...")
                
                # Filter columns for this field
                field_specific_columns = []
                for column in valid_columns:
                    if column.endswith(f'_{field_suffix}'):
                        field_specific_columns.append(column)
                    elif column in field_independent_cols:
                        field_specific_columns.append(column)
                
                if not field_specific_columns:
                    print(f"No valid columns for field {field_suffix}")
                    continue
                
                # Calculate differences for this field
                diff_data_list = []
                for column in field_specific_columns:
                    diff_data = calculate_difference(df1, df2, column)
                    diff_data_list.append(diff_data)
                
                # Create field-specific output filename
                if output_file:
                    base_name = output_file.replace('.pdf', '')
                    field_output_file = f"{base_name}_{field_suffix}.pdf"
                else:
                    field_output_file = None
                
                # Plot comparisons for this field
                plot_multiple_comparisons(diff_data_list, field_specific_columns, 
                                        label1, label2, field_output_file, ylimit_config, ss_map)
                
                # Print summary statistics for this field
                print(f"\nComparison Summary for {field_suffix} ({label1} - {label2}):")
                for column, diff_data in zip(field_specific_columns, diff_data_list):
                    print(f"  {column}:")
                    print(f"    Mean difference: {np.nanmean(diff_data['difference']):.6f}")
                    print(f"    Std difference: {np.nanstd(diff_data['difference']):.6f}")
                    print(f"    Max difference: {np.nanmax(diff_data['difference']):.6f}")
                    print(f"    Min difference: {np.nanmin(diff_data['difference']):.6f}")
        else:
            # Single field comparison (original behavior)
            diff_data_list = []
            for column in valid_columns:
                diff_data = calculate_difference(df1, df2, column)
                diff_data_list.append(diff_data)
            
            # Plot multiple comparisons
            plot_multiple_comparisons(diff_data_list, valid_columns, label1, label2, output_file, ylimit_config, ss_map)
            
            # Print summary statistics
            print(f"\nComparison Summary ({label1} - {label2}):")
            for i, (column, diff_data) in enumerate(zip(valid_columns, diff_data_list)):
                print(f"  {column}:")
                print(f"    Mean difference: {np.nanmean(diff_data['difference']):.6f}")
                print(f"    Std difference: {np.nanstd(diff_data['difference']):.6f}")
                print(f"    Max difference: {np.nanmax(diff_data['difference']):.6f}")
                print(f"    Min difference: {np.nanmin(diff_data['difference']):.6f}")
        
    except Exception as e:
        print(f"Error in comparison: {e}")

def main():
    parser = argparse.ArgumentParser(description='Compare two NMR datasets by subtraction with configurable y-axis limits and secondary structure visualization')
    parser.add_argument('file1', nargs='?', help='First dataset (reference)')
    parser.add_argument('file2', nargs='?', help='Second dataset (comparison)')
    parser.add_argument('--column', '-c', help='Column(s) to compare (comma-separated for multiple)')
    parser.add_argument('--label1', help='Label for first dataset')
    parser.add_argument('--label2', help='Label for second dataset')
    parser.add_argument('--field', '-f', default='f1', choices=['f1', 'f2', 'both'], 
                       help='Field to compare (default: f1, use "both" for both fields)')
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
    
    if args.interactive or not (args.file1 and args.file2):
        interactive_mode()
        return
    
    # Validate files
    if not os.path.exists(args.file1):
        print(f"File not found: {args.file1}")
        return
    if not os.path.exists(args.file2):
        print(f"File not found: {args.file2}")
        return
    
    # Load y-limits configuration
    ylimit_config = load_ylimit_config(args.ylim_config) if args.ylim_config else None
    
    # Load secondary structure
    ss_map = load_secondary_structure(args.ss_file) if args.ss_file else None
    
    # Set defaults
    label1 = args.label1 or "Dataset 1"
    label2 = args.label2 or "Dataset 2"
    
    # Load data
    df1, df2 = load_and_validate_data(args.file1, args.file2, label1, label2)
    if df1 is None or df2 is None:
        return
    
    if not args.column:
        print("Available columns for comparison:")
        available = get_available_columns_for_comparison(df1, df2)
        for category, cols in available.items():
            if cols:
                print(f"\n{category.replace('_', ' ').title()}:")
                for col in cols:
                    print(f"  {col}")
        return
    
    # Parse multiple columns
    columns = [col.strip() for col in args.column.split(',')]
    
    # Get field handling
    field_independent_cols = ['S2', 'te', 'tc']
    rex_cols = ['Rex']  # Rex has field-dependent versions
    has_field_dependent = any(col not in field_independent_cols 
                             and not any(f'_{f}' in col for f in ['f1', 'f2']) 
                             for col in columns)
    
    field = args.field if has_field_dependent else 'auto'
    
    # Convert column names to actual column names in dataset
    actual_columns = []
    for column in columns:
        if column in field_independent_cols:
            # Field-independent columns (S2, te, tc)
            if column in df1.columns and column in df2.columns:
                actual_columns.append(column)
            else:
                print(f"Warning: Column '{column}' not found in both datasets - skipping")
        elif field == 'both' and column not in field_independent_cols:
            # Both fields for field-dependent columns
            for field_suffix in ['f1', 'f2']:
                col_name = f"{column}_{field_suffix}"
                if col_name in df1.columns and col_name in df2.columns:
                    actual_columns.append(col_name)
                else:
                    print(f"Warning: Column '{col_name}' not found in both datasets - skipping")
        elif column in rex_cols or not column.endswith(('_f1', '_f2')):
            # Field-dependent columns (R1, R2, hetNOE, J0, JwN, JwH, Rex, etc.)
            if field in ['f1', 'f2']:
                col_name = f"{column}_{field}"
            else:
                col_name = column
                
            if col_name in df1.columns and col_name in df2.columns:
                actual_columns.append(col_name)
            else:
                print(f"Warning: Column '{col_name}' not found in both datasets - skipping")
        else:
            # Column already has field suffix
            if column in df1.columns and column in df2.columns:
                actual_columns.append(column)
            else:
                print(f"Warning: Column '{column}' not found in both datasets - skipping")
    
    valid_columns = actual_columns
    
    if not valid_columns:
        print("No valid columns found in both datasets")
        return
    
    # Calculate and plot
    try:
        if field == 'both':
            # Create separate comparisons for f1 and f2
            field_independent_cols = ['S2', 'te', 'tc']  # Define here for command line
            for field_suffix in ['f1', 'f2']:
                print(f"\nProcessing field {field_suffix}...")
                
                # Filter columns for this field
                field_specific_columns = []
                for column in valid_columns:
                    if column.endswith(f'_{field_suffix}'):
                        field_specific_columns.append(column)
                    elif column in field_independent_cols:
                        field_specific_columns.append(column)
                
                if not field_specific_columns:
                    print(f"No valid columns for field {field_suffix}")
                    continue
                
                # Calculate differences for this field
                diff_data_list = []
                for column in field_specific_columns:
                    diff_data = calculate_difference(df1, df2, column)
                    diff_data_list.append(diff_data)
                
                # Create field-specific output filename
                if args.output:
                    base_name = args.output.replace('.pdf', '')
                    field_output_file = f"{base_name}_{field_suffix}.pdf"
                else:
                    field_output_file = None
                
                # Plot comparisons for this field
                plot_multiple_comparisons(diff_data_list, field_specific_columns, 
                                        label1, label2, field_output_file, ylimit_config, ss_map)
                
                # Print summary statistics for this field
                print(f"\nComparison Summary for {field_suffix} ({label1} - {label2}):")
                for column, diff_data in zip(field_specific_columns, diff_data_list):
                    print(f"  {column}:")
                    print(f"    Mean difference: {np.nanmean(diff_data['difference']):.6f}")
                    print(f"    Std difference: {np.nanstd(diff_data['difference']):.6f}")
        else:
            # Single field comparison (original behavior)
            diff_data_list = []
            for column in valid_columns:
                diff_data = calculate_difference(df1, df2, column)
                diff_data_list.append(diff_data)
            
            # Plot multiple comparisons
            plot_multiple_comparisons(diff_data_list, valid_columns, label1, label2, args.output, ylimit_config, ss_map)
            
            # Print summary
            print(f"\nComparison Summary ({label1} - {label2}):")
            for column, diff_data in zip(valid_columns, diff_data_list):
                print(f"  {column}:")
                print(f"    Mean difference: {np.nanmean(diff_data['difference']):.6f}")
                print(f"    Std difference: {np.nanstd(diff_data['difference']):.6f}")
        
    except Exception as e:
        print(f"Error in comparison: {e}")

if __name__ == '__main__':
    main()
    
## python compare_datasetsRESlimitsSSPD.py T5D.csv WT.csv -c R1,R2,hetNOE,J0,JwN,JwH,S2,te,Rex --ylim-config ylimits.json --label1 T5D --label2 WT --field both --output T5DminusWT --ss-file SS.txt
## python compare_datasetsRESlimitsSSPD.py T6D.csv WT.csv -c R1,R2,hetNOE,J0,JwN,JwH,S2,te,Rex --ylim-config ylimits.json --label1 T6D --label2 WT --field both --output T6DminusWT --ss-file SS.txt
## python compare_datasetsRESlimitsSSPD.py T5D.csv T6D.csv -c R1,R2,hetNOE,J0,JwN,JwH,S2,te,Rex --ylim-config ylimits.json --label1 T5D --label2 T6D --field both --output T5DminusT6D --ss-file SS.txt

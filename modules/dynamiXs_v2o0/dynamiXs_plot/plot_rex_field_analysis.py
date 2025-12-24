#!/usr/bin/env python3
"""
Standalone script for plotting Rex Field Dependence and Rex Field Scaling plots
from dual-field NMR relaxation data CSV files.

This script creates two plots side-by-side:
1. Rex Field Dependence: Scatter plot of Rex_f1 vs Rex_f2 with diagonal line
2. Rex Field Scaling: Rex ratio vs expected field ratio with filtering

The Rex Field Scaling plot uses Option 3 filtering:
- Minimum Rex threshold (default 0.5 s⁻¹) at both fields
- Outlier removal (ratio must be within 0.3× to 3× expected)

Usage:
    python plot_rex_field_analysis.py <csv_file> [--field1 600] [--field2 700] [--min-rex 0.5] [--output output.pdf]

Example:
    python plot_rex_field_analysis.py ../density_function_macro/087_WT_density_basic.csv --field1 600 --field2 700 --output rex_analysis_WT.pdf
"""

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import argparse
import os
import sys

# Configure matplotlib for Adobe Illustrator compatibility
matplotlib.rcParams['pdf.fonttype'] = 42  # TrueType fonts
matplotlib.rcParams['ps.fonttype'] = 42   # PostScript compatibility
matplotlib.rcParams['font.family'] = 'sans-serif'
matplotlib.rcParams['font.sans-serif'] = ['DejaVu Sans', 'Liberation Sans', 'Arial']


def detect_field_frequencies(df):
    """
    Try to detect field frequencies from column names or data

    Parameters:
    -----------
    df : pd.DataFrame
        Input dataframe

    Returns:
    --------
    tuple : (field1_freq, field2_freq) or (None, None) if cannot detect
    """
    # Common field combinations
    common_fields = {
        ('600', '700'): (600, 700),
        ('500', '600'): (500, 600),
        ('700', '800'): (700, 800),
        ('800', '900'): (800, 900),
    }

    # Check if there's metadata in the dataframe
    # (This is dataset-specific, adjust as needed)

    return None, None  # Default: cannot auto-detect


def load_and_validate_csv(csv_file, rex_col_f1='Rex_f1', rex_col_f2='Rex_f2'):
    """
    Load CSV file and validate it has required Rex columns

    Parameters:
    -----------
    csv_file : str
        Path to CSV file
    rex_col_f1 : str
        Column name for Rex at field 1
    rex_col_f2 : str
        Column name for Rex at field 2

    Returns:
    --------
    pd.DataFrame : Loaded dataframe
    """
    if not os.path.exists(csv_file):
        raise FileNotFoundError(f"CSV file not found: {csv_file}")

    try:
        df = pd.read_csv(csv_file)
    except Exception as e:
        raise ValueError(f"Error reading CSV file: {e}")

    # Check for required columns
    if rex_col_f1 not in df.columns:
        # Try alternative naming
        if 'Rex_field1' in df.columns:
            rex_col_f1 = 'Rex_field1'
        else:
            raise ValueError(f"Required column '{rex_col_f1}' not found in CSV. Available columns: {list(df.columns)}")

    if rex_col_f2 not in df.columns:
        # Try alternative naming
        if 'Rex_field2' in df.columns:
            rex_col_f2 = 'Rex_field2'
        else:
            raise ValueError(f"Required column '{rex_col_f2}' not found in CSV. Available columns: {list(df.columns)}")

    print(f"Loaded CSV file: {csv_file}")
    print(f"  Total rows: {len(df)}")
    print(f"  Using columns: {rex_col_f1}, {rex_col_f2}")

    return df, rex_col_f1, rex_col_f2


def plot_rex_field_analysis(df, field1_freq, field2_freq,
                            rex_col_f1='Rex_f1', rex_col_f2='Rex_f2',
                            min_rex=0.5, output_file=None, title_prefix=''):
    """
    Create Rex Field Dependence and Rex Field Scaling plots

    Parameters:
    -----------
    df : pd.DataFrame
        Input dataframe with Rex data
    field1_freq : float
        Field 1 frequency in MHz
    field2_freq : float
        Field 2 frequency in MHz
    rex_col_f1 : str
        Column name for Rex at field 1
    rex_col_f2 : str
        Column name for Rex at field 2
    min_rex : float
        Minimum Rex threshold for scaling analysis (s⁻¹)
    output_file : str
        Output PDF filename (optional)
    title_prefix : str
        Prefix for plot titles (e.g., dataset name)
    """
    # Filter out NaN values
    good_data = df[[rex_col_f1, rex_col_f2]].dropna()

    if len(good_data) == 0:
        print("Error: No valid Rex data found in dataframe")
        return

    print(f"\nRex data statistics:")
    print(f"  Valid data points: {len(good_data)}")
    print(f"  Rex_f1 range: [{good_data[rex_col_f1].min():.2f}, {good_data[rex_col_f1].max():.2f}] s⁻¹")
    print(f"  Rex_f2 range: [{good_data[rex_col_f2].min():.2f}, {good_data[rex_col_f2].max():.2f}] s⁻¹")

    # Create figure with 2 subplots
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # =================================================================
    # LEFT PLOT: Rex Field Dependence
    # =================================================================
    ax_dep = axes[0]

    # Scatter plot
    ax_dep.scatter(good_data[rex_col_f1], good_data[rex_col_f2], alpha=0.6, s=50)

    # Diagonal line (y=x)
    max_val = max(good_data[[rex_col_f1, rex_col_f2]].max().max(), 1)
    ax_dep.plot([0, max_val], [0, max_val], 'k--', alpha=0.5, linewidth=1.5, label='y=x')

    # Labels and formatting
    ax_dep.set_xlabel(f'Rex {field1_freq} MHz (s⁻¹)', fontsize=11)
    ax_dep.set_ylabel(f'Rex {field2_freq} MHz (s⁻¹)', fontsize=11)

    if title_prefix:
        ax_dep.set_title(f'{title_prefix}: Rex Field Dependence', fontsize=12, fontweight='bold')
    else:
        ax_dep.set_title('Rex Field Dependence', fontsize=12, fontweight='bold')

    ax_dep.grid(True, alpha=0.3)
    ax_dep.legend()

    # =================================================================
    # RIGHT PLOT: Rex Field Scaling (with Option 3 filtering)
    # =================================================================
    ax_scale = axes[1]

    # Calculate expected field ratio squared
    field_ratio = (field2_freq / field1_freq)**2

    # Option 3: Filter by minimum threshold + outlier removal
    # Filter by minimum threshold at both fields
    mask = (good_data[rex_col_f1] > min_rex) & (good_data[rex_col_f2] > min_rex)
    filtered_rex_data = good_data[mask].copy()

    print(f"\nRex Field Scaling filtering:")
    print(f"  Min Rex threshold: {min_rex} s⁻¹")
    print(f"  Data points after threshold filter: {len(filtered_rex_data)}/{len(good_data)}")

    if len(filtered_rex_data) > 0:
        # Calculate ratios without offset (not needed after filtering)
        rex_ratio = filtered_rex_data[rex_col_f2] / filtered_rex_data[rex_col_f1]

        # Remove outliers (ratio > 3× expected or < 0.3× expected)
        ratio_mask = (rex_ratio < 3 * field_ratio) & (rex_ratio > 0.3 * field_ratio)
        final_filtered_data = filtered_rex_data[ratio_mask]
        final_rex_ratio = rex_ratio[ratio_mask]

        print(f"  Data points after outlier removal: {len(final_filtered_data)}/{len(filtered_rex_data)}")
        print(f"  Final data points plotted: {len(final_filtered_data)}/{len(good_data)}")

        if len(final_filtered_data) > 0:
            # Calculate statistics
            mean_ratio = final_rex_ratio.mean()
            std_ratio = final_rex_ratio.std()

            print(f"  Rex ratio statistics:")
            print(f"    Mean: {mean_ratio:.2f} (expected: {field_ratio:.2f})")
            print(f"    Std: {std_ratio:.2f}")
            print(f"    Range: [{final_rex_ratio.min():.2f}, {final_rex_ratio.max():.2f}]")

        # Plot filtered data
        ax_scale.scatter(np.full(len(final_filtered_data), field_ratio), final_rex_ratio,
                        alpha=0.6, s=50, color='steelblue', edgecolors='black', linewidth=0.5)

        # Expected ratio line
        ax_scale.axhline(y=field_ratio, color='red', linestyle='--', linewidth=2,
                        label=f'Expected ratio = {field_ratio:.2f}')

        # Add info about filtering
        n_total = len(good_data)
        n_filtered = len(final_filtered_data)
        info_text = f'n = {n_filtered}/{n_total}\n(Rex > {min_rex} s⁻¹)'
        ax_scale.text(0.05, 0.95, info_text,
                     transform=ax_scale.transAxes, verticalalignment='top',
                     bbox=dict(boxstyle='round', facecolor='white', alpha=0.8), fontsize=10)

        # Set y-axis limits to show reasonable range
        y_margin = 0.2 * field_ratio
        y_min = max(0, field_ratio - y_margin)
        y_max = field_ratio + y_margin
        if len(final_rex_ratio) > 0:
            data_min = final_rex_ratio.min()
            data_max = final_rex_ratio.max()
            y_min = min(y_min, data_min - 0.1 * (data_max - data_min))
            y_max = max(y_max, data_max + 0.1 * (data_max - data_min))
        ax_scale.set_ylim(y_min, y_max)

    else:
        # No data passed filtering
        ax_scale.text(0.5, 0.5, f'No residues with Rex > {min_rex} s⁻¹\nat both fields',
                     transform=ax_scale.transAxes, ha='center', va='center', fontsize=12,
                     bbox=dict(boxstyle='round', facecolor='yellow', alpha=0.5))

        print(f"  WARNING: No data points passed the minimum Rex threshold!")

    # Labels and formatting
    ax_scale.set_xlabel('Field Ratio² (B₂²/B₁²)', fontsize=11)
    ax_scale.set_ylabel('Rex Ratio (Rex₂/Rex₁)', fontsize=11)

    if title_prefix:
        ax_scale.set_title(f'{title_prefix}: Rex Field Scaling', fontsize=12, fontweight='bold')
    else:
        ax_scale.set_title('Rex Field Scaling', fontsize=12, fontweight='bold')

    ax_scale.grid(True, alpha=0.3)
    ax_scale.legend()

    # Overall figure title
    if title_prefix:
        fig.suptitle(f'{title_prefix}: Rex Field Analysis ({field1_freq}/{field2_freq} MHz)',
                    fontsize=14, fontweight='bold', y=1.02)

    plt.tight_layout()

    # Save or show
    if output_file:
        try:
            plt.savefig(output_file, dpi=300, bbox_inches='tight', fonttype=42)
            print(f"\nPlot saved to: {output_file}")
        except TypeError:
            # Fallback if fonttype parameter is not supported
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            print(f"\nPlot saved to: {output_file}")
    else:
        plt.show()


def main():
    parser = argparse.ArgumentParser(
        description='Plot Rex Field Dependence and Rex Field Scaling from dual-field NMR data',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage with auto-detected columns
  python plot_rex_field_analysis.py data.csv --field1 600 --field2 700

  # Specify custom Rex columns
  python plot_rex_field_analysis.py data.csv --field1 600 --field2 700 --rex-col-f1 Rex_field1 --rex-col-f2 Rex_field2

  # Adjust minimum Rex threshold
  python plot_rex_field_analysis.py data.csv --field1 600 --field2 700 --min-rex 1.0

  # Save to PDF
  python plot_rex_field_analysis.py data.csv --field1 600 --field2 700 --output rex_analysis.pdf
"""
    )

    parser.add_argument('csv_file', help='CSV file with dual-field Rex data')
    parser.add_argument('--field1', type=float, required=True,
                       help='Field 1 frequency in MHz (e.g., 600)')
    parser.add_argument('--field2', type=float, required=True,
                       help='Field 2 frequency in MHz (e.g., 700)')
    parser.add_argument('--rex-col-f1', default='Rex_f1',
                       help='Column name for Rex at field 1 (default: Rex_f1)')
    parser.add_argument('--rex-col-f2', default='Rex_f2',
                       help='Column name for Rex at field 2 (default: Rex_f2)')
    parser.add_argument('--min-rex', type=float, default=0.5,
                       help='Minimum Rex threshold in s⁻¹ for scaling analysis (default: 0.5)')
    parser.add_argument('--output', '-o', help='Output PDF filename (optional)')
    parser.add_argument('--title', help='Title prefix for plots (e.g., dataset name)')

    args = parser.parse_args()

    # Load and validate CSV
    try:
        df, rex_col_f1, rex_col_f2 = load_and_validate_csv(
            args.csv_file, args.rex_col_f1, args.rex_col_f2
        )
    except (FileNotFoundError, ValueError) as e:
        print(f"Error: {e}")
        sys.exit(1)

    # Create plots
    plot_rex_field_analysis(
        df=df,
        field1_freq=args.field1,
        field2_freq=args.field2,
        rex_col_f1=rex_col_f1,
        rex_col_f2=rex_col_f2,
        min_rex=args.min_rex,
        output_file=args.output,
        title_prefix=args.title or ''
    )


if __name__ == '__main__':
    main()

# Rex Field Analysis Plotting Script

## Overview

`plot_rex_field_analysis.py` creates two side-by-side plots for analyzing Rex (chemical exchange) field dependence from dual-field NMR relaxation data:

1. **Rex Field Dependence**: Scatter plot showing Rex at field 2 vs Rex at field 1, with a diagonal reference line (y=x)
2. **Rex Field Scaling**: Tests whether Rex scales quadratically with magnetic field strength (expected for chemical exchange)

## Key Features

### Intelligent Filtering (Option 3)

The Rex Field Scaling plot implements robust filtering to avoid issues with small/noisy Rex values:

1. **Minimum Rex Threshold**: Only includes residues with Rex > threshold (default 0.5 s⁻¹) at **both** fields
2. **Outlier Removal**: Excludes ratios that deviate significantly from expected (outside 0.3× to 3× expected ratio)
3. **Visual Feedback**: Displays filtering statistics on the plot (n filtered/n total)

### Why This Filtering Matters

**Problem**: Small Rex values (< 0.5 s⁻¹) are often within noise:
- Dividing two small noisy numbers produces unreliable ratios
- Example: 0.05 / 0.01 = 5.0 (appears to scale incorrectly, but both values are noise)

**Solution**: Filter out residues without significant chemical exchange before calculating ratios

## Usage

### Basic Usage

```bash
python plot_rex_field_analysis.py <csv_file> --field1 600 --field2 700
```

### Common Examples

```bash
# Analyze WT data at 600/700 MHz
python plot_rex_field_analysis.py ../../density_function_macro/087_WT_density_basic.csv \
    --field1 600 --field2 700 --title "WT" --output rex_analysis_WT.pdf

# Analyze T5D variant with stricter threshold
python plot_rex_field_analysis.py ../../density_function_macro/087_T5D_density_basic.csv \
    --field1 600 --field2 700 --title "T5D" --min-rex 1.0 --output rex_analysis_T5D.pdf

# Use alternative column names (if your CSV uses Rex_field1/Rex_field2)
python plot_rex_field_analysis.py data.csv \
    --field1 600 --field2 700 \
    --rex-col-f1 Rex_field1 --rex-col-f2 Rex_field2 \
    --output rex_analysis.pdf
```

### Command-Line Options

| Option | Description | Default |
|--------|-------------|---------|
| `csv_file` | Input CSV file with Rex data (positional) | Required |
| `--field1` | Field 1 frequency in MHz | Required |
| `--field2` | Field 2 frequency in MHz | Required |
| `--rex-col-f1` | Column name for Rex at field 1 | `Rex_f1` |
| `--rex-col-f2` | Column name for Rex at field 2 | `Rex_f2` |
| `--min-rex` | Minimum Rex threshold (s⁻¹) | 0.5 |
| `--output`, `-o` | Output PDF filename | Display on screen |
| `--title` | Title prefix for plots | None |

## Input File Format

The script expects a CSV file with at least two columns containing Rex values:
- Default: `Rex_f1` and `Rex_f2`
- Alternative: `Rex_field1` and `Rex_field2` (auto-detected)

Example CSV structure:
```
Index,Rex_f1,Rex_f2,S2,te,...
7,1.34,0.66,0.797,12.3,...
9,0.0,0.0,0.793,14.5,...
10,4.73,3.23,0.755,14.6,...
```

## Output

### Terminal Output

The script prints detailed statistics:
```
Loaded CSV file: 087_WT_density_basic.csv
  Total rows: 45
  Using columns: Rex_f1, Rex_f2

Rex data statistics:
  Valid data points: 42
  Rex_f1 range: [0.00, 8.45] s⁻¹
  Rex_f2 range: [0.00, 7.23] s⁻¹

Rex Field Scaling filtering:
  Min Rex threshold: 0.5 s⁻¹
  Data points after threshold filter: 18/42
  Data points after outlier removal: 16/18
  Final data points plotted: 16/42
  Rex ratio statistics:
    Mean: 1.38 (expected: 1.36)
    Std: 0.15
    Range: [1.12, 1.65]
```

### Plot Output

- **Screen display** (default): Interactive matplotlib window
- **PDF file** (with `--output`): High-quality vector graphics (300 DPI, Adobe Illustrator compatible)

## Physical Interpretation

### Expected Behavior

For pure chemical exchange, Rex should scale with the **square of the magnetic field**:
- Rex ∝ ω² where ω is the Larmor frequency
- ω ∝ B (magnetic field strength)
- Therefore: Rex ∝ B²

### Field Scaling Plot

- **X-axis**: Expected field ratio² = (B₂/B₁)²
  - Example: (700/600)² = 1.36
- **Y-axis**: Measured Rex ratio = Rex₂/Rex₁
- **Red dashed line**: Expected ratio if Rex scales perfectly as B²
- **Scatter points**: Actual ratios for residues with significant Rex

### Interpreting Results

| Observation | Interpretation |
|-------------|---------------|
| Points cluster near red line | Rex follows expected B² scaling → consistent with chemical exchange |
| Points systematically above line | Rex₂ larger than expected → possible additional processes at higher field |
| Points systematically below line | Rex₂ smaller than expected → check data quality or model assumptions |
| Large scatter around line | Heterogeneous exchange processes or experimental uncertainty |
| Few/no points plotted | Limited chemical exchange in protein (most Rex < threshold) |

## Adjusting the Min-Rex Threshold

The default threshold (0.5 s⁻¹) works well for most datasets. Adjust if:

- **Increase threshold** (e.g., 1.0 s⁻¹) if:
  - You want to focus only on strong exchange
  - Your data has high noise levels
  - You're seeing too much scatter in the scaling plot

- **Decrease threshold** (e.g., 0.2 s⁻¹) if:
  - You have very clean data with low errors
  - Chemical exchange is weak in your system
  - You're getting too few points with default threshold

## Integration with Other Scripts

This standalone script uses the **same filtering logic** as:
- `dynamiXs/dynamiXs_density_functions/ZZ_2fields_density_087.py`
- `dynamiXs/dynamiXs_density_functions/ZZ_multi_2fields_density_087.py`
- `density_function_macro/ZZ_2fields_density_087.py`
- `density_function_macro/ZZ_multi_2fields_density_087.py`

All these scripts now implement Option 3 filtering for consistent Rex field scaling analysis.

## Troubleshooting

### "Required column 'Rex_f1' not found"

**Solution**: Your CSV uses different column names. Use `--rex-col-f1` and `--rex-col-f2`:
```bash
python plot_rex_field_analysis.py data.csv --field1 600 --field2 700 \
    --rex-col-f1 Rex_field1 --rex-col-f2 Rex_field2
```

### "No residues with Rex > X s⁻¹ at both fields"

**Solution**: Either:
- Lower the threshold: `--min-rex 0.2`
- Your protein may have limited chemical exchange (this is valid data!)

### Empty/blank scaling plot

**Solution**: Check that:
- Your CSV has valid Rex values (not all zeros or NaN)
- Field frequencies are correct (`--field1` and `--field2`)
- Rex values are in s⁻¹ units (not ms⁻¹)

## References

- Rex field scaling theory: Palmer, A.G. (2004) Chem. Rev. 104, 3623-3640
- Model-free analysis: Lipari & Szabo (1982) J. Am. Chem. Soc. 104, 4546-4559

## Version History

- **v1.0** (2025-01-24): Initial release with Option 3 filtering

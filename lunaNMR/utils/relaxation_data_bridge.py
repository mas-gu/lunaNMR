# ABOUTME: Bridge between LunaNMR series integration output and DynamiXs input format
# ABOUTME: Converts tidy CSV format to transposed delay-based matrix for T1/T2/hetNOE analysis

import os
import re
import pandas as pd
import numpy as np
from typing import Dict, List, Optional, Any

from .delay_extractor import DelayExtractor


class RelaxationDataBridge:
    """
    Bridge between LunaNMR series integration output and DynamiXs input format.

    Converts:
    - LunaNMR tidy CSV (spectrum_name, assignment, height/volume, ...)
    - To DynamiXs transposed format (delays as header, residues as rows)

    Also handles hetNOE sat/unsat intensity extraction.
    """

    def __init__(self, peak_list_path: Optional[str] = None):
        """
        Initialize the bridge.

        Args:
            peak_list_path: Optional path to peak list CSV for assignment lookup
        """
        self.delay_extractor = DelayExtractor()
        self.peak_list = None
        self.assignment_to_aa = {}

        if peak_list_path and os.path.exists(peak_list_path):
            self._load_peak_list(peak_list_path)

    def _load_peak_list(self, peak_list_path: str):
        """Load peak list for assignment -> amino acid mapping."""
        try:
            self.peak_list = pd.read_csv(peak_list_path)

            # Build assignment to amino acid mapping if available
            if 'Assignment' in self.peak_list.columns:
                for _, row in self.peak_list.iterrows():
                    assignment = str(row['Assignment'])
                    # Try to extract amino acid if present in a column
                    if 'AA' in self.peak_list.columns:
                        self.assignment_to_aa[assignment] = str(row['AA'])
                    elif 'Amino_Acid' in self.peak_list.columns:
                        self.assignment_to_aa[assignment] = str(row['Amino_Acid'])
        except Exception as e:
            print(f"Warning: Could not load peak list: {e}")
            self.peak_list = None

    def format_residue_id(
        self,
        assignment: str,
        include_aa: bool = False
    ) -> str:
        """
        Format a residue ID for DynamiXs.

        Args:
            assignment: LunaNMR assignment (e.g., "142", "142.ALA")
            include_aa: Whether to include amino acid code

        Returns:
            Formatted residue ID (e.g., "142" or "142.ALA")
        """
        # Extract numeric part
        numbers = re.findall(r'\d+', str(assignment))
        residue_num = numbers[0] if numbers else str(assignment)

        if not include_aa:
            return residue_num

        # Look up amino acid if available
        aa = self.assignment_to_aa.get(str(assignment), None)
        if aa:
            return f"{residue_num}.{aa}"

        # Check if assignment already contains AA (e.g., "142.ALA")
        if '.' in str(assignment):
            parts = str(assignment).split('.')
            if len(parts) == 2 and len(parts[1]) == 3:
                return assignment

        return residue_num

    def extract_delays_from_tidy(
        self,
        tidy_df: pd.DataFrame
    ) -> Dict[str, float]:
        """
        Extract delay values from spectrum names in tidy dataframe.

        Handles two formats:
        1. Filenames with delay patterns (e.g., "T1_50ms.ft") -> extracts 50.0
        2. Pre-extracted column names (e.g., "50", "50_2") -> parses directly

        Args:
            tidy_df: LunaNMR tidy format dataframe

        Returns:
            Dict mapping spectrum_name -> delay_ms
        """
        delays = {}

        if 'spectrum_name' not in tidy_df.columns:
            return delays

        unique_spectra = tidy_df['spectrum_name'].unique()
        for spectrum_name in unique_spectra:
            name_str = str(spectrum_name)

            # First try standard delay pattern extraction (e.g., "_50ms")
            delay = self.delay_extractor.extract_delay_ms(name_str)

            if delay is None:
                # Try parsing as pre-extracted column name (e.g., "50" or "50_2")
                delay = self._parse_column_name_delay(name_str)

            if delay is not None:
                delays[spectrum_name] = delay

        return delays

    def _parse_column_name_delay(self, name: str) -> Optional[float]:
        """
        Parse delay from a pre-extracted column name.

        Handles formats like:
        - "50" -> 50.0
        - "50_2" -> 50.0 (duplicate indicator ignored for delay value)
        - "100.5" -> 100.5
        - "100.5_3" -> 100.5

        Args:
            name: Column name string

        Returns:
            Delay value in ms, or None if not parseable
        """
        import re

        # Match number (possibly with decimal) at start, optionally followed by _N
        match = re.match(r'^(\d+(?:\.\d+)?)(?:_\d+)?$', name)
        if match:
            try:
                return float(match.group(1))
            except ValueError:
                return None
        return None

    def convert_tidy_to_dynamixs_format(
        self,
        tidy_csv_path: str,
        output_csv_path: str,
        intensity_column: str = 'height'
    ) -> Dict[str, Any]:
        """
        Convert LunaNMR tidy CSV to DynamiXs transposed format.

        LunaNMR tidy format:
            spectrum_name,assignment,peak_number,ppm_x,ppm_y,height,volume,...
            T1_50ms,142,1,7.234,120.345,1000000,9000000,...

        DynamiXs format:
            ,50,100,200,500,...
            142,1000000,900000,730000,500000,...
            143,900000,810000,650000,420000,...

        Args:
            tidy_csv_path: Path to LunaNMR tidy CSV
            output_csv_path: Path to output DynamiXs format CSV
            intensity_column: Column to use for intensity ('height' or 'volume')

        Returns:
            Dict with conversion statistics
        """
        # Load tidy data
        tidy_df = pd.read_csv(tidy_csv_path)

        if tidy_df.empty:
            raise ValueError("No data found in tidy CSV")

        # Validate required columns
        required_cols = ['spectrum_name', 'assignment', intensity_column]
        missing = [c for c in required_cols if c not in tidy_df.columns]
        if missing:
            raise ValueError(f"Missing required columns: {missing}")

        # Extract delays from spectrum names
        delays = self.extract_delays_from_tidy(tidy_df)

        if not delays:
            raise ValueError("Could not extract delays from spectrum names")

        # Get unique assignments (residues)
        unique_assignments = tidy_df['assignment'].unique()

        # Sort delays
        sorted_spectra = sorted(delays.keys(), key=lambda x: delays[x])
        sorted_delays = [delays[s] for s in sorted_spectra]

        # Build output matrix
        # First row: empty + delays
        # Subsequent rows: residue + intensities

        rows = []

        for assignment in sorted(unique_assignments, key=lambda x: self._extract_number(x)):
            row = [self.format_residue_id(assignment)]

            for spectrum_name in sorted_spectra:
                # Find intensity for this assignment in this spectrum
                mask = (
                    (tidy_df['assignment'] == assignment) &
                    (tidy_df['spectrum_name'] == spectrum_name)
                )
                matching = tidy_df[mask]

                if not matching.empty:
                    intensity = matching[intensity_column].iloc[0]
                    row.append(intensity)
                else:
                    row.append(np.nan)

            rows.append(row)

        # Create dataframe header
        # If spectrum_name is already a column name (e.g., "50", "50_2"), use directly
        # Otherwise, convert delay value to column name
        header_items = []
        for spectrum_name in sorted_spectra:
            name_str = str(spectrum_name)
            # Check if already a pre-extracted column name
            if self._parse_column_name_delay(name_str) is not None:
                header_items.append(name_str)
            else:
                # Convert delay value to column name
                delay = delays[spectrum_name]
                if delay == int(delay):
                    header_items.append(str(int(delay)))
                else:
                    header_items.append(f"{delay:.1f}")

        header = [''] + header_items

        # Write output
        with open(output_csv_path, 'w') as f:
            # Write header
            f.write(','.join(header) + '\n')
            # Write data rows
            for row in rows:
                formatted = [str(v) if pd.notna(v) else '' for v in row]
                f.write(','.join(formatted) + '\n')

        return {
            'residue_count': len(rows),
            'delay_count': len(sorted_delays),
            'delays': sorted_delays,
            'residues': [r[0] for r in rows],
            'output_file': output_csv_path
        }

    def convert_hetnoe_to_dynamixs_format(
        self,
        sat_tidy_csv: str,
        unsat_tidy_csv: str,
        output_sat_csv: str,
        output_unsat_csv: str,
        intensity_column: str = 'height'
    ) -> Dict[str, Any]:
        """
        Convert hetNOE saturated/unsaturated tidy CSVs to DynamiXs format.

        For hetNOE, the output is simpler:
            residue,intensity
            142,800000
            143,750000

        Args:
            sat_tidy_csv: Path to saturated spectrum tidy CSV
            unsat_tidy_csv: Path to unsaturated spectrum tidy CSV
            output_sat_csv: Output path for saturated intensities
            output_unsat_csv: Output path for unsaturated intensities
            intensity_column: Column to use ('height' or 'volume')

        Returns:
            Dict with conversion statistics
        """
        sat_df = pd.read_csv(sat_tidy_csv)
        unsat_df = pd.read_csv(unsat_tidy_csv)

        # Extract intensities per residue
        sat_intensities = {}
        unsat_intensities = {}

        for _, row in sat_df.iterrows():
            assignment = str(row['assignment'])
            residue = self.format_residue_id(assignment)
            sat_intensities[residue] = row[intensity_column]

        for _, row in unsat_df.iterrows():
            assignment = str(row['assignment'])
            residue = self.format_residue_id(assignment)
            unsat_intensities[residue] = row[intensity_column]

        # Find common residues
        common = set(sat_intensities.keys()) & set(unsat_intensities.keys())
        common_sorted = sorted(common, key=self._extract_number)

        # Write saturated
        with open(output_sat_csv, 'w') as f:
            for residue in common_sorted:
                f.write(f"{residue},{sat_intensities[residue]}\n")

        # Write unsaturated
        with open(output_unsat_csv, 'w') as f:
            for residue in common_sorted:
                f.write(f"{residue},{unsat_intensities[residue]}\n")

        return {
            'residue_count': len(common),
            'residues': list(common_sorted),
            'output_sat': output_sat_csv,
            'output_unsat': output_unsat_csv
        }

    def _extract_number(self, value: Any) -> int:
        """Extract numeric value for sorting."""
        try:
            numbers = re.findall(r'\d+', str(value))
            return int(numbers[0]) if numbers else 0
        except (ValueError, IndexError):
            return 0

    def build_series_params(
        self,
        use_voigt: bool = True,
        parallel: bool = True,
        fix_positions: bool = False,
        fix_linewidths: bool = False
    ) -> Dict[str, Any]:
        """
        Build parameter dictionary for MultiSpectrumProcessor.

        This helper creates an NMRParameterManager-compatible dict
        with sensible defaults for relaxation series processing.

        Args:
            use_voigt: Whether to use Voigt fitting (vs simple integration)
            parallel: Whether to use parallel processing
            fix_positions: Whether to fix peak positions during fitting
            fix_linewidths: Whether to fix linewidths during fitting

        Returns:
            Dict compatible with MultiSpectrumProcessor initialization
        """
        return {
            'use_voigt_fitting': use_voigt,
            'use_parallel': parallel,
            'fix_positions': fix_positions,
            'fix_linewidths': fix_linewidths,
            'peak_source_mode': 'reference',
            'snr_threshold': 3.0,
            'baseline_correction': True,
            'output_format': 'tidy',
            # Default fitting parameters
            'fitting_window_x': 0.04,
            'fitting_window_y': 0.4,
            'max_iterations': 100,
            'tolerance': 1e-6,
        }

    def process_t1_t2_folder(
        self,
        spectra_folder: str,
        peak_list: str,
        output_csv: str,
        experiment_type: str = 'T1',
        intensity_column: str = 'height',
        n_bootstrap: int = 100
    ) -> str:
        """
        Process a folder of NMR spectra to extract T1/T2/T1rho relaxation rates.

        This method:
        1. Extracts delays from spectrum filenames
        2. Runs LunaNMR series integration (Fit Series) to get intensities
        3. Fits exponential decay to get R1/R2/R1rho values
        4. Outputs a CSV with relaxation rates and errors

        Args:
            spectra_folder: Folder containing spectra with delays in filenames
            peak_list: Path to peak list CSV
            output_csv: Path for output CSV with fitted rates
            experiment_type: 'T1', 'T2', or 'T1rho'
            intensity_column: Column for intensity ('height' or 'volume')
            n_bootstrap: Number of bootstrap iterations for error estimation

        Returns:
            Path to output CSV
        """
        from scipy.optimize import curve_fit

        # Step 1: Scan folder for spectra with delays
        files_with_delays = self.delay_extractor.scan_folder(spectra_folder)
        if not files_with_delays:
            raise ValueError(f"No spectra with delays found in {spectra_folder}")

        # Sort by delay
        files_with_delays.sort(key=lambda x: x[1])
        spectrum_files = [f for f, _ in files_with_delays]
        delays_ms = np.array([d for _, d in files_with_delays])

        # Step 2: Load and prepare peak list
        peak_df = pd.read_csv(peak_list)

        # Normalize peak list columns to what MultiSpectrumProcessor expects
        col_map = {
            'Assignment': 'assignment',
            'Residue': 'assignment',
            'res': 'assignment',
            '1H': 'ppm_x',
            'H_ppm': 'ppm_x',
            'w2': 'ppm_x',
            '15N': 'ppm_y',
            'N_ppm': 'ppm_y',
            'w1': 'ppm_y',
        }
        peak_df = peak_df.rename(columns={k: v for k, v in col_map.items() if k in peak_df.columns})

        # Step 3: Run LunaNMR series integration (Fit Series in cascade mode)
        from lunaNMR.processors.multi_spectrum_processor import MultiSpectrumProcessor

        # Create output folder for series results
        output_folder = os.path.dirname(output_csv)
        if not output_folder:
            output_folder = os.getcwd()

        # Build default parameters for MultiSpectrumProcessor
        voigt_params = self.build_series_params(
            use_voigt=True,
            parallel=True,
            fix_positions=False,
            fix_linewidths=False
        )

        # Process spectra using Fit Series
        processor = MultiSpectrumProcessor(voigt_params)
        batch_results = processor.process_nmr_series(
            nmr_files=spectrum_files,
            reference_peaks=peak_df,
            output_folder=output_folder,
            peak_source_mode='cascade'  # Use cascade for best tracking
        )

        # Read the tidy results file created by process_nmr_series
        tidy_file = os.path.join(output_folder, "series_analysis_tidy.csv")
        if not os.path.exists(tidy_file):
            raise ValueError(f"Series integration did not produce results file: {tidy_file}")

        tidy_df = pd.read_csv(tidy_file)

        # Validate we have data
        if tidy_df.empty:
            raise ValueError("No peak intensities could be extracted from spectra")

        if 'assignment' not in tidy_df.columns:
            raise ValueError(f"Missing 'assignment' column in results. Available columns: {list(tidy_df.columns)}")

        # Step 3: Fit exponential decay for each residue
        # Group by assignment
        residue_results = []

        def exp_decay(t, I0, R):
            """Exponential decay: I(t) = I0 * exp(-R * t)"""
            return I0 * np.exp(-R * t / 1000.0)  # t in ms, R in s^-1

        unique_assignments = tidy_df['assignment'].unique()

        for assignment in unique_assignments:
            mask = tidy_df['assignment'] == assignment
            residue_data = tidy_df[mask].copy()

            # Get intensities sorted by delay
            intensities = []
            for f, d in files_with_delays:
                fname = os.path.basename(f)
                match = residue_data[residue_data['spectrum_name'].str.contains(fname, na=False)]
                if not match.empty:
                    intensities.append(match[intensity_column].iloc[0])
                else:
                    intensities.append(np.nan)

            intensities = np.array(intensities)

            # Skip if too many NaN values
            valid_mask = ~np.isnan(intensities)
            if np.sum(valid_mask) < 3:
                continue

            valid_delays = delays_ms[valid_mask]
            valid_intensities = intensities[valid_mask]

            try:
                # Initial guess
                I0_guess = np.max(valid_intensities)
                R_guess = 1.0 / (valid_delays[len(valid_delays)//2] / 1000.0)  # Rough T estimate

                # Fit exponential
                popt, pcov = curve_fit(
                    exp_decay, valid_delays, valid_intensities,
                    p0=[I0_guess, R_guess],
                    bounds=([0, 0.01], [np.inf, 100]),
                    maxfev=1000
                )

                I0_fit, R_fit = popt
                R_err = np.sqrt(np.diag(pcov))[1] if pcov is not None else 0.0

                # Bootstrap error estimation if requested
                if n_bootstrap > 0:
                    R_bootstrap = []
                    for _ in range(n_bootstrap):
                        idx = np.random.choice(len(valid_delays), size=len(valid_delays), replace=True)
                        try:
                            popt_bs, _ = curve_fit(
                                exp_decay, valid_delays[idx], valid_intensities[idx],
                                p0=[I0_fit, R_fit],
                                bounds=([0, 0.01], [np.inf, 100]),
                                maxfev=500
                            )
                            R_bootstrap.append(popt_bs[1])
                        except Exception:
                            pass

                    if len(R_bootstrap) > 10:
                        R_err = np.std(R_bootstrap)

                # Format residue ID
                residue_id = self.format_residue_id(assignment, include_aa=True)

                # Store result based on experiment type
                if experiment_type.lower() in ['t1', 't1rho']:
                    rate_col = 'R1' if experiment_type.lower() == 't1' else 'R1rho'
                    residue_results.append({
                        'residue': residue_id,
                        rate_col: R_fit,
                        f'{rate_col}_err': R_err
                    })
                else:  # T2
                    residue_results.append({
                        'residue': residue_id,
                        'R2': R_fit,
                        'R2_err': R_err
                    })

            except Exception as e:
                print(f"Warning: Could not fit {assignment}: {e}")
                continue

        # Step 4: Save results
        if not residue_results:
            raise ValueError("No residues could be fitted")

        result_df = pd.DataFrame(residue_results)
        result_df.to_csv(output_csv, index=False)

        return output_csv

    def _simple_extract_intensities(
        self,
        spectrum_files: List[str],
        peak_df: pd.DataFrame,
        delays_ms: np.ndarray
    ) -> pd.DataFrame:
        """
        Simple fallback for intensity extraction without full series processor.

        This is a simplified approach that reads spectra and extracts peak heights
        at specified positions.
        """
        import nmrglue as ng

        rows = []

        for i, spec_file in enumerate(spectrum_files):
            fname = os.path.basename(spec_file)
            delay = delays_ms[i]

            try:
                # Load spectrum
                dic, data = ng.pipe.read(spec_file)

                # Get spectral parameters
                udic = ng.pipe.guess_udic(dic, data)
                uc0 = ng.pipe.make_uc(dic, data, dim=0)  # 15N
                uc1 = ng.pipe.make_uc(dic, data, dim=1)  # 1H

                # Extract intensities at peak positions
                for _, peak in peak_df.iterrows():
                    assignment = peak.get('assignment', peak.name)
                    ppm_x = peak.get('ppm_x', peak.get('H_ppm', peak.get('1H')))
                    ppm_y = peak.get('ppm_y', peak.get('N_ppm', peak.get('15N')))

                    if pd.isna(ppm_x) or pd.isna(ppm_y):
                        continue

                    # Convert ppm to indices
                    idx_y = uc0.i(ppm_y)
                    idx_x = uc1.i(ppm_x)

                    # Get height at peak position (with small window for max)
                    window = 2
                    y_start = max(0, idx_y - window)
                    y_end = min(data.shape[0], idx_y + window + 1)
                    x_start = max(0, idx_x - window)
                    x_end = min(data.shape[1], idx_x + window + 1)

                    region = data[y_start:y_end, x_start:x_end]
                    height = np.max(region) if region.size > 0 else np.nan

                    rows.append({
                        'spectrum_name': fname,
                        'assignment': assignment,
                        'height': height,
                        'delay_ms': delay
                    })

            except Exception as e:
                print(f"Warning: Could not process {spec_file}: {e}")
                continue

        return pd.DataFrame(rows)

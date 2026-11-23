#!/usr/bin/env python3
"""
File Management Module

This module handles file operations, data loading, and file system interactions
for the NMR Peak Series application.

Classes:
- NMRFileManager: Core file operations for NMR data
- DataValidator: Validation of NMR data files
- FileMetadata: File metadata collection and analysis
"""

import os
import pandas as pd

class NMRFileManager:
    """Core file operations for NMR data and peak lists"""

    # Bruker processed data files (no extension)
    BRUKER_PDATA_NAMES = ('2rr', '2ri', '2ir', '2ii', '1r', '1i',
                          '3rrr', '3rri', '3rir', '3rii', '3irr', '3iri', '3iir', '3iii')

    def __init__(self):
        self.supported_nmr_formats = ['ft', 'ft1', 'ser', 'ft2', 'ft3', 'pipe', 'ucsf']
        self.supported_peak_formats = ['txt', 'csv', 'peaks']

    def _is_bruker_pdata(self, file_path):
        """Check if file is Bruker processed data based on filename."""
        basename = os.path.basename(file_path)
        return basename in self.BRUKER_PDATA_NAMES

    def _is_varian_fid(self, file_path):
        """Check if file is Varian/Agilent FID data.

        Varian format uses a directory structure:
        - Directory named *.fid/
        - Contains file 'fid' (no extension) with the data
        - Contains 'procpar' file with parameters
        """
        basename = os.path.basename(file_path)
        parent_dir = os.path.dirname(file_path)

        # Check if this is a file named 'fid' with no extension
        if basename != 'fid':
            return False

        # Check if parent directory name ends with .fid
        parent_name = os.path.basename(parent_dir)
        if not parent_name.endswith('.fid'):
            return False

        # Check if procpar exists (required companion file)
        procpar_path = os.path.join(parent_dir, 'procpar')
        if not os.path.exists(procpar_path):
            return False

        return True

    def validate_nmr_file(self, file_path):
        """Validate NMR data file"""
        if not os.path.exists(file_path):
            return False, "File does not exist"

        try:
            # Check for Varian/Agilent FID data (file named 'fid' in *.fid/ directory)
            if self._is_varian_fid(file_path):
                size_mb = os.path.getsize(file_path) / (1024 * 1024)
                if size_mb < 0.001:
                    return False, "File too small"
                if size_mb > 1000:
                    return False, "File too large"
                return True, "Valid Varian/Agilent FID file"

            # Check for Bruker pdata files (no extension)
            if self._is_bruker_pdata(file_path):
                size_mb = os.path.getsize(file_path) / (1024 * 1024)
                if size_mb < 0.001:
                    return False, "File too small"
                if size_mb > 1000:
                    return False, "File too large"
                return True, "Valid Bruker pdata file"

            ext = os.path.splitext(file_path)[1].lower().lstrip('.')
            if ext not in self.supported_nmr_formats:
                return False, f"Unsupported format: {ext}"

            size_mb = os.path.getsize(file_path) / (1024 * 1024)
            if size_mb < 0.001:  # Less than 1KB
                return False, "File too small"

            if size_mb > 1000:  # More than 1GB
                return False, "File too large"

            return True, "Valid NMR file"

        except Exception as e:
            return False, f"Validation error: {str(e)}"

    def _read_peak_table(self, file_path):
        """A peak list file as a frame with stripped, standardised column names.

        A .csv is comma-separated; any other extension is tried comma-separated first
        and, when that gives a single column, tab-separated.
        """
        ext = os.path.splitext(file_path)[1].lower().lstrip('.')
        if ext == 'csv':
            df = pd.read_csv(file_path)
        else:
            try:
                df = pd.read_csv(file_path, sep=',')
                if len(df.columns) == 1:
                    df = pd.read_csv(file_path, sep='\t')
            except Exception:
                try:
                    df = pd.read_csv(file_path, sep='\t')
                except Exception:
                    raise ValueError("Could not parse file with comma or tab separators")

        df.columns = df.columns.str.strip()
        return self.standardize_peak_columns(df)

    def validate_peak_file(self, file_path):
        """Validate peak list file"""
        if not os.path.exists(file_path):
            return False, "File does not exist"

        try:
            ext = os.path.splitext(file_path)[1].lower().lstrip('.')
            if ext not in self.supported_peak_formats:
                return False, f"Unsupported format: {ext}"

            df = self._read_peak_table(file_path)

            required_cols = ['Position_X', 'Position_Y']
            missing_cols = [col for col in required_cols if col not in df.columns]

            if missing_cols:
                return False, f"Missing columns: {missing_cols}"

            if len(df) == 0:
                return False, "Empty peak list"

            x_valid = pd.to_numeric(df['Position_X'], errors='coerce').notna().all()
            y_valid = pd.to_numeric(df['Position_Y'], errors='coerce').notna().all()

            if not x_valid or not y_valid:
                return False, "Invalid chemical shift values"

            return True, f"Valid peak list ({len(df)} peaks)"

        except Exception as e:
            return False, f"Validation error: {str(e)}"

    def load_peak_list(self, file_path):
        """Load peak list file"""
        try:
            valid, message = self.validate_peak_file(file_path)
            if not valid:
                raise ValueError(message)

            df = self._read_peak_table(file_path)

            # Strip whitespace from assignments: a list written "3LysH, 8.2, 126.3"
            # otherwise carries the space into every output CSV, and residue
            # matching across datasets is by exact string.
            if 'Assignment' in df.columns:
                df['Assignment'] = df['Assignment'].astype(str).str.strip()

            return df

        except Exception as e:
            raise RuntimeError(f"Failed to load peak list: {str(e)}")

    def standardize_peak_columns(self, df):
        """Standardize peak list column names"""
        # Common column name mappings
        column_mappings = {
            'x': 'Position_X',
            'y': 'Position_Y',
            'ppm_x': 'Position_X',
            'ppm_y': 'Position_Y',
            '1H': 'Position_X',
            '15N': 'Position_Y',
            '13C': 'Position_Y',
            'f2': 'Position_X',
            'f1': 'Position_Y',
            'assignment': 'Assignment',
            'assign': 'Assignment',
            'label': 'Assignment',
            'height': 'Height',
            'intensity': 'Intensity',
            'amp': 'Intensity',
            'amplitude': 'Intensity',
            'int': 'Intensity',
            'volume': 'Volume'
        }

        for old_name, new_name in column_mappings.items():
            for col in df.columns:
                if col.lower() == old_name.lower():
                    df.rename(columns={col: new_name}, inplace=True)
                    break

        if 'Position_X' not in df.columns:
            if 'X' in df.columns:
                df.rename(columns={'X': 'Position_X'}, inplace=True)

        if 'Position_Y' not in df.columns:
            if 'Y' in df.columns:
                df.rename(columns={'Y': 'Position_Y'}, inplace=True)

        if 'Assignment' not in df.columns:
            df['Assignment'] = [f'Peak_{i+1}' for i in range(len(df))]

        return df

#!/usr/bin/env python3
"""
File Management Module

This module handles file operations, data loading, and file system interactions
for the NMR Peak Series application.

Classes:
- NMRFileManager: Core file operations for NMR data
- DataValidator: Validation of NMR data files
- FileMetadata: File metadata collection and analysis

Author: Guillaume Mas
Date: 2025
"""

import os
import glob
import pandas as pd
import numpy as np
from pathlib import Path
from datetime import datetime
import json
import warnings
warnings.filterwarnings('ignore')

try:
    import nmrglue as ng
except ImportError:
    print("Warning: nmrglue not available - some features may be limited")
    ng = None

class NMRFileManager:
    """Core file operations for NMR data and peak lists"""

    def __init__(self):
        self.supported_nmr_formats = ['ft', 'fid', 'ser']
        self.supported_peak_formats = ['txt', 'csv', 'peaks']
        self.recent_files = []
        self.max_recent = 10

    def validate_nmr_file(self, file_path):
        """Validate NMR data file"""
        if not os.path.exists(file_path):
            return False, "File does not exist"

        try:
            # Check file extension
            ext = os.path.splitext(file_path)[1].lower().lstrip('.')
            if ext not in self.supported_nmr_formats:
                return False, f"Unsupported format: {ext}"

            # Check file size (basic validation)
            size_mb = os.path.getsize(file_path) / (1024 * 1024)
            if size_mb < 0.001:  # Less than 1KB
                return False, "File too small"

            if size_mb > 1000:  # More than 1GB
                return False, "File too large"

            return True, "Valid NMR file"

        except Exception as e:
            return False, f"Validation error: {str(e)}"

    def validate_peak_file(self, file_path):
        """Validate peak list file"""
        if not os.path.exists(file_path):
            return False, "File does not exist"

        try:
            # Check file extension
            ext = os.path.splitext(file_path)[1].lower().lstrip('.')
            if ext not in self.supported_peak_formats:
                return False, f"Unsupported format: {ext}"

            # Try to read the file with different separators
            df = None
            if ext == 'csv':
                df = pd.read_csv(file_path)
            else:
                # Try comma-separated first (more common), then tab-separated
                try:
                    df = pd.read_csv(file_path, sep=',')
                    # Check if we got multiple columns (successful comma separation)
                    if len(df.columns) == 1:
                        # Probably not comma-separated, try tab
                        df = pd.read_csv(file_path, sep='\t')
                except:
                    try:
                        df = pd.read_csv(file_path, sep='\t')
                    except:
                        return False, "Could not parse file with comma or tab separators"

            # Clean column names (remove extra spaces)
            df.columns = df.columns.str.strip()

            # Check for required columns
            required_cols = ['Position_X', 'Position_Y']
            missing_cols = [col for col in required_cols if col not in df.columns]

            if missing_cols:
                return False, f"Missing columns: {missing_cols}"

            # Check data validity
            if len(df) == 0:
                return False, "Empty peak list"

            # Check for valid chemical shift values
            x_valid = pd.to_numeric(df['Position_X'], errors='coerce').notna().all()
            y_valid = pd.to_numeric(df['Position_Y'], errors='coerce').notna().all()

            if not x_valid or not y_valid:
                return False, "Invalid chemical shift values"

            return True, f"Valid peak list ({len(df)} peaks)"

        except Exception as e:
            return False, f"Validation error: {str(e)}"

    def load_nmr_data(self, file_path):
        """Load NMR data using nmrglue"""
        if not ng:
            raise ImportError("nmrglue is required for NMR data loading")

        try:
            # Validate file first
            valid, message = self.validate_nmr_file(file_path)
            if not valid:
                raise ValueError(message)

            # Load based on file type
            ext = os.path.splitext(file_path)[1].lower().lstrip('.')

            if ext in ['ft', 'fid']:
                dic, data = ng.varian.read(os.path.dirname(file_path))
            else:
                raise ValueError(f"Unsupported format for loading: {ext}")

            # Add to recent files
            self.add_recent_file(file_path)

            return dic, data

        except Exception as e:
            raise RuntimeError(f"Failed to load NMR data: {str(e)}")

    def load_peak_list(self, file_path):
        """Load peak list file"""
        try:
            # Validate file first
            valid, message = self.validate_peak_file(file_path)
            if not valid:
                raise ValueError(message)

            # Determine file type and load
            ext = os.path.splitext(file_path)[1].lower().lstrip('.')

            if ext == 'csv':
                df = pd.read_csv(file_path)
            else:
                # Try comma-separated first (more common), then tab-separated
                try:
                    df = pd.read_csv(file_path, sep=',')
                    # Check if we got multiple columns (successful comma separation)
                    if len(df.columns) == 1:
                        # Probably not comma-separated, try tab
                        df = pd.read_csv(file_path, sep='\t')
                except:
                    try:
                        df = pd.read_csv(file_path, sep='\t')
                    except:
                        raise ValueError("Could not parse file with comma or tab separators")

            # Clean column names (remove extra spaces)
            df.columns = df.columns.str.strip()

            # Standardize column names
            df = self.standardize_peak_columns(df)

            # Add to recent files
            self.add_recent_file(file_path)

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
            'label': 'Assignment'
        }

        # Apply mappings (case insensitive)
        for old_name, new_name in column_mappings.items():
            for col in df.columns:
                if col.lower() == old_name.lower():
                    df.rename(columns={col: new_name}, inplace=True)
                    break

        # Ensure required columns exist
        if 'Position_X' not in df.columns:
            if 'X' in df.columns:
                df.rename(columns={'X': 'Position_X'}, inplace=True)

        if 'Position_Y' not in df.columns:
            if 'Y' in df.columns:
                df.rename(columns={'Y': 'Position_Y'}, inplace=True)

        # Add Assignment column if missing
        if 'Assignment' not in df.columns:
            df['Assignment'] = [f'Peak_{i+1}' for i in range(len(df))]

        return df

    def save_peak_list(self, df, file_path, format='csv'):
        """Save peak list to file"""
        try:
            if format.lower() == 'csv':
                df.to_csv(file_path, index=False, float_format='%.6f')
            else:
                df.to_csv(file_path, sep='\t', index=False, float_format='%.6f')

            self.add_recent_file(file_path)
            return True

        except Exception as e:
            raise RuntimeError(f"Failed to save peak list: {str(e)}")

    def add_recent_file(self, file_path):
        """Add file to recent files list"""
        abs_path = os.path.abspath(file_path)

        # Remove if already in list
        if abs_path in self.recent_files:
            self.recent_files.remove(abs_path)

        # Add to beginning
        self.recent_files.insert(0, abs_path)

        # Limit size
        if len(self.recent_files) > self.max_recent:
            self.recent_files = self.recent_files[:self.max_recent]

    def get_recent_files(self):
        """Get list of recent files"""
        # Filter out non-existent files
        existing_files = [f for f in self.recent_files if os.path.exists(f)]
        self.recent_files = existing_files
        return existing_files.copy()

    def find_files_in_folder(self, folder_path, file_types=None):
        """Find files of specified types in folder"""
        if not os.path.exists(folder_path):
            return []

        if file_types is None:
            file_types = self.supported_nmr_formats + self.supported_peak_formats

        files = []
        for file_type in file_types:
            pattern = os.path.join(folder_path, f"*.{file_type}")
            files.extend(glob.glob(pattern))

        return sorted(files)

    def get_file_metadata(self, file_path):
        """Get detailed file metadata"""
        if not os.path.exists(file_path):
            return None

        try:
            stat_info = os.stat(file_path)

            metadata = {
                'path': os.path.abspath(file_path),
                'name': os.path.basename(file_path),
                'extension': os.path.splitext(file_path)[1].lower().lstrip('.'),
                'size_bytes': stat_info.st_size,
                'size_mb': stat_info.st_size / (1024 * 1024),
                'created': datetime.fromtimestamp(stat_info.st_ctime),
                'modified': datetime.fromtimestamp(stat_info.st_mtime),
                'accessed': datetime.fromtimestamp(stat_info.st_atime),
                'is_nmr_file': stat_info.st_size > 1000  # Basic heuristic
            }

            # Add validation info
            if metadata['extension'] in self.supported_nmr_formats:
                valid, message = self.validate_nmr_file(file_path)
                metadata['valid'] = valid
                metadata['validation_message'] = message
            elif metadata['extension'] in self.supported_peak_formats:
                valid, message = self.validate_peak_file(file_path)
                metadata['valid'] = valid
                metadata['validation_message'] = message
            else:
                metadata['valid'] = False
                metadata['validation_message'] = "Unsupported file type"

            return metadata

        except Exception as e:
            return {
                'path': file_path,
                'error': str(e),
                'valid': False
            }

class DataValidator:
    """Validation utilities for NMR data integrity"""

    def __init__(self):
        self.validation_rules = {
            'chemical_shift_range_1H': (-5.0, 20.0),
            'chemical_shift_range_15N': (80.0, 140.0),
            'chemical_shift_range_13C': (0.0, 220.0),
            'min_peaks': 1,
            'max_peaks': 10000
        }

    def validate_chemical_shifts(self, df, nucleus_x='1H', nucleus_y='15N'):
        """Validate chemical shift ranges"""
        issues = []

        # Get appropriate ranges
        if nucleus_x == '1H':
            x_range = self.validation_rules['chemical_shift_range_1H']
        elif nucleus_x == '13C':
            x_range = self.validation_rules['chemical_shift_range_13C']
        else:
            x_range = (-50, 300)  # Generous range

        if nucleus_y == '15N':
            y_range = self.validation_rules['chemical_shift_range_15N']
        elif nucleus_y == '13C':
            y_range = self.validation_rules['chemical_shift_range_13C']
        else:
            y_range = (-50, 300)  # Generous range

        # Check X positions
        x_out_of_range = df[(df['Position_X'] < x_range[0]) | (df['Position_X'] > x_range[1])]
        if not x_out_of_range.empty:
            issues.append(f"{len(x_out_of_range)} peaks have {nucleus_x} shifts outside {x_range}")

        # Check Y positions
        y_out_of_range = df[(df['Position_Y'] < y_range[0]) | (df['Position_Y'] > y_range[1])]
        if not y_out_of_range.empty:
            issues.append(f"{len(y_out_of_range)} peaks have {nucleus_y} shifts outside {y_range}")

        return issues

    def validate_peak_list_integrity(self, df):
        """Comprehensive peak list validation"""
        issues = []

        # Check required columns
        required_cols = ['Position_X', 'Position_Y']
        missing_cols = [col for col in required_cols if col not in df.columns]
        if missing_cols:
            issues.append(f"Missing required columns: {missing_cols}")
            return issues  # Can't continue without these

        # Check for empty data
        if len(df) == 0:
            issues.append("Peak list is empty")
            return issues

        # Check for NaN values
        x_nan = df['Position_X'].isna().sum()
        y_nan = df['Position_Y'].isna().sum()
        if x_nan > 0:
            issues.append(f"{x_nan} peaks have missing X positions")
        if y_nan > 0:
            issues.append(f"{y_nan} peaks have missing Y positions")

        # Check for duplicates
        duplicates = df.duplicated(subset=['Position_X', 'Position_Y']).sum()
        if duplicates > 0:
            issues.append(f"{duplicates} duplicate peak positions found")

        # Check peak count
        peak_count = len(df)
        if peak_count < self.validation_rules['min_peaks']:
            issues.append(f"Too few peaks: {peak_count}")
        elif peak_count > self.validation_rules['max_peaks']:
            issues.append(f"Too many peaks: {peak_count}")

        return issues

    def auto_fix_peak_list(self, df):
        """Attempt to automatically fix common peak list issues"""
        fixed_df = df.copy()
        fixes_applied = []

        # Remove rows with NaN coordinates
        initial_count = len(fixed_df)
        fixed_df = fixed_df.dropna(subset=['Position_X', 'Position_Y'])
        removed_nan = initial_count - len(fixed_df)
        if removed_nan > 0:
            fixes_applied.append(f"Removed {removed_nan} peaks with missing coordinates")

        # Remove duplicate positions
        initial_count = len(fixed_df)
        fixed_df = fixed_df.drop_duplicates(subset=['Position_X', 'Position_Y'])
        removed_dupes = initial_count - len(fixed_df)
        if removed_dupes > 0:
            fixes_applied.append(f"Removed {removed_dupes} duplicate peak positions")

        # Ensure Assignment column exists
        if 'Assignment' not in fixed_df.columns:
            fixed_df['Assignment'] = [f'Peak_{i+1}' for i in range(len(fixed_df))]
            fixes_applied.append("Added Assignment column")

        # Reset index
        fixed_df = fixed_df.reset_index(drop=True)

        return fixed_df, fixes_applied

class FileMetadata:
    """File metadata collection and analysis"""

    def __init__(self):
        self.metadata_cache = {}
        self.cache_duration = 300  # 5 minutes

    def get_folder_summary(self, folder_path):
        """Get comprehensive folder summary"""
        if not os.path.exists(folder_path):
            return None

        try:
            summary = {
                'path': os.path.abspath(folder_path),
                'name': os.path.basename(folder_path),
                'total_files': 0,
                'total_size_mb': 0,
                'file_types': {},
                'nmr_files': [],
                'peak_files': [],
                'other_files': [],
                'last_modified': None,
                'scan_time': datetime.now()
            }

            file_manager = NMRFileManager()
            latest_mod_time = None

            for item in os.listdir(folder_path):
                item_path = os.path.join(folder_path, item)

                if os.path.isfile(item_path):
                    summary['total_files'] += 1

                    # Get file info
                    stat_info = os.stat(item_path)
                    size_mb = stat_info.st_size / (1024 * 1024)
                    summary['total_size_mb'] += size_mb

                    # Track latest modification
                    mod_time = datetime.fromtimestamp(stat_info.st_mtime)
                    if latest_mod_time is None or mod_time > latest_mod_time:
                        latest_mod_time = mod_time

                    # Categorize file
                    ext = os.path.splitext(item)[1].lower().lstrip('.')
                    summary['file_types'][ext] = summary['file_types'].get(ext, 0) + 1

                    if ext in file_manager.supported_nmr_formats:
                        summary['nmr_files'].append(item)
                    elif ext in file_manager.supported_peak_formats:
                        summary['peak_files'].append(item)
                    else:
                        summary['other_files'].append(item)

            summary['last_modified'] = latest_mod_time
            summary['nmr_file_count'] = len(summary['nmr_files'])
            summary['peak_file_count'] = len(summary['peak_files'])

            return summary

        except Exception as e:
            return {
                'path': folder_path,
                'error': str(e),
                'scan_time': datetime.now()
            }

    def compare_folders(self, folder1, folder2):
        """Compare two folders for compatibility"""
        summary1 = self.get_folder_summary(folder1)
        summary2 = self.get_folder_summary(folder2)

        if not summary1 or not summary2:
            return None

        comparison = {
            'folder1': summary1,
            'folder2': summary2,
            'compatibility': {},
            'recommendations': []
        }

        # Compare file counts
        if summary1['nmr_file_count'] != summary2['peak_file_count']:
            comparison['recommendations'].append(
                f"File count mismatch: {summary1['nmr_file_count']} NMR files vs "
                f"{summary2['peak_file_count']} peak files"
            )

        # Check for matching names
        nmr_names = [os.path.splitext(f)[0] for f in summary1['nmr_files']]
        peak_names = [os.path.splitext(f)[0] for f in summary2['peak_files']]

        matches = set(nmr_names) & set(peak_names)
        comparison['compatibility']['matching_names'] = len(matches)
        comparison['compatibility']['match_ratio'] = len(matches) / max(len(nmr_names), len(peak_names), 1)

        if comparison['compatibility']['match_ratio'] < 0.5:
            comparison['recommendations'].append("Low name matching - verify file correspondence")

        return comparison

import re
import pandas as pd
import glob
import os

"""
Transforms NMR peak list files in data/ matching 'T1_A1_WT_0o0*.txt'.
- Reads tab-separated files with headers: Assignment, Position_X, Position_Y, Height.
- Transforms assignments like '(<NA:A.10.VAL.H>, <NA:A.10.VAL.N>)' to '10.VAL'.
- Sorts by the number before the first '.' in the transformed Assignment (e.g., 10 in '10.VAL').
- Saves as comma-separated files with spaces after commas, including headers.
"""

# Define input directory and file pattern
input_directory = "data"
file_pattern = "*"  # Matches T1_A1_WT_0o0s, T1_A1_WT_0o0, etc.
input_files = glob.glob(os.path.join(cdinput_directory, file_pattern))

# Check if input directory exists and contains files
if not os.path.exists(input_directory):
    raise FileNotFoundError(f"Input directory '{input_directory}' does not exist")
if not input_files:
    print(f"No files found matching '{file_pattern}' in '{input_directory}'")
else:
    for file_path in input_files:
        try:
            # Read the file, recognizing the header
            data = pd.read_csv(file_path, sep="\t", header=0)
            
            # Verify expected columns
            expected_columns = ["Assignment", "Position_X", "Position_Y", "Height"]
            if not all(col in data.columns for col in expected_columns):
                print(f"Warning: '{file_path}' does not have expected columns, skipping")
                continue
            
            # Function to transform assignments
            def transform_assignment(entry):
                # Handle None or non-string entries
                if pd.isna(entry) or not isinstance(entry, str):
                    return "None"
                # Extract first <NA:A.<number>.<residue>.<atom>
                match = re.search(r"A\.(\d+)\.(\w+)\.\w+", entry)
                if match:
                    return f"{match.group(1)}.{match.group(2)}"
                # Handle cases like 'None' or unmatched formats
                return "None"
            
            # Function to extract number before first '.' for sorting
            def extract_sort_key(assignment):
                # Handle None or non-string entries
                if pd.isna(assignment) or not isinstance(assignment, str):
                    return float('inf')
                # Extract number before first '.'
                match = re.search(r'^(\d+)\.', assignment)
                return int(match.group(1)) if match else float('inf')
            
            # Apply transformation
            data['Assignment'] = data['Assignment'].apply(transform_assignment)
            
            # Create sort key from transformed Assignment
            data['SortKey'] = data['Assignment'].apply(extract_sort_key)
            
            # Sort by SortKey and drop it
            data = data.sort_values(by='SortKey', ascending=True).drop(columns=['SortKey'])
            
            # Generate output path
            output_path = file_path.replace('.txt', '_transformed.txt')
            
            # Save as comma-separated, then post-process to add spaces
            temp_path = output_path + '.tmp'
            data.to_csv(temp_path, sep=",", index=False)
            
            # Read temporary file and add spaces after commas
            with open(temp_path, 'r') as temp_file, open(output_path, 'w') as out_file:
                for line in temp_file:
                    # Replace commas with comma+space, avoiding trailing space
                    formatted_line = ', '.join(field.strip() for field in line.strip().split(','))
                    out_file.write(formatted_line + '\n')
            
            # Remove temporary file
            os.remove(temp_path)
            
            print(f"Transformed data saved to '{output_path}'")
        
        except pd.errors.ParserError:
            print(f"Warning: Failed to parse '{file_path}' (invalid format), skipping")
        except PermissionError:
            print(f"Warning: Cannot write to '{output_path}' (permission denied), skipping")
        except Exception as e:
            print(f"Warning: Error processing '{file_path}': {str(e)}, skipping")

print("Processing complete!")
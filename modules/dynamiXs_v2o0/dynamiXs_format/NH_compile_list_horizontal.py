import glob
import os
import re

# Define input directory and file pattern
input_directory = "data_comp"
#file_pattern = "*_peak_list_2_transformed.txt"

file_pattern = "*"

input_files = glob.glob(os.path.join(input_directory, file_pattern))

# Output file
output_file = "compiled_peak_lists.txt"

# Check if input directory exists and contains files
if not os.path.exists(input_directory):
    raise FileNotFoundError(f"Input directory '{input_directory}' does not exist")

if not input_files:
    print(f"No files found matching '{file_pattern}' in '{input_directory}'")
else:
    # Function to extract the number before "_peak" from filename
    def extract_number(filename):
        # Extract just the filename without path
        basename = os.path.basename(filename)
        # Look for number before "_peak"
        match = re.search(r'(\d+)_peak', basename)
        return int(match.group(1)) if match else float('inf')
    
    # Sort files by the number before "_peak"
    sorted_files = sorted(input_files, key=extract_number)
    
    # Read all files and store their content
    file_data = {}
    max_lines = 0
    
    for file_path in sorted_files:
        filename = os.path.basename(file_path)
        with open(file_path, 'r') as infile:
            lines = infile.readlines()
            # Strip newlines but keep the content
            lines = [line.rstrip('\n') for line in lines]
            file_data[filename] = lines
            max_lines = max(max_lines, len(lines))
    
    # Write the compiled output
    with open(os.path.join(input_directory, output_file), 'w') as outfile:
        # Write headers (filenames)
        headers = list(file_data.keys())
        header_line = "\t,\t,\t,\t,\t".join(headers)  # Tab-comma-tab for separation
        outfile.write(header_line + "\n")
        
        # Write data rows
        for row_idx in range(max_lines):
            row_data = []
            for filename in headers:
                if row_idx < len(file_data[filename]):
                    row_data.append(file_data[filename][row_idx])
                else:
                    row_data.append("")  # Empty cell if file has fewer lines
            
            # Join with tab-comma-tab for separation
            row_line = "\t,\t,\t,\t,\t".join(row_data)
            outfile.write(row_line + "\n")
    
    print(f"Compiled data saved to '{os.path.join(input_directory, output_file)}'")
    print(f"Files processed in order:")
    for i, file_path in enumerate(sorted_files, 1):
        filename = os.path.basename(file_path)
        number = extract_number(file_path)
        print(f"  {i}. {filename} (number: {number})")

print("Processing complete!")
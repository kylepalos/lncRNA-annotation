#!/bin/bash
touch paired_end_library_type.txt
output_file="paired_end_library_type.txt"

# Clear the output file if it exists
> "$output_file"

# Specify the directory with all of the salmon quant files
base_directory="/home/kpalos/new_arabidopsis_lncRNAs/salmon_files/paired_end_quants"

# Specify the name of the file from which you want to extract a line
file_name="salmon_quant.log"

# Specify the line number you want to extract
target_line_number=12

# Iterate through each subdirectory
for directory in "$base_directory"/*/; do
    # Get the directory name by removing the path
    dir_name=$(basename "$directory")

    # Print the directory name
    echo -n "$dir_name" >> "$output_file"
    # Construct the path to the target file
    target_file_path="$directory/logs/$file_name"

    # Check if the target file exists
    if [ -f "$target_file_path" ]; then
        # Extract the specified line from the target file
        specific_line=$(sed -n "${target_line_number}p" "$target_file_path")
        last_three_chars="${specific_line: -3}"

        # Print the extracted line
        echo -e ' \t ' "$last_three_chars" >> "$output_file"
    else
        echo "Target file not found in $directory"
    fi

done

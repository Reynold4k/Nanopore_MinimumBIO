#!/bin/bash

## Overview

#The main purpose of this script is to analyze the gene coverage between experimental and control groups, calculate differential coverage. The specific steps are as follows:

#1. **Load Required Modules**: Use the `module load` command to load the necessary Python and computational tools.
#2. **Define Paths**: Specify the paths for the experimental group, control group, and reference files.
#3. **Find the Latest Experimental Directory**: Retrieve the latest experimental and control round directories based on the folder structure.
#4. **Check Coverage Files**: Confirm that the coverage files exist.
#5. **Read Coverage Data**: Read the coverage data for the control and experimental groups into associative arrays.
#6. **Calculate Differences**: Compare the coverage between the control and experimental groups, calculate the difference, and write the results to an output file.




# What you what to modify:
# 1.Change the line 25 EXPERIMENTAL_FOLDER to your/exp/fastq parent path
# 2.Change the line 26 CONTROL_FOLDER to your/control/fastq parent path




# Define paths
EXPERIMENTAL_FOLDER="/srv/scratch/z3546698/tutorial/Bait_Glue/VHL/MB015/TON/230827"
CONTROL_FOLDER="/srv/scratch/z3546698/tutorial/Bait_Glue/VHL/No glue/TON/240427"


# Function to get all round directory names
function get_round_directories {
    local folder="$1"
    find "$folder" -mindepth 1 -maxdepth 1 -type d -name "R*"
}

# Get all round directories for experimental and control groups
round_dirs_exp=($(get_round_directories "$EXPERIMENTAL_FOLDER"))
round_dirs_ctrl=($(get_round_directories "$CONTROL_FOLDER"))

# Check if any round directories were found
if [ ${#round_dirs_exp[@]} -eq 0 ] || [ ${#round_dirs_ctrl[@]} -eq 0 ]; then
    echo "Error: No round directories found in one or both folders."
    exit 1
fi

# Loop through each round directory in the experimental group
for round_exp in "${round_dirs_exp[@]}"; do
    round_name=$(basename "$round_exp")

    # Corresponding control round directory based on the round name
    last_round_ctrl="${CONTROL_FOLDER}/${round_name}"

    # Check if control round directory exists
    if [ ! -d "$last_round_ctrl" ]; then
        echo "Warning: Control round directory $last_round_ctrl does not exist. Skipping."
        continue
    fi

    # Paths for top 1000 positions for experimental and control
    top_1000_exp="${round_exp}/step2/top_1000_positions.txt"
    top_1000_ctrl="${last_round_ctrl}/step2/top_1000_positions.txt"

    # Check if top 1000 files exist
    if [ ! -f "$top_1000_exp" ] || [ ! -f "$top_1000_ctrl" ]; then
        echo "Error: One or both top 1000 position files do not exist in $round_name. Skipping."
        continue
    fi

    # Declare associative arrays for coverage values
    declare -A control_coverage
    declare -A experimental_coverage

    # Read control coverage into an associative array
    while read -r gene_id start end coverage; do
        control_coverage["$gene_id"]="$coverage"
    done < <(cut -f 1,2,3,4 "$top_1000_ctrl")  # Assuming format: gene_id start end coverage

    # Read experimental coverage for top 1000 genes
    while read -r gene_id start end exp_coverage; do
        experimental_coverage["$gene_id"]=$(echo "${experimental_coverage["$gene_id"]:-0} + $exp_coverage" | bc)
    done < <(cut -f 1,2,3,4 "$top_1000_exp")  # Assuming format: gene_id start end coverage

    # Output file for differential coverage
    output_diff_coverage="${round_exp}/differential_coverage.txt"

    # Header for the output file
    echo -e "Gene\tControl_Coverage\tExperimental_Coverage\tDifference" > "$output_diff_coverage"

    # Write the results, summing the control and experimental coverage for the same gene ID
    for gene_id in "${!experimental_coverage[@]}"; do
        ctrl_coverage=${control_coverage[$gene_id]:-0}  # Default to 0 if not found in control
        exp_coverage=${experimental_coverage[$gene_id]:-0}  # Get the summed experimental coverage
        
        # Calculate the difference
        difference=$(awk "BEGIN {print $exp_coverage - $ctrl_coverage}")

        # Only output if the difference is non-zero
        if (( $(echo "$difference != 0" | bc -l) )); then
            # Output to file
            echo -e "$gene_id\t$ctrl_coverage\t$exp_coverage\t$difference" >> "$output_diff_coverage"
        fi
    done

    # Sort the output file based on the difference
    sort -k4,4gr -S 75% --parallel=6 --temporary-directory=/tmp -o "$output_diff_coverage" "$output_diff_coverage"

    echo "Differential coverage analysis for $round_name complete. Results saved in: $output_diff_coverage"
done


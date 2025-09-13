#!/bin/bash

# Help function to display usage instructions
show_help() {
    cat << EOF
Usage: $(basename "$0") [OPTIONS]

This script traverses LRT and Wald directories, checking for significant.tsv files
and executes a pipeline command when found.

Required arguments:
    -p, --pipeline PATH    Path to the pipeline directory containing fetch_1.sh
    -i, --input PATH      Path to the input data directory

Optional arguments:
    -t, --threads NUM     Number of threads to use (default: 40)
    -h, --help           Show this help message and exit

Example:
    $(basename "$0") --pipeline /path/to/pipeline --input /path/to/input
EOF
}

# Initialize default values
threads=40
pipeline_path=""
input_path=""

# Parse command line arguments using a while loop
while [[ $# -gt 0 ]]; do
    case $1 in
        -p|--pipeline)
            pipeline_path="$2"
            shift 2
            ;;
        -i|--input)
            input_path="$2"
            shift 2
            ;;
        -t|--threads)
            threads="$2"
            shift 2
            ;;
        -h|--help)
            show_help
            exit 0
            ;;
        *)
            echo "Error: Unknown argument $1"
            show_help
            exit 1
            ;;
    esac
done

# Validate required arguments
if [ -z "$pipeline_path" ] || [ -z "$input_path" ]; then
    echo "Error: Missing required arguments"
    show_help
    exit 1
fi

# Validate that pipeline directory exists and contains fetch_1.sh
if [ ! -d "$pipeline_path" ]; then
    echo "Error: Pipeline directory does not exist: $pipeline_path"
    exit 1
fi

if [ ! -f "$pipeline_path/fetch_1.sh" ]; then
    echo "Error: fetch_1.sh not found in pipeline directory: $pipeline_path"
    exit 1
fi

# Validate that input path exists
if [ ! -d "$input_path" ]; then
    echo "Error: Input directory does not exist: $input_path"
    exit 1
fi

# Store the starting directory
starting_dir=$(pwd)

# Define the main test types and result categories as arrays for easy iteration
test_types=("LRT" "Wald")
result_categories=("positive" "negative")

# Print execution settings for confirmation
echo "Execution settings:"
echo "Pipeline path: $pipeline_path"
echo "Input path: $input_path"
echo "Number of threads: $threads"
echo "Starting directory: $starting_dir"
echo "------------------------"

# Loop through each test type (LRT and Wald)
for test in "${test_types[@]}"; do
    # Loop through each result category (positive and negative)
    for result in "${result_categories[@]}"; do
        # Construct the full directory path for current iteration
        dir_path="${test}/${result}"
        echo "Checking directory: ${dir_path}"
        # Check if significant.tsv exists in current directory
        if [ -f "${dir_path}/significant.tsv" ]; then
            echo "Found significant.tsv in ${dir_path}"
            
            # Get file size for logging purposes
            size=$(ls -lh "${dir_path}/significant.tsv" | awk '{print $5}')
            echo "File size: ${size}"
            
            # Change to the directory containing significant.tsv
            echo "Changing to directory: ${dir_path}"
            if ! pushd "${dir_path}" &>/dev/null; then
                echo "Error: Failed to change to directory: ${dir_path}"
                continue
            fi
            
            # Get absolute path of significant.tsv
            significant_abs_path="$(pwd)/significant.tsv"
            
            # Execute the pipeline command
            echo "Executing pipeline..."
            if ! bash "$pipeline_path/fetch_1.sh" \
                --significant "$significant_abs_path" \
                --data "$input_path" \
                --threads "$threads"; then
                echo "Warning: Pipeline execution failed for ${dir_path}"
            fi
            # Return to the previous directory
            echo "Returning to previous directory"
            popd &>/dev/null
            echo "Pipeline execution completed for ${dir_path}"
        else
            echo "No significant.tsv found in ${dir_path}"
        fi
        echo "------------------------"
    done
done

# Extra safety check to ensure we're back in the starting directory
if [ "$(pwd)" != "$starting_dir" ]; then
    echo "Warning: Directory mismatch detected, returning to starting directory"
    cd "$starting_dir" || exit 1
fi

#!/bin/bash
# Function to show usage/help
usage() {
    echo "Usage: $0 -b <base_path> -i <input_path> -p <pipeline_path> -m <metadata_file> -s <strand> -t <num_threads>"
    echo "  -b: Base path for the project"
    echo "  -i: Path to input directory containing processed fastq files"
    echo "  -p: Path to pipeline scripts"
    echo "  -m: Path to metadata TSV file"
    echo "  -s: Strand (default: 'forward')"
    echo "  -t: Number of CPU threads to use"
    exit 1
}

# Parse command line arguments with defaults
strand="forward"
while getopts "b:i:p:m:s:t:h" opt; do
    case $opt in
        b) base_path="$OPTARG" ;;
        i) input_path="$OPTARG" ;;
        p) pipeline="$OPTARG" ;;
        m) metadata="$OPTARG" ;;
        s) strand="$OPTARG" ;;
        t) num_cpus="$OPTARG" ;;
        h) usage ;;
        ?) usage ;;
    esac
done

# Check if required arguments are provided
if [ -z "$base_path" ] || [ -z "$input_path" ] || [ -z "$pipeline" ] || [ -z "$metadata" ] || [ -z "$num_cpus" ]; then
    usage
fi

# Convert paths to absolute paths
base_path=$(realpath "$base_path")
input_path=$(realpath "$input_path")
pipeline=$(realpath "$pipeline")
metadata=$(realpath "$metadata")

# Create CGR directory
CGR_path="${base_path}/02_CGR"
mkdir -p "$CGR_path"

# Create logs directory
logs_dir="${CGR_path}/cgr_logs"
mkdir -p "$logs_dir"

# Ensure the logs directory exists and is writable
if [[ ! -d "$logs_dir" ]]; then
    echo "Error: Could not create logs directory: $logs_dir"
    exit 1
fi

# Count total number of samples for progress bar
total_samples=$(awk 'NR>1' "$metadata" | wc -l)
current_sample=0

# Process each sample
while IFS=$'\t' read -r SRA description sample strain clone; do
    # Skip header
    [ "$SRA" = "SRA" ] && continue

    # Create progress bar
    current_sample=$((current_sample + 1))
    progress=$((current_sample * 100 / total_samples))
    progress_bar=$(printf "[%-50s] %d%%\r" $(printf "#%.0s" $(seq 1 $((progress/2)))) $progress)
    echo -ne "$progress_bar"

    # Create sample directory
    sample_dir="${CGR_path}/${sample}"
    mkdir -p "$sample_dir"
    
    # Create log file early to check for permissions
    log_file="${logs_dir}/${sample}_cgr.log"
    touch "$log_file" || {
        echo "Error: Cannot create log file: $log_file"
        exit 1
    }

    cd "$sample_dir" || exit 1

    # Run CGR generation with redirected output
    python "${pipeline}/generate_CGR_1.py" \
        --merge_inpre "${input_path}/${sample}/${description}_merge" \
        --r1_inpre "${input_path}/${sample}/${description}_R1" \
        --r2_inpre "${input_path}/${sample}/${description}_R2" \
        --outpre "$sample" \
        --strand "$strand" \
        --num_cpus "$num_cpus" \
        2>&1 > "$log_file"

    cd "$base_path" || exit 1
done < "$metadata"

echo -e "\nCGR generation complete. Logs are stored in: $logs_dir"

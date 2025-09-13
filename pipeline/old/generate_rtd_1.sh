#!/bin/bash
# Function to show usage/help
usage() {
    echo "Usage: $0 -b <base_path> -i <input_path> -p <pipeline_path> -m <metadata_file> -k <kmer_length> -a <alphabet> -s <strand> -t <num_threads>"
    echo "  -b: Base path for the project"
    echo "  -i: Path to input directory containing processed fastq files"
    echo "  -p: Path to pipeline scripts"
    echo "  -m: Path to metadata TSV file"
    echo "  -k: K-mer length (default: 2)"
    echo "  -a: Alphabet (default: A C G T)"
    echo "  -s: Strand (default: 'forward')"
    echo "  -t: Number of CPU threads to use"
    exit 1
}

# Parse command line arguments with defaults
k_len=2
alphabet=("A" "C" "G" "T")  # Store alphabet as array
strand="forward"
while getopts "b:i:p:m:k:a:s:t:h" opt; do
    case $opt in
        b) base_path="$OPTARG" ;;
        i) input_path="$OPTARG" ;;
        p) pipeline="$OPTARG" ;;
        m) metadata="$OPTARG" ;;
        k) k_len="$OPTARG" ;;
        a) # Split alphabet input into array
           IFS=' ' read -r -a alphabet <<< "$OPTARG" ;;
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

# Create RTD directory
RTD_path="${base_path}/01_RTD"
mkdir -p "$RTD_path"

# Create logs directory
logs_dir="${RTD_path}/rtd_logs"
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
    sample_dir="${RTD_path}/${sample}"
    mkdir -p "$sample_dir"
    
    # Create log file early to check for permissions
    log_file="${logs_dir}/${sample}_rtd.log"
    touch "$log_file" || {
        echo "Error: Cannot create log file: $log_file"
        exit 1
    }

    cd "$sample_dir" || exit 1

    # Run RTD generation with redirected output
    # Construct alphabet argument by joining array elements with spaces
    alphabet_arg=""
    for letter in "${alphabet[@]}"; do
        alphabet_arg+="$letter "
    done

    python "${pipeline}/generate_RTD_1.py" \
        --k "$k_len" \
        --alphabet ${alphabet_arg} \
        --merge_inpre "${input_path}/${sample}/${description}_merge" \
        --r1_inpre "${input_path}/${sample}/${description}_R1" \
        --r2_inpre "${input_path}/${sample}/${description}_R2" \
        --outpre "$sample" \
        --strand "$strand" \
        --num_cpus "$num_cpus" \
        2>&1 > "$log_file"

    cd "$base_path" || exit 1
done < "$metadata"

echo -e "\nRTD generation complete. Logs are stored in: $logs_dir"

#!/bin/bash
# Function to show usage/help
usage() {
    echo "Usage: $0 -b <base_path> -i <input_path> -p <pipeline_path> -m <metadata_file> -k <kmer_length> -a <alphabet> -s <strand> -t <num_threads> -o <output_dir> [--paired TRUE/FALSE]"
    echo "  -b: Base path for the project"
    echo "  -i: Path to input directory containing processed fastq files"
    echo "  -p: Path to pipeline scripts"
    echo "  -m: Path to metadata TSV file"
    echo "  -k: K-mer length (default: 2)"
    echo "  -a: Alphabet (default: A C G T)"
    echo "  -s: Strand (default: 'forward')"
    echo "  -t: Number of CPU threads to use"
    echo "  -o: Output directory name (default: '01_RTD')"
    echo "  --paired: Whether to process paired-end reads (default: TRUE)"
    exit 1
}

# Parse command line arguments with defaults
k_len=2
alphabet=("A" "C" "G" "T")  # Store alphabet as array
strand="forward"
output_dir="01_RTD"  # Default output directory name
paired="TRUE"  # Default to paired-end

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -b) base_path="$2"; shift 2 ;;
        -i) input_path="$2"; shift 2 ;;
        -p) pipeline="$2"; shift 2 ;;
        -m) metadata="$2"; shift 2 ;;
        -k) k_len="$2"; shift 2 ;;
        -a) IFS=' ' read -r -a alphabet <<< "$2"; shift 2 ;;
        -s) strand="$2"; shift 2 ;;
        -t) num_cpus="$2"; shift 2 ;;
        -o) output_dir="$2"; shift 2 ;;
        --paired) paired="$2"; shift 2 ;;
        -h) usage ;;
        *) usage ;;
    esac
done

# Convert paired to uppercase for consistency
paired=${paired^^}

# Validate paired argument
if [[ "$paired" != "TRUE" && "$paired" != "FALSE" ]]; then
    echo "Error: --paired must be TRUE or FALSE"
    exit 1
fi

# Check if required arguments are provided
if [ -z "$base_path" ] || [ -z "$input_path" ] || [ -z "$pipeline" ] || [ -z "$metadata" ] || [ -z "$num_cpus" ]; then
    usage
fi

# Convert paths to absolute paths
base_path=$(realpath "$base_path")
input_path=$(realpath "$input_path")
pipeline=$(realpath "$pipeline")
metadata=$(realpath "$metadata")

# Create RTD directory using the specified output directory name
RTD_path="${base_path}/${output_dir}"
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

    # Construct alphabet argument by joining array elements with spaces
    alphabet_arg=""
    for letter in "${alphabet[@]}"; do
        alphabet_arg+="$letter "
    done

    # Set r1_inpre based on whether data is paired or single-ended
    if [[ "$paired" == "TRUE" ]]; then
        r1_path="\"${input_path}/${sample}/${description}_R1\""
    else
        r1_path="\"${input_path}/${sample}/${description}\""
    fi

    # Construct the base command
    cmd="python \"${pipeline}/generate_RTD.py\" \
        --k \"$k_len\" \
        --alphabet ${alphabet_arg} \
        --r1_inpre ${r1_path} \
        --outpre \"$sample\" \
        --strand \"$strand\" \
        --num_cpus \"$num_cpus\" \
        --paired $paired"

    # Add paired-end specific arguments if paired is TRUE
    if [[ "$paired" == "TRUE" ]]; then
        cmd+=" --r2_inpre \"${input_path}/${sample}/${description}_R2\" \
              --merge_inpre \"${input_path}/${sample}/${description}_merge\""
    fi

    # Execute the command with redirected output
    eval "$cmd" 2>&1 > "$log_file"

    cd "$base_path" || exit 1
done < "$metadata"

echo -e "\nRTD generation complete. Logs are stored in: $logs_dir"

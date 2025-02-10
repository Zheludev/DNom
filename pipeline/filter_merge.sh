#!/bin/bash

# Function to show usage/help
usage() {
    echo "Usage: $0 -p <input_path> -m <metadata_file> -t <num_threads> [--paired TRUE|FALSE]"
    echo "  -p: Path to input directory containing fastq files"
    echo "  -m: Path to metadata TSV file"
    echo "  -t: Number of CPU threads to use"
    echo "  --paired: Whether the input is paired-end (TRUE) or single-end (FALSE) data [default: TRUE]"
    exit 1
}

# Parse command line arguments
paired="TRUE"  # Default value
while [ "$#" -gt 0 ]; do
    case "$1" in
        -p) input_path="$2"; shift 2 ;;
        -m) metadata="$2"; shift 2 ;;
        -t) num_cpus="$2"; shift 2 ;;
        --paired) paired="$2"; shift 2 ;;
        -h) usage ;;
        *) usage ;;
    esac
done

# Check if required arguments are provided
if [ -z "$input_path" ] || [ -z "$metadata" ] || [ -z "$num_cpus" ]; then
    usage
fi

# Validate paired argument
if [ "$paired" != "TRUE" ] && [ "$paired" != "FALSE" ]; then
    echo "Error: --paired must be either TRUE or FALSE"
    exit 1
fi

input_path=$(realpath "$input_path")
metadata=$(realpath "$metadata")
logs_dir="${input_path}/fastp_logs"

# Create logs directory
logs_dir="${input_path}/fastp_logs"
mkdir -p "$logs_dir"

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
    mkdir -p "${input_path}/${sample}"
    cd "${input_path}/${sample}" || exit 1
    
    # Run fastp based on paired-end status
    if [ "$paired" = "TRUE" ]; then
        # Paired-end processing
        fastp --in1="${input_path}/${SRA}/${SRA}_1.fastq.gz" \
              --out1="${description}_R1.fastq.gz" \
              --in2="${input_path}/${SRA}/${SRA}_2.fastq.gz" \
              --out2="${description}_R2.fastq.gz" \
              --merged_out="${description}_merge.fastq.gz" \
              --merge \
              --average_qual=30 \
              --cut_front \
              --cut_tail \
              --thread "$num_cpus" \
              --json "${logs_dir}/${sample}_fastp.json" \
              --html "${logs_dir}/${sample}_fastp.html" \
              2>&1 > "${logs_dir}/${sample}_fastp.log"
    else
        # Single-end processing
        fastp --in1="${input_path}/${SRA}/${SRA}.fastq.gz" \
              --out1="${description}.fastq.gz" \
              --average_qual=30 \
              --cut_front \
              --cut_tail \
              --thread "$num_cpus" \
              --json "${logs_dir}/${sample}_fastp.json" \
              --html "${logs_dir}/${sample}_fastp.html" \
              2>&1 > "${logs_dir}/${sample}_fastp.log"
    fi
    
    cd "$input_path" || exit 1
    
done < "$metadata"

echo -e "\nMoving raw read directories to raw_reads/..."

# Create raw_reads directory and move SRR directories
raw_reads_dir="${input_path}/raw_reads"
mkdir -p "$raw_reads_dir"

# Move all SRR directories using the metadata file
while IFS=$'\t' read -r SRA description sample strain clone; do
    # Skip header
    [ "$SRA" = "SRA" ] && continue
    
    # Move SRR directory if it exists
    if [ -d "${input_path}/${SRA}" ]; then
        mv "${input_path}/${SRA}" "$raw_reads_dir/"
    fi
done < "$metadata"

echo "Processing complete. Logs are stored in: $logs_dir"
echo "Raw reads moved to: $raw_reads_dir"

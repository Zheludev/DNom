#!/bin/bash

# Function to show usage/help
usage() {
    echo "Usage: $0 -p <input_path> -m <metadata_file> -t <num_threads>"
    echo "  -p: Path to input directory containing fastq files"
    echo "  -m: Path to metadata TSV file"
    echo "  -t: Number of CPU threads to use"
    exit 1
}

# Parse command line arguments
while getopts "p:m:t:h" opt; do
    case $opt in
        p) input_path="$OPTARG" ;;
        m) metadata="$OPTARG" ;;
        t) num_cpus="$OPTARG" ;;
        h) usage ;;
        ?) usage ;;
    esac
done

# Check if required arguments are provided
if [ -z "$input_path" ] || [ -z "$metadata" ] || [ -z "$num_cpus" ]; then
    usage
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
    
    # Run fastp with redirected output
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

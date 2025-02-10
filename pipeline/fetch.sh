#!/bin/bash

# Default values for arguments
significant=""
data_path=""
threads=4  # Default is 4 threads
paired="TRUE"  # Default is TRUE for paired-end data

# Function to parse arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --significant) significant="$2"; shift ;;
        --data) data_path="$2"; shift ;;
        --threads) threads="$2"; shift ;;
        --paired) paired="$2"; shift ;;
        *) echo "Unknown parameter: $1"; exit 1 ;;
    esac
    shift
done

# Check if required arguments are provided
if [[ -z "$significant" || -z "$data_path" ]]; then
    echo "Error: --significant and --data arguments are required."
    exit 1
fi

# Validate paired argument
if [[ "$paired" != "TRUE" && "$paired" != "FALSE" ]]; then
    echo "Error: --paired must be either TRUE or FALSE"
    exit 1
fi

# Ensure significant.tsv file exists
if [[ ! -f "$significant" ]]; then
    echo "Error: Significant file '$significant' does not exist."
    exit 1
fi

# Ensure data_path is valid
if [[ ! -d "$data_path" ]]; then
    echo "Error: Data path '$data_path' does not exist."
    exit 1
fi

# Divide the input threads by 4 (real_threads)
real_threads=$((threads / 4))

# Create 'reads' directory if it doesn't exist
mkdir -p reads
cd reads

# Function to process each cluster
process_cluster() {
    local cluster=$1
    local significant=$2
    local data_path=$3
    local paired=$4

    # Create cluster-specific directory
    mkdir -p "$cluster"
    cd "$cluster"

    if [[ "$paired" == "TRUE" ]]; then
        # Paired-end processing
        awk -F'\t' -v cluster="$cluster" '($3 == cluster) {print $1}' "$significant" | \
        seqkit grep --quiet -f - --infile-list <(awk -F'\t' -v cluster="$cluster" '($3 == cluster) {split($1,SRA,/\./); print SRA[1]}' ../../significant.tsv | uniq | grep -Ff - <(find "$data_path" -name "*_1.fastq.gz")) | pigz --best > "$cluster"_1.fastq.gz

        awk -F'\t' -v cluster="$cluster" '($3 == cluster) {print $1}' "$significant" | \
        seqkit grep --quiet -f - --infile-list <(awk -F'\t' -v cluster="$cluster" '($3 == cluster) {split($1,SRA,/\./); print SRA[1]}' ../../significant.tsv | uniq | grep -Ff - <(find "$data_path" -name "*_2.fastq.gz")) | pigz --best > "$cluster"_2.fastq.gz
    else
        # Single-end processing
        awk -F'\t' -v cluster="$cluster" '($3 == cluster) {print $1}' "$significant" | \
        seqkit grep --quiet -f - --infile-list <(awk -F'\t' -v cluster="$cluster" '($3 == cluster) {split($1,SRA,/\./); print SRA[1]}' ../../significant.tsv | uniq | grep -Ff - <(find "$data_path" -name "*.fastq.gz")) | pigz --best > "$cluster".fastq.gz
    fi

    cd ..
    echo ''
}

export -f process_cluster  # Export function to be used with parallel

# Count the total number of clusters
total_clusters=$(awk -F'\t' 'FNR>1 {print $3}' "$significant" | uniq | wc -l)

# Read unique clusters from significant.tsv and run them in parallel, pipe the progress to tqdm
awk -F'\t' 'FNR>1 {print $3}' "$significant" | uniq | \
    parallel -j "$real_threads" process_cluster {} "$significant" "$data_path" "$paired" | \
    tqdm --total="$total_clusters" --desc="Processing clusters" --unit="cluster" --unit_scale=1 > /dev/null

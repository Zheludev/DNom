#!/bin/bash

# Exit on error
set -e

# Function to display usage information
usage() {
    echo "Usage: $0 -p PIPELINE_DIR -e EMBEDDING_TYPE -m PARAMETER [OPTIONS]"
    echo ""
    echo "Required arguments:"
    echo "  -p PIPELINE_DIR     Path to pipeline directory containing scripts"
    echo "  -e EMBEDDING_TYPE   Type of embedding (RTD, CGR, or merge)"
    echo "  -m PARAMETER        Column in metadata.tsv for Sleuth analysis"
    echo ""
    echo "Optional arguments:"
    echo "  -d BASE_DIR         Base directory (default: current directory)"
    echo "  -t THREADS         Number of threads (default: 16)"
    echo "  -k K_SIZE          k size for RTD embedding (default: 2)"
    echo "  -a ALPHABET        Nucleotide alphabet for RTD (default: 'A C G T')"
    echo "  -s STRAND          Read orientation (default: forward)"
    echo "  -c MIN_CLUSTER     Minimum cluster size (default: 400)"
    echo "  -n MIN_SAMPLES     Minimum samples for HDBSCAN (default: 5)"
    echo "  -g MAX_GMM         Maximum GMMs to fit (default: 5)"
    echo "  -r SEED            Random seed (default: 42)"
    echo "  -o OVERLAP         Reads needed for overlap (default: 10)"
    echo "  -b N_BOOTSTRAPS    Number of bootstraps (default: 100)"
    echo "  -u PSEUDOCOUNT     Pseudocount for Sleuth (default: 0.5)"
    echo "  -q ALPHA           Significance threshold (default: 0.05)"
    echo "  -f FOLDCHANGE      Foldchange threshold (default: 10)"
    exit 1
}

# Function to log major steps
log_step() {
    echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"
}

# Function to validate required files and directories
validate_requirements() {
    local pipeline_dir=$1
    local base_dir=$2

    # Check pipeline directory and scripts
    if [[ ! -d "$pipeline_dir" ]]; then
        echo "Error: Pipeline directory not found: $pipeline_dir"
        exit 1
    fi

    required_scripts=(
        "filter_merge_1.sh"
        "generate_rtd_1.sh"
        "generate_cgr_1.sh"
        "cluster_1.py"
        "merge_1.py"
        "bootstrap_2.py"
        "sleuth_1.R"
        "extract_1.py"
        "fetch_loop_1.sh"
    )

    for script in "${required_scripts[@]}"; do
        if [[ ! -f "$pipeline_dir/$script" ]]; then
            echo "Error: Required script not found: $pipeline_dir/$script"
            exit 1
        fi
    done

    # Check metadata.tsv
    if [[ ! -f "$base_dir/metadata.tsv" ]]; then
        echo "Error: metadata.tsv not found in $base_dir"
        exit 1
    fi

    # Validate metadata.tsv columns
    required_columns=("SRA" "description" "sample")
    header=$(head -n 1 "$base_dir/metadata.tsv")
    for col in "${required_columns[@]}"; do
        if ! echo "$header" | grep -q "$col"; then
            echo "Error: Required column '$col' not found in metadata.tsv"
            exit 1
        fi
    done
}

# Parse command line arguments
while getopts "p:e:m:d:t:k:a:s:c:n:g:r:o:b:u:q:f:h" opt; do
    case $opt in
        p) PIPELINE_DIR="$OPTARG" ;;
        e) EMBEDDING_TYPE="$OPTARG" ;;
        m) PARAMETER="$OPTARG" ;;
        d) BASE_DIR="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
        k) K_SIZE="$OPTARG" ;;
        a) ALPHABET="$OPTARG" ;;
        s) STRAND="$OPTARG" ;;
        c) MIN_CLUSTER="$OPTARG" ;;
        n) MIN_SAMPLES="$OPTARG" ;;
        g) MAX_GMM="$OPTARG" ;;
        r) SEED="$OPTARG" ;;
        o) OVERLAP="$OPTARG" ;;
        b) N_BOOTSTRAPS="$OPTARG" ;;
        u) PSEUDOCOUNT="$OPTARG" ;;
        q) ALPHA="$OPTARG" ;;
        f) FOLDCHANGE="$OPTARG" ;;
        h) usage ;;
        \?) usage ;;
    esac
done

# Validate required arguments
if [[ -z "$PIPELINE_DIR" || -z "$EMBEDDING_TYPE" || -z "$PARAMETER" ]]; then
    echo "Error: Missing required arguments"
    usage
fi

# Validate embedding type
if [[ ! "$EMBEDDING_TYPE" =~ ^(RTD|CGR|merge)$ ]]; then
    echo "Error: Invalid embedding type. Must be RTD, CGR, or merge"
    exit 1
fi

# Set defaults for optional arguments
BASE_DIR="${BASE_DIR:-$(pwd)}"
THREADS="${THREADS:-16}"
K_SIZE="${K_SIZE:-2}"
ALPHABET="${ALPHABET:-'A C G T'}"
STRAND="${STRAND:-forward}"
MIN_CLUSTER="${MIN_CLUSTER:-400}"
MIN_SAMPLES="${MIN_SAMPLES:-5}"
MAX_GMM="${MAX_GMM:-5}"
SEED="${SEED:-42}"
OVERLAP="${OVERLAP:-10}"
N_BOOTSTRAPS="${N_BOOTSTRAPS:-100}"
PSEUDOCOUNT="${PSEUDOCOUNT:-0.5}"
ALPHA="${ALPHA:-0.05}"
FOLDCHANGE="${FOLDCHANGE:-10}"

# Validate requirements
validate_requirements "$PIPELINE_DIR" "$BASE_DIR"

# Setup directory structure
INPUT_PATH="$BASE_DIR/00_input"
RTD_PATH="$BASE_DIR/01_RTD"
CGR_PATH="$BASE_DIR/02_CGR"
RTD_CLUST_PATH="$BASE_DIR/03_RTD_clust"
CGR_CLUST_PATH="$BASE_DIR/04_CGR_clust"
MERGE_PATH="$BASE_DIR/05_merge"
RTD_BOOTSTRAP_PATH="$BASE_DIR/06_RTD_bootstrap"
CGR_BOOTSTRAP_PATH="$BASE_DIR/07_CGR_bootstrap"
MERGE_BOOTSTRAP_PATH="$BASE_DIR/08_merge_bootstrap"
RTD_SLEUTH_PATH="$BASE_DIR/09_RTD_Sleuth"
CGR_SLEUTH_PATH="$BASE_DIR/10_CGR_Sleuth"
MERGE_SLEUTH_PATH="$BASE_DIR/11_merge_Sleuth"
RTD_SIG_PATH="$BASE_DIR/12_RTD_significant"
CGR_SIG_PATH="$BASE_DIR/13_CGR_significant"
MERGE_SIG_PATH="$BASE_DIR/14_merge_significant"

# Function to run RTD analysis
run_rtd() {
    log_step "Starting RTD analysis"
    
    # Generate RTD embedding
    mkdir -p "$RTD_PATH"
    "$PIPELINE_DIR/generate_rtd_1.sh" -b "$BASE_DIR" -i "$INPUT_PATH" -p "$PIPELINE_DIR" \
        -m metadata.tsv -k "$K_SIZE" -a "$ALPHABET" -s "$STRAND" -t "$THREADS"
    
    # Run clustering
    mkdir -p "$RTD_CLUST_PATH"
    cd "$RTD_CLUST_PATH"
    python "$PIPELINE_DIR/cluster_1.py" --file_pattern "$RTD_PATH/*/*.RTD.gz" \
        --output_file clustering_results.tsv --min_cluster_size "$MIN_CLUSTER" \
        --min_samples "$MIN_SAMPLES" --num_cpus "$THREADS" --max_GMM "$MAX_GMM" --seed "$SEED"
    cd "$BASE_DIR"
    
    # Run bootstrapping
    mkdir -p "$RTD_BOOTSTRAP_PATH"
    cd "$RTD_BOOTSTRAP_PATH"
    python "$PIPELINE_DIR/bootstrap_2.py" --counts "$RTD_CLUST_PATH/table.otu" \
        --volumes "$RTD_CLUST_PATH/est_volumes.tsv" --n_bootstraps "$N_BOOTSTRAPS"
    cd "$BASE_DIR"
    
    # Run Sleuth
    mkdir -p "$RTD_SLEUTH_PATH"
    cd "$RTD_SLEUTH_PATH"
    Rscript "$PIPELINE_DIR/sleuth_1.R" --root "$BASE_DIR" --metadata metadata.tsv \
        --bootstrap "06_RTD_bootstrap" --pseudocount "$PSEUDOCOUNT" --alpha "$ALPHA" \
        --fc "$FOLDCHANGE" --parameter "$PARAMETER"
    cd "$BASE_DIR"
    
    # Extract significant reads
    mkdir -p "$RTD_SIG_PATH"
    cd "$RTD_SIG_PATH"
    python "$PIPELINE_DIR/extract_1.py" --clusters "$RTD_CLUST_PATH/clustering_results.tsv" \
        --sleuth "$RTD_SLEUTH_PATH/plotting.tsv" --alpha "$ALPHA" --foldchange "$FOLDCHANGE"
    "$PIPELINE_DIR/fetch_loop_1.sh" --pipeline "$PIPELINE_DIR" --input "$INPUT_PATH" --threads "$THREADS"
    cd "$BASE_DIR"
}

# Function to run CGR analysis
run_cgr() {
    log_step "Starting CGR analysis"
    
    # Generate CGR embedding
    mkdir -p "$CGR_PATH"
    "$PIPELINE_DIR/generate_cgr_1.sh" -b "$BASE_DIR" -i "$INPUT_PATH" -p "$PIPELINE_DIR" \
        -m metadata.tsv -s "$STRAND" -t "$THREADS"
    
    # Run clustering
    mkdir -p "$CGR_CLUST_PATH"
    cd "$CGR_CLUST_PATH"
    python "$PIPELINE_DIR/cluster_1.py" --file_pattern "$CGR_PATH/*/*.tsv.gz" \
        --output_file clustering_results.tsv --min_cluster_size "$MIN_CLUSTER" \
        --min_samples "$MIN_SAMPLES" --num_cpus "$THREADS" --max_GMM "$MAX_GMM" --seed "$SEED"
    cd "$BASE_DIR"
    
    # Run bootstrapping
    mkdir -p "$CGR_BOOTSTRAP_PATH"
    cd "$CGR_BOOTSTRAP_PATH"
    python "$PIPELINE_DIR/bootstrap_2.py" --counts "$CGR_CLUST_PATH/table.otu" \
        --volumes "$CGR_CLUST_PATH/est_volumes.tsv" --n_bootstraps "$N_BOOTSTRAPS"
    cd "$BASE_DIR"
    
    # Run Sleuth
    mkdir -p "$CGR_SLEUTH_PATH"
    cd "$CGR_SLEUTH_PATH"
    Rscript "$PIPELINE_DIR/sleuth_1.R" --root "$BASE_DIR" --metadata metadata.tsv \
        --bootstrap "07_CGR_bootstrap" --pseudocount "$PSEUDOCOUNT" --alpha "$ALPHA" \
        --fc "$FOLDCHANGE" --parameter "$PARAMETER"
    cd "$BASE_DIR"
    
    # Extract significant reads
    mkdir -p "$CGR_SIG_PATH"
    cd "$CGR_SIG_PATH"
    python "$PIPELINE_DIR/extract_1.py" --clusters "$CGR_CLUST_PATH/clustering_results.tsv" \
        --sleuth "$CGR_SLEUTH_PATH/plotting.tsv" --alpha "$ALPHA" --foldchange "$FOLDCHANGE"
    "$PIPELINE_DIR/fetch_loop_1.sh" --pipeline "$PIPELINE_DIR" --input "$INPUT_PATH" --threads "$THREADS"
    cd "$BASE_DIR"
}

# Function to run merge analysis
run_merge() {
    log_step "Starting merge analysis"
    
    # Run merge
    mkdir -p "$MERGE_PATH"
    cd "$MERGE_PATH"
    python "$PIPELINE_DIR/merge_1.py" --clusterings "$BASE_DIR/*_clust/clustering_results.tsv" \
        --volumes "$BASE_DIR/*_clust/est_volumes.tsv" --overlap "$OVERLAP"
    cd "$BASE_DIR"
    
    # Run bootstrapping
    mkdir -p "$MERGE_BOOTSTRAP_PATH"
    cd "$MERGE_BOOTSTRAP_PATH"
    python "$PIPELINE_DIR/bootstrap_2.py" --counts "$MERGE_PATH/merged_table.otu" \
        --volumes "$MERGE_PATH/merged_est_volumes.tsv" --n_bootstraps "$N_BOOTSTRAPS"
    cd "$BASE_DIR"
    
    # Run Sleuth
    mkdir -p "$MERGE_SLEUTH_PATH"
    cd "$MERGE_SLEUTH_PATH"
    Rscript "$PIPELINE_DIR/sleuth_1.R" --root "$BASE_DIR" --metadata metadata.tsv \
        --bootstrap "08_merge_bootstrap" --pseudocount "$PSEUDOCOUNT" --alpha "$ALPHA" \
        --fc "$FOLDCHANGE" --parameter "$PARAMETER"
    cd "$BASE_DIR"
    
    # Extract significant reads
    mkdir -p "$MERGE_SIG_PATH"
    cd "$MERGE_SIG_PATH"
    python "$PIPELINE_DIR/extract_1.py" --clusters "$MERGE_PATH/merged_clustering_results.tsv" \
        --sleuth "$MERGE_SLEUTH_PATH/plotting.tsv" --alpha "$ALPHA" --foldchange "$FOLDCHANGE"
    "$PIPELINE_DIR/fetch_loop_1.sh" --pipeline "$PIPELINE_DIR" --input "$INPUT_PATH" --threads "$THREADS"
    cd "$BASE_DIR"
}

# Main execution
log_step "Starting pipeline with embedding type: $EMBEDDING_TYPE"

# Step 1: Set up directory structure and filter/merge reads
mkdir -p "$INPUT_PATH"
log_step "Setting up input directory structure"

# Create initial directory structure based on metadata
while IFS=$'\t' read -r SRA description sample strain clone; do
    # Skip header
    [ "$SRA" = "SRA" ] && continue
    
    # Create SRA directory
    mkdir -p "$INPUT_PATH/$SRA"
    
    # Check if fastq files exist
    if [[ ! -f "$INPUT_PATH/$SRA/${SRA}_1.fastq.gz" || ! -f "$INPUT_PATH/$SRA/${SRA}_2.fastq.gz" ]]; then
        echo "Error: Required fastq files not found for $SRA"
        echo "Please ensure ${SRA}_1.fastq.gz and ${SRA}_2.fastq.gz exist in $INPUT_PATH/$SRA"
        exit 1
    fi
done < "metadata.tsv"

log_step "Running filter/merge step"
"$PIPELINE_DIR/filter_merge_1.sh" -p "$INPUT_PATH" -m "$BASE_DIR/metadata.tsv" -t "$THREADS"

# Run appropriate analysis based on embedding type
case "$EMBEDDING_TYPE" in
    "RTD")
        run_rtd
        ;;
    "CGR")
        run_cgr
        ;;
    "merge")
        run_rtd
        run_cgr
        run_merge
        ;;
esac

log_step "Pipeline completed successfully"

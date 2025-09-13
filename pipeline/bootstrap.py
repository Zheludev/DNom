import os
import h5py
import numpy as np
import pandas as pd
from sklearn.utils import resample
from tqdm import tqdm
import argparse
import time

def calculate_tpm(counts, lengths):
    """Convert raw read counts to TPM."""
    rpk = counts / lengths  # Reads Per Kilobase
    tpm = (rpk / np.sum(rpk)) * 1e6  # Transcripts Per Million
    return tpm

def calculate_eff_lengths(lengths):
    """Calculate effective lengths for Sleuth. For now, assume it's the same as gene lengths."""
    return lengths

def write_abundance_tsv(cluster_ids, read_counts, tpm_values, eff_lengths, sample_name):
    """
    Write a tab-delimited abundance.tsv file in Kallisto format.

    Parameters:
    -----------
    cluster_ids : array-like
        The target/gene identifiers
    read_counts : array-like
        The estimated counts for each target
    tpm_values : array-like
        The TPM values for each target
    eff_lengths : array-like
        The effective lengths for each target
    sample_name : str
        Name of the sample directory where the file will be written
    """
    # Create a DataFrame with the required columns
    abundance_df = pd.DataFrame({
        'target_id': cluster_ids,
        'length': eff_lengths,  # Including original length for reference
        'eff_length': eff_lengths,
        'est_counts': read_counts,
        'tpm': tpm_values
    })

    # Ensure the columns are in the correct order as per Kallisto format
    ordered_columns = ['target_id', 'length', 'eff_length', 'est_counts', 'tpm']
    abundance_df = abundance_df[ordered_columns]

    # Create the output directory if it doesn't exist
    if not os.path.exists(sample_name):
        os.makedirs(sample_name)

    # Write the file with tab separation and no index
    output_path = os.path.join(sample_name, 'abundance.tsv')
    abundance_df.to_csv(output_path, sep='\t', index=False)
    print(f"Abundance TSV file created at: {output_path}")

def generate_hdf5(cluster_ids, read_counts, gene_lengths, sample_name, n_bootstraps):
    """Generate an HDF5 file for a given sample with bootstrapped counts, TPM, lengths, and effective lengths."""
    # Adjust counts to TPM using gene lengths
    tpm_values = calculate_tpm(read_counts, gene_lengths)

    # Calculate effective lengths (for now, assume it's the same as gene lengths)
    eff_lengths = calculate_eff_lengths(gene_lengths)

    # Write the abundance.tsv file first
    write_abundance_tsv(cluster_ids, read_counts, tpm_values, eff_lengths, sample_name)

    # Perform bootstrapping for this sample with tqdm progress bar
    bootstrap_samples = []
    for _ in tqdm(range(n_bootstraps), desc=f"Bootstrapping {sample_name}", unit="iter"):
        bootstrap_sample = resample(read_counts, replace=True, n_samples=len(read_counts))
        bootstrap_samples.append(bootstrap_sample)
    bootstrap_samples = np.array(bootstrap_samples)  # Shape: (n_bootstraps, n_genes)

    # Capture the start time at the beginning of the process
    start_time = time.time()

    # Define the HDF5 file paths
    hdf5_file_path = os.path.join(sample_name, 'abundance.h5')

    # Create and write the HDF5 file
    with h5py.File(hdf5_file_path, 'w') as f:
        # Store cluster IDs (as bytes)
        f.create_dataset('/aux/ids', data=np.array(cluster_ids, dtype='S'))

        # Store read counts (these are also the "estimated counts")
        f.create_dataset('/est_counts', data=read_counts)

        # Store TPM values
        f.create_dataset('/tpm', data=tpm_values)

        # Store effective lengths (required by Sleuth)
        f.create_dataset('/aux/eff_lengths', data=eff_lengths)

        # Store raw lengths (required by Sleuth)
        f.create_dataset('/aux/lengths', data=gene_lengths)

        # Store the index version (required by Sleuth)
        f.create_dataset('/aux/index_version', data=np.bytes_('v0.50.1'))

        # Store the kallisto version (required by Sleuth)
        f.create_dataset('/aux/kallisto_version', data=np.bytes_('v0.50.1'))

        # Store the system start time (in seconds since the epoch)
        f.create_dataset('/aux/start_time', data=start_time)

        # Store the number of bootstrap samples
        f.create_dataset('/aux/num_bootstrap', data=n_bootstraps)

        # Store bootstrap replicates
        bootstrap_group = f.create_group('/bootstrap')
        for i in range(n_bootstraps):
            bootstrap_group.create_dataset(f'bs{i}', data=bootstrap_samples[i])

    print(f"HDF5 file for {sample_name} created at: {hdf5_file_path}")

def main(counts_file, volumes_file, n_bootstraps):
    """Main function to process the input files and generate HDF5 files."""
    # Load the counts file (tab-delimited)
    counts_data = pd.read_csv(counts_file, sep='\t')

    # Load the volumes file (tab-delimited)
    volumes_data = pd.read_csv(volumes_file, sep='\t')
    
    # Extract the cluster IDs from both files
    counts_cluster_ids = counts_data.iloc[:, 0].values
    volumes_cluster_ids = volumes_data['clusterID'].values
    
    # Debug: Print sample cluster IDs to diagnose the issue
    print(f"First 5 counts cluster IDs: {counts_cluster_ids[:5]}")
    print(f"First 5 volumes cluster IDs: {volumes_cluster_ids[:5]}")
    
    # Function to normalize cluster IDs to facilitate matching
    def normalize_cluster_id(cluster_id):
        cluster_id = str(cluster_id).strip()
        # Remove type prefixes if present (e.g., "RTD_" or "CGR_")
        if "_" in cluster_id and cluster_id.split("_")[0] in ["RTD", "CGR"]:
            cluster_id = "_".join(cluster_id.split("_")[1:])
        return cluster_id
    
    # Normalize the IDs before comparison
    normalized_counts_ids = [normalize_cluster_id(cid) for cid in counts_cluster_ids]
    normalized_volumes_ids = [normalize_cluster_id(cid) for cid in volumes_cluster_ids]
    
    # Create mapping dictionaries to original IDs
    counts_id_map = dict(zip(normalized_counts_ids, counts_cluster_ids))
    volumes_id_map = dict(zip(normalized_volumes_ids, volumes_cluster_ids))
    
    # Check if normalized sets of IDs have overlap
    common_normalized_ids = set(normalized_counts_ids) & set(normalized_volumes_ids)
    print(f"After normalization: Found {len(common_normalized_ids)} common clusters")
    
    if len(common_normalized_ids) == 0:
        # Try stripping all prefixes as a last resort
        def strip_all_prefixes(cluster_id):
            # Keep only the numeric part of the ID
            return ''.join(c for c in str(cluster_id) if c.isdigit())
        
        normalized_counts_ids = [strip_all_prefixes(cid) for cid in counts_cluster_ids]
        normalized_volumes_ids = [strip_all_prefixes(cid) for cid in volumes_cluster_ids]
        
        # Create new mapping dictionaries
        counts_id_map = dict(zip(normalized_counts_ids, counts_cluster_ids))
        volumes_id_map = dict(zip(normalized_volumes_ids, volumes_cluster_ids))
        
        # Check if we have a match now
        common_normalized_ids = set(normalized_counts_ids) & set(normalized_volumes_ids)
        print(f"After stripping all prefixes: Found {len(common_normalized_ids)} common clusters")
    
    # If we still have no overlap, raise a more informative error
    if len(common_normalized_ids) == 0:
        raise ValueError("No common clusters found between counts and volumes files after normalization attempts. Check that both files refer to the same clusters.")
    
    # Map the normalized IDs back to their original form in each file
    common_counts_ids = [counts_id_map[nid] for nid in common_normalized_ids if nid in counts_id_map]
    common_volumes_ids = [volumes_id_map[nid] for nid in common_normalized_ids if nid in volumes_id_map]
    
    # Now filter the data to include only these matched clusters
    counts_data = counts_data[counts_data.iloc[:, 0].isin(common_counts_ids)]
    volumes_data = volumes_data[volumes_data['clusterID'].isin(common_volumes_ids)]
    
    # Create mappings between the IDs to ensure consistent ordering
    id_mappings = dict(zip(common_counts_ids, common_volumes_ids))
    
    # Adjust the cluster IDs in the counts file to match those in volumes
    counts_data.iloc[:, 0] = counts_data.iloc[:, 0].map(lambda x: id_mappings.get(x, x))
    
    # Extract the gene IDs
    cluster_ids = counts_data.iloc[:, 0].values
    
    # Handle zero values in estVol by replacing them with the smallest non-zero value
    gene_lengths = volumes_data['estVol'].values
    min_nonzero_length = np.min(gene_lengths[gene_lengths > 0])
    
    # Log information about zero values being replaced
    zero_count = np.sum(gene_lengths == 0)
    if zero_count > 0:
        print(f"Found {zero_count} genes with zero volume. Replacing with minimum non-zero volume: {min_nonzero_length}")
        gene_lengths[gene_lengths == 0] = min_nonzero_length

    # Iterate over each sample (columns after the first one)
    for sample_col in counts_data.columns[1:]:
        # Extract the sample name (everything before the first '.')
        sample_name = sample_col.split('.')[0]

        # Extract counts for this sample
        read_counts = counts_data[sample_col].values

        # Generate HDF5 file for this sample
        generate_hdf5(cluster_ids, read_counts, gene_lengths, sample_name, n_bootstraps)

if __name__ == "__main__":
    # Argument parsing
    parser = argparse.ArgumentParser(description='Generate kallisto-like HDF5 files with bootstrapping.')
    parser.add_argument('--counts', type=str, required=True, help='Path to the input tab-delimited gene count file')
    parser.add_argument('--volumes', type=str, required=True, help='Path to the input gene volumes file (for length estimates)')
    parser.add_argument('--n_bootstraps', type=int, default=100, help='Number of bootstrap replicates (default: 100)')

    args = parser.parse_args()

    # Call main function with parsed arguments
    main(args.counts, args.volumes, args.n_bootstraps)


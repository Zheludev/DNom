## this code takes as input, groups of two files:
    ## 1. 'clustering_results.tsv'
        ## table of all reads, file origin, and clusterID (including noise)
            ## tab-delimited
    ## 2. 'est_volumes.tsv'
        ## clusterID, estimated cluster volumes, and GMM components
            ## tab-delimited

## the code then also takes the parameter 'overlap', which specifies a (read) count

## over >1 groups of files (corresponding to individual HDBSCAN clusterings),
    ## the code then identifies clusters between groups that have shared reads
        ## if the shared reads overlap is >= 'overlap', then a new merged cluster is created
            ## this iterates until all overlaps are accounted for
        ## only cluster-cluster overlaps are considered (cluster-noise overlaps are ignored)

## the resulting merged clusters are then written to a new merged 'clustering_results.tsv'
    ## their estimated volumes / components are summed -> new 'est_volumes.tsv'
    ## and from 'clustering_results.tsv', a tab-delimted 'table.otu' is written

import gzip
import pandas as pd
import numpy as np
from tqdm import tqdm
import argparse
import os
import glob

def sort_cluster_ids(cluster_ids):
    """Sort cluster IDs numerically, ignoring any text prefixes."""
    def extract_number(cluster_id):
        # Extract numerical part from cluster ID
        num = ''.join(filter(str.isdigit, str(cluster_id)))
        return int(num) if num else float('inf')
    
    return sorted(cluster_ids, key=extract_number)

def extract_type(path):
    """Extract clustering type from directory path."""
    # Get the directory name containing clustering_results.tsv
    dirname = os.path.basename(os.path.dirname(path))
    return dirname

def load_data(clustering_patterns, volume_patterns):
    """Load and combine clustering and volume data from multiple files."""
    clustering_dfs = []
    volume_dfs = []
    clustering_paths = []
    volume_paths = []

    # Handle the clustering patterns
    if isinstance(clustering_patterns, list):
        for pattern in clustering_patterns:
            clustering_paths.extend(sorted(glob.glob(pattern)))
    else:
        clustering_paths = sorted(glob.glob(clustering_patterns))

    # Handle the volume patterns
    if isinstance(volume_patterns, list):
        for pattern in volume_patterns:
            volume_paths.extend(sorted(glob.glob(pattern)))
    else:
        volume_paths = sorted(glob.glob(volume_patterns))

    # Verify we found files
    if not clustering_paths:
        raise FileNotFoundError(f"No clustering files found matching pattern: {clustering_patterns}")
    if not volume_paths:
        raise FileNotFoundError(f"No volume files found matching pattern: {volume_patterns}")

    # Store actual paths for later use
    load_data.actual_paths = clustering_paths

    # Verify equal number of files
    if len(clustering_paths) != len(volume_paths):
        raise ValueError(f"Found {len(clustering_paths)} clustering files but {len(volume_paths)} volume files")

    print(f"Found {len(clustering_paths)} clustering files:")
    for path in clustering_paths:
        print(f"  {path}")
    print(f"\nFound {len(volume_paths)} volume files:")
    for path in volume_paths:
        print(f"  {path}")
    print()

    for clust_path, vol_path in zip(clustering_paths, volume_paths):
        # Extract clustering type from directory
        clust_type = extract_type(clust_path)

        # Check if files exist (redundant but safe)
        if not os.path.exists(clust_path):
            raise FileNotFoundError(f"Clustering file not found: {clust_path}")
        if not os.path.exists(vol_path):
            raise FileNotFoundError(f"Volume file not found: {vol_path}")

        # Load clustering data
        print(f"Loading clustering data from: {clust_path}")
        clust_df = pd.read_csv(clust_path, sep='\t')
        # Add clustering type as a column
        clust_df['type'] = clust_type
        clustering_dfs.append(clust_df)

        # Load volume data
        print(f"Loading volume data from: {vol_path}")
        vol_df = pd.read_csv(vol_path, sep='\t')
        # Remove the column renaming
        volume_dfs.append(vol_df)

    print("Clustering DFs:", [df.shape for df in clustering_dfs])
    print("Volume DFs:", [df.shape for df in volume_dfs])

    # Combine all clustering data
    combined_clustering = pd.concat(clustering_dfs, ignore_index=True)
    combined_volumes = pd.concat(volume_dfs, ignore_index=True)

    return combined_clustering, combined_volumes

def find_overlapping_clusters(clustering_df, min_overlap, clustering_paths):
    """Find clusters that share reads above the minimum overlap threshold."""
    print("\nFinding overlapping clusters...")
    
    # Get actual paths if they were stored during load_data
    actual_paths = getattr(load_data, 'actual_paths', None)
    if actual_paths:
        clustering_paths = actual_paths
    
    # Get clustering types from the actual paths
    clustering_types = sorted({extract_type(path) for path in clustering_paths})
    print(f"Identified clustering types: {clustering_types}")
    
    # Verify clustering types in DataFrame
    print("\nClustering types found in data:")
    value_counts = clustering_df['type'].value_counts()
    for ctype, count in value_counts.items():
        print(f"  {ctype}: {count:,} reads")
    
    # Verify we have the expected number of types
    if len(clustering_types) < 2:
        print("\nError: Not enough clustering types found!")
        print(f"Clustering paths: {clustering_paths}")
        print(f"Found types: {clustering_types}")
        print(f"Types in data: {sorted(clustering_df['type'].unique())}")
        raise ValueError(f"Expected at least 2 clustering types, found: {clustering_types}")
    
    # Create mappings
    read_to_cluster = {}
    cluster_to_type = {}
    cluster_reads = {}
    
    # Only process non-noise clusters and store initial cluster counts
    first_type = clustering_types[0]
    second_type = clustering_types[1]
    
    # Get initial cluster counts before filtering
    initial_counts = {
        ctype: len(clustering_df[(clustering_df['type'] == ctype) & 
                               (clustering_df['clusterID'] != 'noise')]['clusterID'].unique())
        for ctype in clustering_types
    }
    
    # Store initial counts for summary
    find_overlapping_clusters.initial_counts = initial_counts
    
    # Only process non-noise clusters
    non_noise_df = clustering_df[clustering_df['clusterID'] != 'noise']
    
    # Build mappings
    print("\nBuilding read and cluster mappings...")
    for _, row in tqdm(non_noise_df.iterrows(), total=len(non_noise_df)):
        read_id = row['readID']
        cluster_id = row['clusterID']
        ctype = row['type']
        
        if read_id not in read_to_cluster:
            read_to_cluster[read_id] = set()
        read_to_cluster[read_id].add(cluster_id)
        
        cluster_to_type[cluster_id] = ctype
        if cluster_id not in cluster_reads:
            cluster_reads[cluster_id] = set()
        cluster_reads[cluster_id].add(read_id)
    
    # Find clusters for each type
    first_clusters = {cluster_id for cluster_id, ctype in cluster_to_type.items() 
                     if ctype == first_type}
    second_clusters = {cluster_id for cluster_id, ctype in cluster_to_type.items() 
                      if ctype == second_type}
    
    print(f"\nComparing clusters between {first_type} ({initial_counts[first_type]} clusters) and "
          f"{second_type} ({initial_counts[second_type]} clusters)")
    
    # Find overlapping clusters
    overlapping_pairs = []
    
    # For each cluster in first type
    for cluster1 in tqdm(first_clusters):
        reads1 = cluster_reads[cluster1]
        # Find all reads from this cluster that appear in any cluster from second type
        for cluster2 in second_clusters:
            reads2 = cluster_reads[cluster2]
            overlap = len(reads1 & reads2)
            if overlap >= min_overlap:
                overlapping_pairs.append((cluster1, cluster2))
    
    print(f"\nFound {len(overlapping_pairs)} overlapping cluster pairs")
    return overlapping_pairs

def merge_clusters(clustering_df, volume_df, overlapping_pairs):
    """Merge overlapping clusters and update their volumes."""
    print("\nStarting cluster merging process...")

    if not overlapping_pairs:
        print("No overlapping clusters found. Renaming existing clusters...")
        non_noise_clusters = clustering_df[clustering_df['clusterID'] != 'noise']['clusterID'].unique()
        rename_dict = {cluster: f's_{i}' for i, cluster in enumerate(sorted(non_noise_clusters))}

        clustering_df.loc[clustering_df['clusterID'] != 'noise', 'clusterID'] = (
            clustering_df.loc[clustering_df['clusterID'] != 'noise', 'clusterID'].map(rename_dict)
        )
        volume_df['clusterID'] = volume_df['clusterID'].map(rename_dict)

        print(f"Renamed {len(rename_dict)} clusters.")
        return clustering_df, volume_df

    print(f"Found {len(overlapping_pairs)} overlapping pairs. Building cluster groups...")

    # Create connected components
    from collections import defaultdict
    graph = defaultdict(set)
    for c1, c2 in overlapping_pairs:
        graph[c1].add(c2)
        graph[c2].add(c1)

    # Find connected components using DFS
    def dfs(node, component, visited):
        visited.add(node)
        component.add(node)
        for neighbor in graph[node]:
            if neighbor not in visited:
                dfs(neighbor, component, visited)

    visited = set()
    merged_groups = []
    print("Finding connected components...")
    for node in tqdm(graph):
        if node not in visited:
            component = set()
            dfs(node, component, visited)
            merged_groups.append(component)

    merged_clusters = set().union(*merged_groups) if merged_groups else set()

    print("Identifying singleton clusters...")
    all_non_noise_clusters = set(clustering_df[clustering_df['clusterID'] != 'noise']['clusterID'].unique())
    singleton_clusters = sorted(all_non_noise_clusters - merged_clusters)

    # Prepare renaming operations
    print("Preparing cluster renaming...")
    rename_dict = {}

    # Handle merged clusters
    for i, group in enumerate(merged_groups):
        new_cluster_id = f'm_{i}'
        for old_id in group:
            rename_dict[old_id] = new_cluster_id

    # Handle singleton clusters
    for i, cluster in enumerate(singleton_clusters):
        rename_dict[cluster] = f's_{i}'

    print("Applying cluster renaming...")
    clustering_df.loc[clustering_df['clusterID'] != 'noise', 'clusterID'] = (
        clustering_df.loc[clustering_df['clusterID'] != 'noise', 'clusterID'].map(rename_dict)
    )

    print("Updating volume information...")
    new_volumes = []

    # Process merged clusters
    print("Processing merged clusters...")
    for i, group in tqdm(enumerate(merged_groups), total=len(merged_groups)):
        merged_volume = volume_df[volume_df['clusterID'].isin(group)]['estVol'].sum()
        merged_components = volume_df[volume_df['clusterID'].isin(group)]['components'].sum()
        new_volumes.append({
            'clusterID': f'm_{i}',
            'estVol': merged_volume,
            'components': merged_components
        })

    # Process singleton clusters
    print("Processing singleton clusters...")
    for i, cluster in tqdm(enumerate(singleton_clusters), total=len(singleton_clusters)):
        vol_row = volume_df[volume_df['clusterID'] == cluster].iloc[0]
        new_volumes.append({
            'clusterID': f's_{i}',
            'estVol': vol_row['estVol'],
            'components': vol_row['components']
        })

    # Create new volume DataFrame
    volume_df = pd.DataFrame(new_volumes)
    
    # Sort using our consistent sorting function
    sorted_cluster_ids = sort_cluster_ids(volume_df['clusterID'].unique())
    volume_df['clusterID'] = pd.Categorical(
        volume_df['clusterID'],
        categories=sorted_cluster_ids,
        ordered=True
    )
    volume_df = volume_df.sort_values('clusterID').reset_index(drop=True)
    
    return clustering_df, volume_df

def print_summary(clustering_paths, clustering_df, volume_df, overlap, overlapping_pairs):
    """Print summary statistics of the clustering process."""
    print("\nMERGING SUMMARY")
    print("-----------------")

    # Total reads
    total_reads = len(clustering_df['readID'].unique())
    print(f"\nTotal reads processed: {total_reads:,}")

    # Get clustering types from input paths
    actual_paths = getattr(load_data, 'actual_paths', None)
    if actual_paths:
        clustering_paths = actual_paths
    clustering_types = sorted({extract_type(path) for path in clustering_paths})

    # Get initial cluster counts stored at the start
    initial_counts = getattr(main, 'initial_counts', {})

    # Reads per clustering type
    print("\nReads per clustering type:")
    for ctype in clustering_types:
        reads = len(clustering_df[clustering_df['type'] == ctype]['readID'].unique())
        print(f"  {ctype}: {reads:,}")

    # Overlap threshold
    print(f"\nOverlap threshold used: {overlap}")

    # Initial clusters per type
    print("\nInitial clusters per clustering type:")
    for ctype in clustering_types:
        print(f"  {ctype}: {initial_counts.get(ctype, 0)}")

    # Number of overlapping clusters
    unique_clusters = set()
    for c1, c2 in overlapping_pairs:
        unique_clusters.add(c1)
        unique_clusters.add(c2)
    print(f"\nNumber of overlapping clusters identified: {len(unique_clusters)}")

    # Final cluster counts
    merged_clusters = len(volume_df[volume_df['clusterID'].str.startswith('m_')])
    singleton_clusters = len(volume_df[volume_df['clusterID'].str.startswith('s_')])
    print(f"\nFinal clusters generated:")
    print(f"  Merged clusters (m_X): {merged_clusters}")
    print(f"  Singleton clusters (s_X): {singleton_clusters}")
    print(f"  Total non-noise clusters: {merged_clusters + singleton_clusters}")

    # Volume per clustering type and total volumes
    merged_volume = volume_df[volume_df['clusterID'].str.startswith('m_')]['estVol'].sum()
    singleton_volume = volume_df[volume_df['clusterID'].str.startswith('s_')]['estVol'].sum()
    total_volume = merged_volume + singleton_volume
    print(f"\nFinal estimated volumes:")
    print(f"  Merged clusters (m_X): {merged_volume:,.2e}")
    print(f"  Singleton clusters (s_X): {singleton_volume:,.2e}")
    print(f"  Total: {total_volume:,.2e}")

def sort_cluster_ids(cluster_ids):
    """Sort cluster IDs in a consistent order: all 'm_' clusters first, then all 's_' clusters."""
    def extract_info(cluster_id):
        # Split into prefix (m_ or s_) and number
        prefix = cluster_id[0]  # 'm' or 's'
        number = int(''.join(filter(str.isdigit, cluster_id)))
        # Make 'm' come before 's' by using 0 for 'm' and 1 for 's'
        prefix_value = 0 if prefix == 'm' else 1
        return (prefix_value, number)
    
    return sorted(cluster_ids, key=extract_info)

def save_otu(read_ids, filenames, cluster_labels, output_file="merged_table.otu"):
    """Create and save OTU table from clustering results."""
    result_df = pd.DataFrame({
        'readID': read_ids,
        'filename': filenames,
        'clusterID': cluster_labels
    })
    
    # Process filename
    result_df['filename'] = result_df['filename'].str.split('.').str[0]
    
    # Filter out noise and create the grouped table
    non_noise_df = result_df[result_df['clusterID'] != 'noise']
    grouped = non_noise_df.groupby(['filename', 'clusterID']).size().unstack(fill_value=0)
    
    # Sort columns using our consistent sorting function
    sorted_cluster_ids = sort_cluster_ids(grouped.columns.tolist())
    grouped = grouped[sorted_cluster_ids]
    
    # Transpose and save
    transposed_grouped = grouped.T
    transposed_grouped.to_csv(output_file, sep='\t', index_label='', header=True)
    print(f"OTU table saved as '{output_file}'")

def main(clustering_paths, volume_paths, overlap, output_clustering, output_volume, output_otu):
    """Main function to merge overlapping clusters across multiple clusterings."""
    # Load all data
    clustering_df, volume_df = load_data(clustering_paths, volume_paths)
    
    # Get clustering types
    actual_paths = getattr(load_data, 'actual_paths', None)
    if actual_paths:
        clustering_paths = actual_paths
    clustering_types = sorted({extract_type(path) for path in clustering_paths})
    
    # Store initial cluster counts before any merging
    initial_counts = {}
    for ctype in clustering_types:
        count = len(clustering_df[(clustering_df['type'] == ctype) & 
                                (clustering_df['clusterID'] != 'noise')]['clusterID'].unique())
        initial_counts[ctype] = count
    
    # Store for use in summary
    main.initial_counts = initial_counts
    
    # Track all overlapping pairs for summary
    all_overlapping_pairs = []
    
    while True:
        # Find overlapping clusters
        overlapping_pairs = find_overlapping_clusters(clustering_df, overlap, clustering_paths)
        
        if not overlapping_pairs:
            break
        
        all_overlapping_pairs.extend(overlapping_pairs)
        
        # Merge overlapping clusters
        clustering_df, volume_df = merge_clusters(clustering_df, volume_df, overlapping_pairs)
    
    # Print summary
    print_summary(clustering_paths, clustering_df, volume_df, overlap, all_overlapping_pairs)
    
    # Remove type column before saving
    if 'type' in clustering_df.columns:
        clustering_df = clustering_df.drop('type', axis=1)
    
    # Ensure consistent sorting of cluster IDs
    sorted_cluster_ids = sort_cluster_ids(clustering_df[clustering_df['clusterID'] != 'noise']['clusterID'].unique())
    
    # Sort volume_df
    volume_df['clusterID'] = pd.Categorical(
        volume_df['clusterID'],
        categories=sorted_cluster_ids,
        ordered=True
    )
    volume_df = volume_df.sort_values('clusterID').reset_index(drop=True)
    
    # Save results
    clustering_df.to_csv(output_clustering, sep='\t', index=False)
    volume_df.to_csv(output_volume, sep='\t', index=False)
    
    # Create and save OTU table (already uses sort_cluster_ids internally)
    save_otu(
        clustering_df['readID'].values,
        clustering_df['filename'].values,
        clustering_df['clusterID'].values,
        output_otu
    )
    
    return clustering_df, volume_df

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Merge overlapping clusters from multiple clusterings.')
    parser.add_argument('--clusterings', type=str, required=True, help='Glob pattern for clustering result files (e.g. "*_clust/clustering_results.tsv")')
    parser.add_argument('--volumes', type=str, required=True, help='Glob pattern for volume estimation files (e.g. "*_clust/est_volumes.tsv")')
    parser.add_argument('--overlap', type=int, default=10, help='Minimum overlap threshold')
    parser.add_argument('--out-clustering', default='merged_clustering_results.tsv',
                      help='Output path for merged clustering results')
    parser.add_argument('--out-volume', default='merged_est_volumes.tsv',
                      help='Output path for merged volume estimations')
    parser.add_argument('--out-otu', default='merged_table.otu',
                      help='Output path for OTU table')
    
    args = parser.parse_args()
    
    main(args.clusterings, args.volumes, args.overlap,
         args.out_clustering, args.out_volume, args.out_otu)

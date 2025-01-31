import pandas as pd
import argparse
import os
import math

def create_directories():
    """
    Creates the directory structure for organizing results by test type and direction.
    Creates: Wald/positive, Wald/negative, LRT/positive, LRT/negative
    """
    for test_type in ['Wald', 'LRT']:
        for direction in ['positive', 'negative']:
            os.makedirs(os.path.join(test_type, direction), exist_ok=True)

def filter_significant_clusters(sleuth_file, clusters_file, alpha_threshold, foldchange, test_type):
    """
    Filters clusters based on significance criteria for either Wald or LRT tests.
    
    Parameters:
    sleuth_file (str): Path to the plotting.tsv file
    clusters_file (str): Path to clustering results file
    alpha_threshold (float): Alpha threshold for significance
    foldchange (float): Fold change threshold (will be converted to log2)
    test_type (str): Either 'Wald' or 'LRT' to determine which test results to use
    """
    # Load the Sleuth file
    sleuth_df = pd.read_csv(sleuth_file, sep='\t')
    
    # Calculate log2 fold change threshold
    log2_fc_threshold = math.log2(foldchange)
    
    # Get total number of clusters
    total_clusters = sleuth_df.shape[0]
    
    # Process both positive and negative directions
    for direction in ['positive', 'negative']:
        # Set up the output directory path
        output_dir = os.path.join(test_type, direction)
        
        # Define filter conditions based on test type and direction
        if test_type == 'Wald':
            alpha_column = 'wt_qval'
        else:  # LRT
            alpha_column = 'lrt_qval'
            
        # First filter based on alpha threshold
        significant_df = sleuth_df[sleuth_df[alpha_column] <= alpha_threshold].copy()
        
        # Then filter based on fold change direction
        if direction == 'positive':
            significant_df = significant_df[significant_df['mean_log2fc'] >= log2_fc_threshold]
        else:  # negative
            significant_df = significant_df[significant_df['mean_log2fc'] <= -log2_fc_threshold]
        
        # Get number of significant clusters
        num_significant_clusters = len(significant_df)
        
        # Calculate percentage of significant clusters
        percentage_significant = (num_significant_clusters / total_clusters) * 100
        
        # Write summary statistics to a file in the appropriate directory
        with open(os.path.join(output_dir, 'summary_statistics.txt'), 'w') as f:
            f.write(f"Total clusters analyzed: {total_clusters}\n")
            f.write(f"Significant clusters ({test_type} test, {direction} regulation): {num_significant_clusters}\n")
            f.write(f"Percentage significant: {percentage_significant:.2f}%\n")
            f.write(f"Alpha threshold: {alpha_threshold}\n")
            f.write(f"Fold change threshold: {foldchange}x (log2: {log2_fc_threshold:.2f})\n")
        
        if num_significant_clusters == 0:
            # Create empty file indicating no significant clusters
            with open(os.path.join(output_dir, 'NO_SIGNIFICANT_CLUSTERS_FOUND'), 'w') as f:
                pass
            continue
        
        # Extract the list of significant cluster IDs
        significant_clusters = significant_df['target_id'].tolist()
        
        # Load and filter the clustering results
        clusters_df = pd.read_csv(clusters_file, sep='\t')
        filtered_clusters_df = clusters_df[clusters_df['clusterID'].isin(significant_clusters)].copy()
        
        # Extract numeric part of clusterID and cast to int for sorting
        filtered_clusters_df['numeric_clusterID'] = filtered_clusters_df['clusterID'].str.extract(r'(\d+)').astype(int)
        filtered_clusters_df['numeric_readID'] = filtered_clusters_df['readID'].str.extract(r'\.(\d+)$').astype(int)
        
        # Sort by clusterID (numerically), filename (alphabetically), then readID (numerically)
        filtered_clusters_df = filtered_clusters_df.sort_values(
            by=['numeric_clusterID', 'filename', 'numeric_readID']
        )
        
        # Drop the temporary numeric columns
        filtered_clusters_df = filtered_clusters_df.drop(columns=['numeric_clusterID', 'numeric_readID'])
        
        # Save to output file in the appropriate directory
        output_file = os.path.join(output_dir, 'significant.tsv')
        filtered_clusters_df.to_csv(output_file, sep='\t', index=False)

def main(sleuth_file, clusters_file, alpha_threshold, foldchange):
    # Create directory structure
    create_directories()
    
    # Process Wald test results
    filter_significant_clusters(sleuth_file, clusters_file, alpha_threshold, foldchange, 'Wald')
    
    # Process LRT test results
    filter_significant_clusters(sleuth_file, clusters_file, alpha_threshold, foldchange, 'LRT')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Filter and sort clustering results based on statistical significance and fold change.')
    parser.add_argument('--sleuth', default='plotting.tsv', help='Path to plotting.tsv file (default: plotting.tsv)')
    parser.add_argument('--clusters', required=True, help='Path to clustering_results.tsv file')
    parser.add_argument('--alpha', type=float, default=0.05, help='Alpha threshold for filtering (default: 0.05)')
    parser.add_argument('--foldchange', type=float, default=10.0, help='Fold change threshold (default: 10.0)')
    
    args = parser.parse_args()
    
    # Call main function with parsed arguments
    main(args.sleuth, args.clusters, args.alpha, args.foldchange)

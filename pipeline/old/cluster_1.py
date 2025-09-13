## draft clustering algorithm based on HDBSCAN co-written with Claude/ChatGPT4o

import warnings
warnings.filterwarnings('ignore', message='.*force_all_finite.*', category=FutureWarning)
from sklearn.utils.deprecation import _is_deprecated
warnings.filterwarnings('ignore', category=FutureWarning, module='sklearn.utils.deprecation')
warnings.filterwarnings('ignore', category=FutureWarning, module='sklearn.utils.validation')
import gzip
import pandas as pd
import fast_hdbscan as hdbscan
import glob
import os
import argparse
import numba
import joblib
from joblib import Parallel, delayed
import time
from tqdm import tqdm
import numpy as np
from sklearn.mixture import GaussianMixture
from scipy.stats import chi2
import scipy.special
from concurrent.futures import ThreadPoolExecutor
import logging
import sys
from sklearn.mixture import GaussianMixture
import multiprocessing

# Initiate number of threads
def setup_parallel_env(num_cpus):
	if num_cpus == -1:
		num_cpus = multiprocessing.cpu_count()
	
	# Set Numba threads
	numba.set_num_threads(num_cpus)
	
	# Set OpenMP threads (used by some numerical libraries)
	os.environ['OMP_NUM_THREADS'] = str(num_cpus)
	
	# Set MKL threads (if using Intel MKL)
	try:
		import mkl
		mkl.set_num_threads(num_cpus)
	except ImportError:
		pass
	
	return num_cpus

# Time measurement decorator
def time_function(func):
	def timed_func(*args, **kwargs):
		start_time = time.time()
		result = func(*args, **kwargs)
		end_time = time.time()
		elapsed_time = end_time - start_time
		print(f"Function '{func.__name__}' took {elapsed_time:.4f} seconds to complete.")
		return result
	return timed_func

def setup_logging(logfile='cluster.log'):
	# Create a logger
	logger = logging.getLogger()
	logger.setLevel(logging.INFO)

	# Create a file handler to write log output to a file
	file_handler = logging.FileHandler(logfile)
	file_handler.setLevel(logging.INFO)

	# Create a stream handler to write log output to the terminal
	console_handler = logging.StreamHandler(sys.stdout)
	console_handler.setLevel(logging.INFO)

	# Set format for the log messages
	formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
	file_handler.setFormatter(formatter)
	console_handler.setFormatter(formatter)

	# Add handlers to the logger
	logger.addHandler(file_handler)
	logger.addHandler(console_handler)

	return logger

def load_data(file_pattern):
	files = glob.glob(file_pattern)
	data_frames = []
	for file in files:
		with gzip.open(file, 'rt') as f:
			df = pd.read_csv(f, sep='\t')
			df['filename'] = os.path.basename(file)
			data_frames.append(df)
	return pd.concat(data_frames, ignore_index=True)

@time_function
def cluster_data(data, min_cluster_size, min_samples, num_cpus):
	clusterer = hdbscan.HDBSCAN(
		min_cluster_size=min_cluster_size,
		min_samples=min_samples,
		metric='euclidean',
		core_dist_num_cpus=num_cpus,
		n_jobs=num_cpus
	)
	cluster_labels = clusterer.fit_predict(data)
	
	cluster_labels = [f'c_{label}' if label != -1 else 'noise' for label in cluster_labels]
	joblib.dump((clusterer, cluster_labels), 'HDBSCAN.pkl', compress=9)
	print("HDBSCAN model saved to 'HDBSCAN.pkl'")
	
	return cluster_labels, clusterer

@time_function
def precomputed(hdbscan_file):
	print(f"Loading precomputed HDBSCAN model and cluster labels from '{hdbscan_file}'")
	clusterer, cluster_labels = joblib.load(hdbscan_file)
	return cluster_labels, clusterer

def sort_cluster_ids(cluster_labels):
	def sort_key(cluster_id):
		if cluster_id == 'noise':
			return float('inf')
		return int(cluster_id.split('_')[1])
	return sorted(cluster_labels, key=sort_key)

def save_results(read_ids, filenames, cluster_labels, output_file):
	result_df = pd.DataFrame({
		'readID': read_ids,
		'filename': filenames,
		'clusterID': cluster_labels
	})
	result_df.to_csv(output_file, sep='\t', index=False)

def save_abundances(read_ids, filenames, cluster_labels):
	result_df = pd.DataFrame({
		'readID': read_ids,
		'filename': filenames,
		'clusterID': cluster_labels
	})
	non_noise_df = result_df[result_df['clusterID'] != 'noise']
	grouped = non_noise_df.groupby(['filename', 'clusterID']).size().reset_index(name='readCount')
	grouped['clusterID'] = pd.Categorical(grouped['clusterID'], sorted(set(sort_cluster_ids(grouped['clusterID']))), ordered=True)

	for filename in grouped['filename'].unique():
		file_abundances = grouped[grouped['filename'] == filename].drop(columns=['filename'])
		file_abundances = file_abundances[['clusterID', 'readCount']]
		output_abundance_file = f"{filename.split('.')[0]}_abundances.tsv"
		file_abundances.sort_values(by='clusterID').to_csv(output_abundance_file, sep='\t', index=False)
		print(f"Abundances saved to {output_abundance_file}")

def save_otu(read_ids, filenames, cluster_labels):
	result_df = pd.DataFrame({
		'readID': read_ids,
		'filename': filenames,
		'clusterID': cluster_labels
	})
	result_df['filename'] = result_df['filename'].str.split('.').str[0]
	non_noise_df = result_df[result_df['clusterID'] != 'noise']
	grouped = non_noise_df.groupby(['filename', 'clusterID']).size().unstack(fill_value=0)
	sorted_cluster_ids = sort_cluster_ids(grouped.columns.tolist())
	grouped = grouped[sorted_cluster_ids]
	transposed_grouped = grouped.T
	transposed_grouped.to_csv("table.otu", sep='\t', index_label='', header=True)
	print("OTU table saved as 'table.otu'")

def estimate_optimal_components(cluster_data, max_components, random_state):
	n_components_range = range(1, max_components + 1)
	bic_scores = []

	for n_components in n_components_range:
		gmm = GaussianMixture(n_components=n_components, random_state=random_state)
		gmm.fit(cluster_data)
		bic_scores.append(gmm.bic(cluster_data))

	optimal_components = n_components_range[np.argmin(bic_scores)]
	return optimal_components

# Function to estimate GMM-based volume for a cluster, factoring out duplicate points
def estimate_GMM_volume(cluster_data, max_components, random_state):
	# Remove duplicate points and keep track of their frequency
	unique_data, counts = np.unique(cluster_data, axis=0, return_counts=True)
	n_samples = len(unique_data)

	if n_samples < 2:
		# Too few points to estimate a volume, return 0
		return 0, 0

	# If the number of samples is greater than the max_components, find the optimal number of components
	if n_samples > max_components:
		n_components = estimate_optimal_components(unique_data, max_components, random_state)
	else:
		# If not enough points, set n_components to be the number of samples
		n_components = min(max_components, n_samples)

	# Fit a Gaussian Mixture Model to the cluster data
	gmm = GaussianMixture(n_components=n_components, random_state=random_state)
	gmm.fit(unique_data)

	# Calculate the total volume as the sum of volumes for each component
	total_volume = 0
	n_dims = unique_data.shape[1]
	confidence_interval = 0.95  # You can adjust this value
	chi2_val = chi2.ppf(confidence_interval, n_dims)

	for weight, covariance in zip(gmm.weights_, gmm.covariances_):
		# Calculate the determinant of the covariance matrix
		cov_det = np.linalg.det(covariance)

		# Calculate the volume of the n-dimensional ellipsoid for this component
		component_volume = weight * np.power(chi2_val, n_dims / 2) * np.sqrt(cov_det) * np.power(np.pi, n_dims / 2) / scipy.special.gamma(n_dims / 2 + 1)
		
		# Factor in the number of duplicate points for each component
		component_volume *= np.sum(counts)
		
		total_volume += component_volume

	return total_volume, n_components

# Function to process a single cluster, with duplicate handling
def process_cluster(cluster, df, max_GMM, seed):
	# Extract data points belonging to the current cluster
	cluster_data = df[df['cluster_label'] == cluster].drop(columns=['cluster_label']).values

	# Compute the GMM-based volume for this cluster, factoring out duplicates
	volume, n_components = estimate_GMM_volume(cluster_data, max_GMM, seed)
	
	return cluster, volume, n_components

def estimate_volumes(numeric_data, cluster_labels, num_cpus, max_GMM, seed):
	# Combine the numeric data and cluster labels into a dataframe
	df = pd.DataFrame(numeric_data)
	df['cluster_label'] = cluster_labels

	# Get unique non-noise clusters (those starting with 'c_')
	clusters = df['cluster_label'].unique()
	clusters = [c for c in clusters if str(c).startswith('c_')]

	# Process clusters in parallel using joblib
	results = Parallel(n_jobs=num_cpus)(
		delayed(process_cluster)(cluster, df, max_GMM, seed)
		for cluster in tqdm(clusters, desc="Processing Clusters", unit="clusters", unit_scale=True)
	)

	# Convert results into a DataFrame
	result_df = pd.DataFrame(results, columns=['clusterID', 'estVol', 'components'])

	# Sort by numerical order of clusterID
	result_df['clusterID_num'] = result_df['clusterID'].str.extract(r'(\d+)').astype(int)
	result_df = result_df.sort_values('clusterID_num').drop(columns=['clusterID_num'])

	# Write the output to a TSV file
	output_file = 'est_volumes.tsv'
	result_df.to_csv(output_file, sep='\t', index=False)
	print(f"Volume estimates saved to {output_file}")

	return result_df

def print_summary(clusterer, cluster_labels, volume_df):
	cluster_labels_series = pd.Series(cluster_labels)
	total_points = len(cluster_labels_series)
	non_noise_labels = cluster_labels_series[cluster_labels_series != 'noise']
	cluster_sizes = non_noise_labels.value_counts()

	n_clusters = len(cluster_sizes)
	min_cluster_size = cluster_sizes.min()
	mean_cluster_size = cluster_sizes.mean()
	max_cluster_size = cluster_sizes.max()
	median_cluster_size = cluster_sizes.median()
	mode_cluster_size = cluster_sizes.mode().values[0]
	n_noise = list(cluster_labels).count('noise')

	# Count noise points and calculate the percentage
	n_noise = total_points - len(non_noise_labels)
	noise_percentage = (n_noise / total_points) * 100

	print(f'Number of clusters: {n_clusters}')
	print(f'Min non-noise cluster size: {min_cluster_size}')
	print(f'Mean non-noise cluster size: {mean_cluster_size}')
	print(f'Mode non-noise cluster size: {mode_cluster_size}')
	print(f'Median non-noise cluster size: {median_cluster_size}')
	print(f'Max non-noise cluster size: {max_cluster_size}')
	print(f'Number of noise points: {n_noise}')
	print(f'Noise points percentage: {noise_percentage:.2f}%')

	min_vol = volume_df['estVol'].min()
	mean_vol = volume_df['estVol'].mean()
	max_vol = volume_df['estVol'].max()
	median_vol = volume_df['estVol'].median()
	mode_vol = volume_df['estVol'].mode().values[0]

	print(f'Min estimated volume: {min_vol}')
	print(f'Mean estimated volume: {mean_vol}')
	print(f'Mode estimated volume: {mode_vol}')
	print(f'Median estimated volume: {median_vol}')
	print(f'Max estimated volume: {max_vol}')

def main(file_pattern, output_file, min_cluster_size, min_samples, num_cpus, hdbscan_file, max_GMM, seed):
	# Set up logging to both terminal and file
	logger = setup_logging()

	# Set up parallel environment
	num_cpus = setup_parallel_env(num_cpus)

	logger.info('Loading embeddings...')
	data = load_data(file_pattern)
	read_ids = data.iloc[:, 0]  # Extract readID (first column)
	filenames = data['filename']  # Extract filename
	numeric_data = data.iloc[:, 1:-1]  # Numeric data (excluding readID and filename)

	if hdbscan_file:
		logger.info(f"Using precomputed HDBSCAN model from {hdbscan_file}")
		cluster_labels, clusterer = precomputed(hdbscan_file)
	else:
		logger.info('Starting HDBSCAN clustering...')
		cluster_labels, clusterer = cluster_data(numeric_data, min_cluster_size, min_samples, num_cpus)

	logger.info('Saving clustering results...')
	save_results(read_ids, filenames, cluster_labels, output_file)

	logger.info('Saving cluster abundances...')
	save_abundances(read_ids, filenames, cluster_labels)

	logger.info('Saving OTU table...')
	save_otu(read_ids, filenames, cluster_labels)

	logger.info('Estimating cluster volumes...')
	volume_df = estimate_volumes(numeric_data, cluster_labels, num_cpus, max_GMM, seed)

	logger.info('Clustering summary:')
	print_summary(clusterer, cluster_labels, volume_df)

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="Clustering and volume estimation using HDBSCAN and GMM.")
	parser.add_argument('--file_pattern', type=str, required=True, help="Pattern for input files (e.g. 'data/*.gz')")
	parser.add_argument('--output_file', type=str, required=True, help="Output file for cluster labels")
	parser.add_argument('--min_cluster_size', type=int, default=50, help="Minimum cluster size for HDBSCAN")
	parser.add_argument('--min_samples', type=int, default=10, help="Minimum samples for HDBSCAN")
	parser.add_argument('--num_cpus', type=int, default=-1, help="Number of parallel jobs")
	parser.add_argument('--HDBSCAN', type=str, default=None, help="Filepath to precomputed HDBSCAN model")
	parser.add_argument('--max_GMM', type=int, default=10, help="Maximum number of components for GMM")
	parser.add_argument('--seed', type=int, default=42, help="Random seed for GMM")
	args = parser.parse_args()

	main(args.file_pattern, args.output_file, args.min_cluster_size, args.min_samples, args.num_cpus, args.HDBSCAN, args.max_GMM, args.seed)


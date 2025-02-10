## code for computing 3D GCRs from fastq based on https://pubmed.ncbi.nlm.nih.gov/39433242/
    ## INZ, ClaudeSonnet3.5new, ChatGPT4o

import argparse
import gzip
import h5py
import numpy as np
import math
from tqdm import tqdm
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from multiprocessing import Pool, cpu_count, Manager
import os
from scipy.stats import skew

# Define the base coordinates for A, C, G, T
base_coords = {
    'A': [(0 - 0.5) * 2 * math.sqrt(1.0 / 3), (0 - 0.5) * 2 * math.sqrt(1.0 / 3), (0 - 0.5) * 2 * math.sqrt(1.0 / 3)],
    'C': [(1 - 0.5) * 2 * math.sqrt(1.0 / 3), (0 - 0.5) * 2 * math.sqrt(1.0 / 3), (1 - 0.5) * 2 * math.sqrt(1.0 / 3)],
    'G': [(0 - 0.5) * 2 * math.sqrt(1.0 / 3), (1 - 0.5) * 2 * math.sqrt(1.0 / 3), (1 - 0.5) * 2 * math.sqrt(1.0 / 3)],
    'T': [(1 - 0.5) * 2 * math.sqrt(1.0 / 3), (1 - 0.5) * 2 * math.sqrt(1.0 / 3), (0 - 0.5) * 2 * math.sqrt(1.0 / 3)],
}

def calculate_oriented_angle(a, b, v_0):
    """Calculate oriented angle between vectors a and b along orientation v_0"""
    dot_product = np.dot(a, b)
    cross_product = np.cross(a, b)
    magnitudes = np.linalg.norm(a) * np.linalg.norm(b)
    
    if magnitudes == 0:
        return 0
        
    sign = np.sign(np.dot(cross_product, v_0))
    angle = sign * np.arccos(np.clip(dot_product / magnitudes, -1.0, 1.0))
    return angle

def calculate_oriented_distance(p, q, v_0):
    """Calculate oriented distance between points p and q along orientation v_0"""
    diff_vector = p - q
    euclidean_distance = np.linalg.norm(diff_vector)
    sign = np.sign(np.dot(diff_vector, v_0))
    return sign * euclidean_distance

def extract_features(coordinates):
    """Extract 24-dimensional feature vector from CGR coordinates"""
    features = []
    orientations = [
        np.array([1, 0, 0]),
        np.array([0, 1, 0]),
        np.array([0, 0, 1])
    ]
    
    # 1. Basic statistics along each axis (9 values)
    for axis in range(3):
        features.append(np.mean(coordinates[:, axis]))
        features.append(np.var(coordinates[:, axis]))
        features.append(skew(coordinates[:, axis]))
    
    # 2. Covariances (3 values)
    cov_matrix = np.cov(coordinates.T)
    features.append(cov_matrix[0, 1])  # xy covariance
    features.append(cov_matrix[1, 2])  # yz covariance
    features.append(cov_matrix[0, 2])  # xz covariance
    
    # 3 & 4. Oriented angles and distances for each orientation
    for orientation in orientations:
        # Window size 2 angles
        angles_w2 = []
        for i in range(1, len(coordinates)-1):
            vec1 = coordinates[i] - np.zeros(3)
            vec2 = coordinates[i+1] - np.zeros(3)
            angle = calculate_oriented_angle(vec1, vec2, orientation)
            angles_w2.append(angle)
        features.append(np.mean(angles_w2))
        
        # Window size 3 angles
        angles_w3 = []
        for i in range(len(coordinates)-2):
            vec1 = coordinates[i] - coordinates[i+1]
            vec2 = coordinates[i+2] - coordinates[i+1]
            angle = calculate_oriented_angle(vec1, vec2, orientation)
            angles_w3.append(angle)
        features.append(np.mean(angles_w3))
        
        # Window size 2 distances
        distances_w2 = []
        for i in range(len(coordinates)-1):
            dist = calculate_oriented_distance(coordinates[i], coordinates[i+1], orientation)
            distances_w2.append(dist)
        features.append(np.mean(distances_w2))
        
        # Window size 3 distances
        distances_w3 = []
        for i in range(len(coordinates)-2):
            dist = calculate_oriented_distance(coordinates[i], coordinates[i+2], orientation)
            distances_w3.append(dist)
        features.append(np.mean(distances_w3))
    
    return np.array(features)

def parse_sequences(file):
    """Parse sequences from a gzipped FASTQ file."""
    return list(SeqIO.parse(gzip.open(file, "rt"), 'fastq'))

def next_coords(base, prev_coords):
    """Calculate the next coordinates in 3D CGR based on the DNA base."""
    if base not in base_coords:
        raise ValueError(f"Invalid base encountered: {base}")
    
    base_coord = base_coords[base]
    return [(prev + base_val) / 2.0 for prev, base_val in zip(prev_coords, base_coord)]

def chaos_game(sequence):
    """Compute 3D CGR coordinates for a DNA sequence."""
    coords = [[0.0, 0.0, 0.0]]
    for base in sequence:
        coords.append(next_coords(base, coords[-1]))
    return np.array(coords)

def process_record(record):
    """Process a single sequence record."""
    try:
        coords = chaos_game(str(record.seq))
        coords = np.round(coords, 6)
        features = extract_features(coords)
        return record.id, coords, features
    except Exception as e:
        return None, None, f"Error processing {record.id}: {str(e)}"

def chunks(lst, n):
    """Yield successive n-sized chunks from lst."""
    chunk_size = max(1, len(lst) // n)
    for i in range(0, len(lst), chunk_size):
        yield lst[i:i + chunk_size]

def process_chunk(chunk_data):
    """Process a chunk of records with proper initialization for parallel processing."""
    chunk, output_queue = chunk_data
    results = []
    
    # Each process gets its own h5py file handle
    temp_results = []
    for record in chunk:
        try:
            coords = chaos_game(str(record.seq))
            coords = np.round(coords, 6)
            features = extract_features(coords)
            temp_results.append((record.id, coords, features))
        except Exception as e:
            print(f"Error processing {record.id}: {str(e)}")
            continue
    
    return temp_results

def get_feature_headers():
    """Generate headers for the feature columns."""
    headers = []
    # Basic statistics
    for axis in ['x', 'y', 'z']:
        headers.extend([f'{axis}_mean', f'{axis}_var', f'{axis}_skew'])
    # Covariances
    headers.extend(['xy_cov', 'yz_cov', 'xz_cov'])
    # Oriented features
    for orientation in ['x', 'y', 'z']:
        headers.extend([
            f'{orientation}_angle_w2',
            f'{orientation}_angle_w3',
            f'{orientation}_dist_w2',
            f'{orientation}_dist_w3'
        ])
    return headers

def process_sequences_parallel(sequences, output_hdf5, output_tsv, num_cpus):
    """Process sequences using improved multiprocessing approach."""
    total_records = len(sequences)
    
    # Split sequences into chunks based on CPU count
    sequence_chunks = list(chunks(sequences, num_cpus))
    print(f"Processing {total_records} records in {len(sequence_chunks)} chunks using {num_cpus} CPUs...")
    
    # Create process pool with maxtasksperchild to prevent memory leaks
    pool = Pool(processes=num_cpus, maxtasksperchild=100)
    
    # Prepare chunk data with output queue
    from multiprocessing import Manager
    manager = Manager()
    output_queue = manager.Queue()
    chunk_data = [(chunk, output_queue) for chunk in sequence_chunks]
    
    # Write headers to TSV file
    feature_headers = get_feature_headers()
    with gzip.open(output_tsv, 'wt') as tsv_file:
        header = ['record_id'] + feature_headers
        tsv_file.write('\t'.join(header) + '\n')
        
        # Process chunks in parallel and write results
        with h5py.File(output_hdf5, 'w') as h5file:
            with tqdm(total=total_records, desc="Processing reads") as pbar:
                # Use imap_unordered for better load balancing
                for chunk_results in pool.imap_unordered(process_chunk, chunk_data):
                    for read_id, coords, features in chunk_results:
                        # Write to HDF5
                        h5file.create_dataset(f"/{read_id}", data=coords, dtype='float32')
                        
                        # Write to TSV
                        feature_str = '\t'.join(map(str, features))
                        tsv_file.write(f"{read_id}\t{feature_str}\n")
                        pbar.update(1)
    
    # Clean up
    pool.close()
    pool.join()

def main(merge_inpre, r1_inpre, r2_inpre, outpre, strand, num_cpus, paired):
    """Modified main function with improved sequence loading."""
    num_cpus = min(num_cpus, cpu_count())
    print(f"Processing with {num_cpus} CPUs...")
    
    # Use a generator approach to load sequences to reduce memory usage
    def sequence_generator():
        if merge_inpre:
            merge_infile = merge_inpre + '.fastq.gz'
            print("Loading merged sequences...")
            yield from SeqIO.parse(gzip.open(merge_infile, "rt"), 'fastq')
        
        r1_infile = r1_inpre + '.fastq.gz'
        print("Loading R1 sequences...")
        r1_sequences = list(SeqIO.parse(gzip.open(r1_infile, "rt"), 'fastq'))
        
        if paired:
            if not r2_inpre:
                raise ValueError("R2 input prefix is required for paired-end processing")
            
            r2_infile = r2_inpre + '.fastq.gz'
            print("Loading R2 sequences...")
            r2_sequences = list(SeqIO.parse(gzip.open(r2_infile, "rt"), 'fastq'))
            
            # Process paired-end reads
            for r1, r2 in zip(r1_sequences, r2_sequences):
                r2.seq = r2.seq.reverse_complement()
                combined_seq = r1.seq + r2.seq
                yield SeqRecord(
                    combined_seq,
                    id=r1.id,
                    description=f"Concatenated {r1.id} and reverse complement of {r2.id}"
                )
        else:
            yield from r1_sequences
    
    # Convert generator to list for processing
    sequences = list(sequence_generator())
    
    # Apply strand-specific transformations
    print("Applying strand transformations...")
    if strand == 'reverse':
        for record in sequences:
            record.seq = record.seq.reverse_complement()
    elif strand == 'unstranded':
        for record in sequences:
            forward_seq = record.seq
            reverse_complement_seq = record.seq.reverse_complement()
            if str(forward_seq) < str(reverse_complement_seq):
                record.seq = forward_seq
            else:
                record.seq = reverse_complement_seq
    
    # Setup output files
    output_hdf5 = outpre + '.CGR.hdf5'
    output_tsv = outpre + '.CGR.tsv.gz'
    
    # Create output directory if it doesn't exist
    os.makedirs(os.path.dirname(output_hdf5) if os.path.dirname(output_hdf5) else '.', exist_ok=True)
    
    # Process sequences with improved parallelization
    print("Generating 3D CGR coordinates and features...")
    process_sequences_parallel(sequences, output_hdf5, output_tsv, num_cpus)
    print(f"\nOutputs written to:\n{output_hdf5}\n{output_tsv}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert FASTQ sequences to 3D Chaos Game Representation coordinates and features")
    parser.add_argument('--merge_inpre', type=str, help='Input merged file prefix (optional)')
    parser.add_argument('--r1_inpre', type=str, required=True, help='Input R1 file prefix')
    parser.add_argument('--r2_inpre', type=str, help='Input R2 file prefix (required if --paired is TRUE)')
    parser.add_argument('--outpre', type=str, required=True, help='Output file prefix')
    parser.add_argument('--strand', type=str, default='forward', 
                      choices=['forward', 'reverse', 'unstranded'],
                      help='Strandedness of data (default: forward)')
    parser.add_argument('--num_cpus', type=int, default=cpu_count(), 
                      help='Number of CPUs to use for parallel processing (default: all available)')
    parser.add_argument('--paired', type=str, default='TRUE', choices=['TRUE', 'FALSE'],
                      help='Whether to process paired-end reads (default: TRUE)')

    args = parser.parse_args()
    
    # Validate arguments
    if args.paired == 'TRUE' and not args.r2_inpre:
        parser.error("--r2_inpre is required when --paired is TRUE")
    
    main(args.merge_inpre, args.r1_inpre, args.r2_inpre, args.outpre, 
         args.strand, args.num_cpus, args.paired == 'TRUE')

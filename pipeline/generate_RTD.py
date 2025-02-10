## a very rough initial go for computing RTD means from sequencing data
## heavily based on https://github.com/pandurang-kolekar/rtd-phylogeny/blob/master/src/rtd_phylogeny.py

## the goal is to take in reads, and for a given k, compute the return times for each k-mer
    ## then compute means for each k-mer for each read - normalized by read length
        ## then output an .avro file with the format:
            ## readAlias    k_mer_means

import gzip
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import numpy as np
import re
from tqdm import tqdm
import argparse
import itertools
from concurrent.futures import ProcessPoolExecutor, as_completed
import multiprocessing as mp
from functools import partial
import math

## function for getting k-mers

def getKmers(k, alphabet):
    kmers = [''.join(p) for p in itertools.product(alphabet, repeat=k)]
    return kmers

## function for computing RTDs
## by first making an rtvector and then passing it to getMean to get an RTD
## add a length correction

def getRTD(sequence, kmer):
    modSeq = re.sub(kmer, "*", sequence)
    rt = modSeq.split("*")
    rtvector = list(map(len, rt))
    if len(rtvector) > 1:
        del rtvector[0]
        del rtvector[-1]
    else:
        rtvector = []
    length = len(sequence)
    m = getMean(rtvector, length)
    return m

## function for computing mean from rtvector
## adjust to account for sequence length

def getMean(vector, length):
    if len(vector) > 0:
        mean = np.mean(vector) / length
    else:
        mean = 0
    m = [mean]
    return m

## function for loading files

def parse_sequences(file):
    return list(SeqIO.parse(gzip.open(file, "rt"), 'fastq'))

def chunks(lst, n):
    """Split list into n roughly equal chunks"""
    chunk_size = max(1, math.ceil(len(lst) / n))
    return [lst[i:i + chunk_size] for i in range(0, len(lst), chunk_size)]

def process_chunk(chunk_data):
    """Process a chunk of sequences with their kmers"""
    sequences, kmers = chunk_data
    results = []
    
    for record in sequences:
        kmer_means = {kmer: '0' for kmer in kmers}
        for kmer in kmers:
            RTD = getRTD(str(record.seq), kmer)
            kmer_means[kmer] = str(RTD[0])
        results.append((record.id, kmer_means))
    
    return results

def write_results(results_queue, output_file, kmers, total_sequences):
    """Write results to output file with progress bar"""
    with gzip.open(output_file, "wt") as fout:
        # Write header
        header = 'readID' + '\t' + '\t'.join(kmers)
        fout.write(header + '\n')
        
        # Process results as they come in
        with tqdm(total=total_sequences, desc="Processing reads") as pbar:
            completed = 0
            while completed < total_sequences:
                result_chunk = results_queue.get()
                if result_chunk is None:
                    break
                
                for record_id, kmer_means in result_chunk:
                    line = record_id + '\t' + '\t'.join(kmer_means[kmer] for kmer in kmers)
                    fout.write(line + '\n')
                
                completed += len(result_chunk)
                pbar.update(len(result_chunk))

def process_worker(chunk_data, results_queue):
    """Worker function for processing chunks"""
    try:
        results = process_chunk(chunk_data)
        results_queue.put(results)
    except Exception as e:
        print(f"Error in worker process: {str(e)}")
        results_queue.put([])

def process_sequences_parallel(sequences, kmers, num_cpus, output_file):
    """Process sequences in parallel with working progress bar and file writing"""
    # Create a manager for sharing data between processes
    manager = mp.Manager()
    results_queue = manager.Queue()
    
    # Split sequences into chunks
    sequence_chunks = chunks(sequences, num_cpus * 4)  # Create more chunks than CPUs for better load balancing
    chunk_data = [(chunk, kmers) for chunk in sequence_chunks]
    
    # Start the writer process
    writer_process = mp.Process(
        target=write_results,
        args=(results_queue, output_file, kmers, len(sequences))
    )
    writer_process.start()
    
    # Process chunks using a process pool
    with ProcessPoolExecutor(max_workers=num_cpus) as executor:
        # Submit all chunks for processing
        futures = [
            executor.submit(process_worker, chunk, results_queue)
            for chunk in chunk_data
        ]
        
        # Wait for all futures to complete
        for future in futures:
            future.result()  # This ensures we catch any exceptions
    
    # Signal completion to writer process
    results_queue.put(None)
    writer_process.join()

def main(k, alphabet, merge_inpre, r1_inpre, r2_inpre, outpre, strand, num_cpus, paired):
    # Set up multiprocessing
    mp.set_start_method('spawn', force=True)
    
    print(f"Initializing with {num_cpus} CPUs...")
    
    # Initialize sequences list
    sequences = []
    
    # Load merged reads if provided
    if merge_inpre:
        merge_infile = merge_inpre + '.fastq.gz'
        print("Loading merged sequences...")
        sequences.extend(parse_sequences(merge_infile))
    
    # Load and process R1 reads
    r1_infile = r1_inpre + '.fastq.gz'
    print("Loading R1 sequences...")
    r1_sequences = parse_sequences(r1_infile)
    
    if paired:
        # Load and process R2 reads if in paired mode
        if not r2_inpre:
            raise ValueError("R2 input prefix is required for paired-end processing")
            
        r2_infile = r2_inpre + '.fastq.gz'
        print("Loading R2 sequences...")
        r2_sequences = parse_sequences(r2_infile)
        
        # Process paired-end reads
        print("Processing paired reads...")
        for r1, r2 in zip(r1_sequences, r2_sequences):
            r2.seq = r2.seq.reverse_complement()
            combined_seq = r1.seq + r2.seq
            combined_record = SeqRecord(
                combined_seq,
                id=r1.id,
                description=f"Concatenated {r1.id} and reverse complement of {r2.id}"
            )
            sequences.append(combined_record)
    else:
        sequences.extend(r1_sequences)
    
    # Apply strand transformations
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
    
    # Generate kmers
    kmers = getKmers(k, alphabet)
    print(f"Generated {len(kmers)} {k}-mers")
    
    # Setup output file
    outfile = outpre + '.RTD.gz'
    
    # Process sequences
    total_sequences = len(sequences)
    print(f"Processing {total_sequences} sequences...")
    process_sequences_parallel(sequences, kmers, num_cpus, outfile)
    print(f"\nOutput written to: {outfile}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compute mean and std RTDs from FASTQ file.")
    parser.add_argument('--k', type=int, default=2, help='Length of k-mers')
    parser.add_argument('--alphabet', type=str, nargs='+', default=['A', 'C', 'G', 'T'], help='Alphabet to use for k-mers')
    parser.add_argument('--merge_inpre', type=str, help='Input merged file prefix (optional)')
    parser.add_argument('--r1_inpre', type=str, required=True, help='Input R1 file prefix')
    parser.add_argument('--r2_inpre', type=str, help='Input R2 file prefix (required if --paired is TRUE)')
    parser.add_argument('--outpre', type=str, required=True, help='Output file prefix')
    parser.add_argument('--strand', type=str, default='forward', help='Strandedness of data, forward, reverse, or unstranded (= canonical counting)')
    parser.add_argument('--num_cpus', type=int, default=1, help='Number of CPUs to use for parallel processing')
    parser.add_argument('--paired', type=str, default='TRUE', choices=['TRUE', 'FALSE'], 
                      help='Whether to process paired-end reads (default: TRUE)')

    args = parser.parse_args()
    
    # Validate arguments
    if args.paired == 'TRUE' and not args.r2_inpre:
        parser.error("--r2_inpre is required when --paired is TRUE")
    
    main(args.k, args.alphabet, args.merge_inpre, args.r1_inpre, args.r2_inpre, 
         args.outpre, args.strand, args.num_cpus, args.paired == 'TRUE')

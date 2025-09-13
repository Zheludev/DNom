## a very rough initial go for computing RTD means from sequencing data
	## heavily based on https://github.com/pandurang-kolekar/rtd-phylogeny/blob/master/src/rtd_phylogeny.py

## the goal is to take in reads, and for a given k, compute the return times for each k-mer
	## then compute means for each k-mer for each read - normalized by read length
		## then output an .avro file with the format:
			## readAlias	k_mer_means

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

## function to process a single sequence and compute RTDs

def process_sequence(record, kmers):
	kmer_means = {kmer: '0' for kmer in kmers}
	for kmer in kmers:
		RTD = getRTD(str(record.seq), kmer)
		kmer_means[kmer] = str(RTD[0])
	return record.id, kmer_means

## main function

def main(k, alphabet, merge_inpre, r1_inpre, r2_inpre, outpre, strand, num_cpus):
	merge_infile = merge_inpre + '.fastq.gz'
	r1_infile = r1_inpre + '.fastq.gz'
	r2_infile = r2_inpre + '.fastq.gz'
	outfile = outpre + '.RTD.gz'
	
	## compute k-mers
	kmers = getKmers(k, alphabet)
	
	## count number of reads in the input files
	with gzip.open(merge_infile, "rt") as f:
		num_reads = sum(1 for _ in SeqIO.parse(f, 'fastq'))
	with gzip.open(r1_infile, "rt") as f:
		num_reads = num_reads + sum(1 for _ in SeqIO.parse(f, 'fastq'))
	
	## load files
	if merge_inpre:
		print("Loading merged sequences...")
		merge_sequences = parse_sequences(merge_infile)
	
	if r1_inpre and r2_inpre:
		print("Loading R1 and R2 sequences...")
		r1_sequences = parse_sequences(r1_infile)
		r2_sequences = parse_sequences(r2_infile)
	
	## fix orientation of r2 reads
	print("Processing R2 reads...")
	r2_reverse_complement = []
	for record in r2_sequences:
		record.seq = record.seq.reverse_complement()
		r2_reverse_complement.append(record)
	
	## concatenate the R1 and (reverse complemented R2) reads
	print("Concatenating R1 and R2 reads...")
	concatenated_sequences = []
	for r1, r2 in zip(r1_sequences, r2_reverse_complement):
		combined_seq = r1.seq + r2.seq
		combined_record = SeqRecord(
			combined_seq,
			id=r1.id,
			description=f"Concatenated {r1.id} and reverse complement of {r2.id}"
		)
		concatenated_sequences.append(combined_record)
	
	## Combine all sequences in the desired order
	sequences = merge_sequences + concatenated_sequences
	
	## Apply transformations based on the strand variable
	print("Applying strand transformations...")
	if strand == 'forward':
		## Keep sequences unchanged
		pass
	elif strand == 'reverse':
		## Transform to reverse complement
		for record in sequences:
			record.seq = record.seq.reverse_complement()
	elif strand == 'unstranded':
		## Transform to canonical representation
		for record in sequences:
			forward_seq = record.seq
			reverse_complement_seq = record.seq.reverse_complement()
			if str(forward_seq) < str(reverse_complement_seq):
				record.seq = forward_seq
			else:
				record.seq = reverse_complement_seq
	
	## process sequences in parallel and write output
	print("Generating RTDs...")
	print(f"Processing {num_reads} reads with {num_cpus} CPUs...")
	
	with gzip.open(outfile, "wt") as fout:
		## Write column headers
		header = 'readID' + '\t' + '\t'.join(kmers)
		fout.write(header + '\n')
		
		## Create a ProcessPoolExecutor
		with ProcessPoolExecutor(max_workers=num_cpus) as executor:
			futures = {executor.submit(process_sequence, record, kmers): record.id for record in sequences}
			for future in tqdm(as_completed(futures), desc="Processing reads", total=num_reads):
				record_id, kmer_means = future.result()
				line = record_id
				for kmer in kmers:
					line += '\t' + kmer_means[kmer]
				fout.write(line + '\n')

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="Compute mean and std RTDs from FASTQ file.")
	parser.add_argument('--k', type=int, default=2, help='Length of k-mers')
	parser.add_argument('--alphabet', type=str, nargs='+', default=['A', 'C', 'G', 'T'], help='Alphabet to use for k-mers')
	parser.add_argument('--merge_inpre', type=str, required=True, help='Input merged file prefix')
	parser.add_argument('--r1_inpre', type=str, required=True, help='Input R1 file prefix')
	parser.add_argument('--r2_inpre', type=str, required=True, help='Input R2 file prefix')
	parser.add_argument('--outpre', type=str, required=True, help='Output file prefix')
	parser.add_argument('--strand', type=str, default='forward', help='Strandedness of data, forward, reverse, or unstranded (= canonical counting)')
	parser.add_argument('--num_cpus', type=int, default=1, help='Number of CPUs to use for parallel processing')
	args = parser.parse_args()
	main(args.k, args.alphabet, args.merge_inpre, args.r1_inpre, args.r2_inpre, args.outpre, args.strand, args.num_cpus)

# `denim`
*reference-free, differential abundance analysis of NGS data in a `snakemake` pipeline*

# overview
`denim` is an end-to-end [`snakemake`](https://snakemake.readthedocs.io/en/stable/) pipeline that <ins>takes in:</ins>
```
1. NGS sequencing data from multiple (arbitrary number) samples

      and

2. a metadata table describing these samples
```
and <ins>produces:</ins>
```
1. clusters of reads that are significantly differentially present across a metadata parameter of the samples

      and

2. volcano plots summarizing this differential presence
```
### <ins>*__critically, `denim` does not require any reference (such as a database or genome) to perform this analysis__*</ins>

# pipeline
at a high level, `denim` achieves this by:

1. [embedding](https://www.ibm.com/think/topics/embedding) all the NGS reads (by sample) into numerical vectors
   - *the (Euclidian) distance between these vectors is related to the sequence similarity between their underlying reads*

3. performing all-versus-all (across samples) clustering of these vectors
   - *essentially grouping all the reads based on sequence similarity*
  
4. using this clustering for [`Sleuth`](https://pachterlab.github.io/sleuth/about)-based statistical analysis
   - *performed over a metadata characteristic of the samples*

6. extracting the read clusters that satisfy (user defined) significance thresholds

the user can then analyse these reads as they see fit with the knowledge that according to `denim`, something is interesting about them.

# details
__preprocessing__
1. paired-end reads are first quality filtered and then merged (if overlaps are present) using [`fastp`](https://github.com/OpenGene/fastp)

__embedding__

2. these pre-processed reads are then embedded using:
   - [return time distributions](https://pubmed.ncbi.nlm.nih.gov/22820020/) (RTDs)
     
       *and/or*
     
   - [3D chaos game representations](https://pubmed.ncbi.nlm.nih.gov/39433242/) (CGRs)
  
   *both methods deterministically convert the reads into vectors that preserve a notion of what sequence begat them*

__clustering__

3. each embedded dataset is then clustered (across all samples simultaneously) using [`fast_HDBSCAN`](https://github.com/TutteInstitute/fast_hdbscan)
   - `fast_HDBSCAN` is able to identify saliant clusters of vectors while also being able to identify noise
   - an estimation of cluster hypervolume is also performed for `Sleuth` analysis

__merging__

4. *optionally*, if more than one embedding was used, their resulting clusterings can be merged based on read/vector overlap

__statistical analysis *via* `Sleuth`__

5. the results of the clustering(s) are then converted into file formats that `Sleuth` can accept for analysis
   - both [Wald tests](https://en.wikipedia.org/wiki/Wald_test) and [Likelihood-ratio tests](https://en.wikipedia.org/wiki/Likelihood-ratio_test) are performed
   - the user defines the *q-value* and *fold-change* thresholds required for signficance as well as the metadata parameter over which the tests are performed
   - automated [Volcano plotting](https://en.wikipedia.org/wiki/Volcano_plot_(statistics)) is also attempted

__retrieval__

6. if significantly differentially present clusters are identified, the corresponding reads are extracted for *post-hoc* analysis

# design
`denim` is purposefully written in a modular way. Each step (*preprocessing*, *embedding*, *clustering*, *statistical analysis*, and *retrieval*) can be altered, swapped-out, and added-to from the modular `snakemake` pipeline.
Additionally, `denim` makes use of CPU parallelization wherever possible for acceleration while also trying to limit its RAM footprint.

# promises
I make no promises about the robustness, statistical-validity, or longevity of `denim`. At this point, I see this tool as more of an intellectual exercize than anything else.

# future
I have no current plans to expand on `denim`, however, implementation of hardware (GPU) acceleration, would likely speed up the whole pipeline while also opening up opportunities for more exotic embedding and clustering approaches.

alternative *embedding* strategies could involve:
- LLMs like [`DNABERT-S`](https://arxiv.org/abs/2402.08777)

alternative *clustering* strategies could involve:
- [locality sensitive hashing](https://www.pinecone.io/learn/series/faiss/locality-sensitive-hashing-random-projection/)
- [GPU-accelerated `HDBSCAN`](https://developer.nvidia.com/blog/faster-hdbscan-soft-clustering-with-rapids-cuml/)

alternative *statistical analysis* strategies could involve:
- [`LinDA`](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-022-02655-5)
- [`MaAsLin3`](https://www.biorxiv.org/content/10.1101/2024.12.13.628459v1)

# alternatives to `denim`
a number of folks smarter than me have already approached this problem. fundementally, these approaches rely on `k-mer`-based analysis (`k-mer` counting / assembly):
- [`IMP`](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-1116-8) (assembly-based)
- [`kmdiff`](https://academic.oup.com/bioinformatics/article/38/24/5443/6782954?login=false) (`k-mer` counting)
- [`KMerFineCompare`](https://github.com/FireLabSoftware/KMerFineCompare) (`k-mer` counting)
- [`SPLASH`](https://www.cell.com/cell/fulltext/S0092-8674(23)01179-0) (`k-mer`-pair counting)

# installation
currently, installation requires you to manually create a new [`mamba`](https://mamba.readthedocs.io/en/latest/) environment (this will take a while):
   
   `mamba env create -f denim_mamba_env.yml`

# running `denim`

__setup__

`denim` requires a tab-delimited metadata file that must have the columns: `SRA`, `description`, and `sample`. A fourth column is likely needed for your statistical test (e.g. `strain`).

`metadata.tsv` example:
```
SRA	description	sample	strain	clone
SRR29668074	ObP1_C_RNaseR	ObP1_C	ObP1	C
SRR29668075	ObP1_B_RNaseR	ObP1_B	ObP1	B
SRR29668076	ObP1_A_RNaseR	ObP1_A	ObP1	A
SRR29668077	ObN1_C_RNaseR	ObN1_C	ObN1	C
SRR29668078	ObN1_B_RNaseR	ObN1_B	ObN1	B
SRR29668079	ObN1_A_RNaseR	ObN1_A	ObN1	A
```

raw, paired-end reads must be in a directory titled `00_input` within corresponding sub-directories titled after rows in `SRA`. These `SRA` sub-directories must contain reads titled `SRA_[12].fastq.gz`

`00_input` example:
```
00_input/
├── SRR29668074
│   ├── SRR29668074_1.fastq.gz
│   └── SRR29668074_2.fastq.gz
├── SRR29668075
│   ├── SRR29668075_1.fastq.gz
│   └── SRR29668075_2.fastq.gz
├── SRR29668076
│   ├── SRR29668076_1.fastq.gz
│   └── SRR29668076_2.fastq.gz
├── SRR29668077
│   ├── SRR29668077_1.fastq.gz
│   └── SRR29668077_2.fastq.gz
├── SRR29668078
│   ├── SRR29668078_1.fastq.gz
│   └── SRR29668078_2.fastq.gz
└── SRR29668079
    ├── SRR29668079_1.fastq.gz
    └── SRR29668079_2.fastq.gz
```

you should also configure the `config.yaml` file to run per your requirements. because of the `snakemake` architecture, you can re-run partially/completely run pipelines with alterations to the `config.yaml` file at a later point.

`config.yaml` example:
```
# Initialization
base_dir: "."
metadata: "metadata.tsv"
pipeline_dir: "<absolute path to pipeline dir>"

# Pipeline mode
embedding: "merge"  # can be: RTD, CGR, or merge

# Computational parameters
threads: 64

# Embedding paramters
strand: "forward"

# RTD parameters
k: 2
alphabet: "A C G T"

# Clustering parameters
min_cluster_size: 400
min_samples: 5
max_GMM: 5
seed: 42

# Merging parameters
overlap: 10

# Bootstrapping parameters
n_bootstraps: 1000

# Sleuth parameters
parameter: "strain"
pseudocount: 0.5
alpha: 0.05
foldchange: 10
```

__actually running `denim`__

*be sure to use the absolute path to the pipeline directory both in command-line command and the `snakemake` config file*

```
# navigate to base directory
cd <base dir containing 00_input>

# create pipeline directory path variable
pipeline=$(echo '<absolute path to pipeline dir>')

# dry run test
# be sure to set '--cores 1' as multitheading is handled within the config.yaml
snakemake --snakefile "$pipeline"/denim.sm --configfile config.yaml --cores 1 -n

# run denim
snakemake --snakefile "$pipeline"/denim.sm --configfile config.yaml --cores 1
```

__a note on analytical design__

`denim` is written for trivial analytical designs with only one variable (e.g. `~strain`). in reality, you often want to account for other variables in your analysis (e.g. `~strain + sex + tissue`). Writing one 'one-size-fits-all' solution felt clunky, so instead, I've included an example [`jupyter`](https://docs.jupyter.org/en/latest/) notebook (`sleuth_plotting_notebook.ipynb`) to assist in this more nuanced analysis.

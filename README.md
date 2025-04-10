# nascent RNA simulator

These scripts are designed for _in silico_ generation of mRNAs or long/short-read RNA-sequencing datasets for use in optimization and validation of nascent RNA analyses. 

In brief, these simulations model transcription elongation for genes of interest and subsequent RNA-sequencing interrogation of nucleotide recoding- and enrichment-based techniques. User-defined simulation parameters are tunable to enable modeling of a variety of experimental conditions and allow for robust validation of downstream nascent RNA analysis pipelines.

### Clustering by gene features

We decided to simulate reads from the human genome (rather than simulating random sequences) so that the nucleotide content is more authentic. seq_and_clustering.R will analyze a variety of metrics (like # exons, # introns, exon/intron/UTR length, nucleotide composition, etc.) for all annotated genes. The longest transcript per gene will be used. Then it will perform hierarchical clustering on the metrics to define groups of genes/transcripts with similar metrics. It outputs a gtf(ish) file where each row is a feature (exon or intron) and the last column contains the sequence for it

Input
`-g, --gtf
-f,--genome_fasta
-o; --dir_out
-t; --threads 
-c; n_gene_clusters
-n; number_of_genes
`

Usage example:

`Rscript seq_and_clustering.R -g Homo_sapiens.GRCh38.95.gtf -f Homo_sapiens.GRCh38.dna.primary_assembly.fa -o out_dir`

Output:
metrics_per_transcript_df.tsv -> containing different metrics for the longest annotated transcripts per gene
longest_transcript_and_features.bed -> coordinates for the longest transcript, its introns and exons
longest_transcript_and_features.fasta -> sequence for the longest transcript, its introns and exons
A tsv file per gene where each row is a feature (exon or intron) and the last column contains the sequence for it


### Assigning an RNA Pol II elongation rate per nucleotide

This script will assign a rate to each nucleotide in the transcript.


Input
`--tsv
--region_size_range
--elong_rate_range
--pause_elong_rate
--pause_occur_chance
--o
--flat_rates
--gene_level`

Usage example
`python rates_per_region.py --tsv ENSG00000237672.tsv --region_size_range 2,100 --pause_occur_chance 0 --o .`

Output:
A rate_per_gene directory containing a RatesandTraversalTimes and a VariableElongationRateRegions file per gene. RatesandTraversalTimes contains the elongation rate per nucleotide while VariableElongationRateRegions contains the elongation rate variability information per region.

### mRNA generation

mRNAs can be generated in two different ways using mRNA_generator_monolabel.py or mRNA_generator_kinetic_barcoding.py. The user can opt to simulate a DRB-treated (all RNA Pol II are synchronized at the TSS at the beginning of the experiment) labeling or no DRB. For the monolabel, it can be specified if and which nucleotide substitution is induced.

### Read generation

This pipeline can simulate both short and long-read sequencing based on the generated mRNA molecules. The scripts for that are fastq_generator_long_read.py and fastq_generator_short_read.py that also offer options to the desired library preparation and strandness type.

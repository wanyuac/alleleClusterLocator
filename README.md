# Validating Inferred Co-transferred Allele Clusters through Nucleotide-sequence Clustering

Yu Wan

23 Aug 2020



This repository consists of a helper script (`acl.py`) of R package [GeneMates](https://github.com/wanyuac/GeneMates) and the script is developed for validation of interred horizontally co-transferred allele clusters. This Python script locates co-localised alleles in contigs of genome assemblies, extracts genomic regions that harbour all the alleles, and clusters sequences of these regions. More specifically, it carries out the following major steps:

1. Extracting allele sequences from a FASTA file in accordance with clusters defined in a TSV file.
2. Searching for these allele sequences in several genome assemblies (FASTA format) using megaBLAST or BLASTn.
3. Concatenating BLAST results and extracting perfect matches (namely, exact hits, which display 100% nucleotide identity and query coverage).
4. For every cluster in each genome, determining the smallest region (called a minimal region) harbouring all member alleles in a contig when all alleles have exact matches in the same contig.
5. Extracting these minimal regions from contigs and performing sequence clustering with CD-HIT-EST.



Arguments of this script:

```bash
python acl.py --help
usage: acl.py [-h] --clusters CLUSTERS --allele_db ALLELE_DB --assemblies ASSEMBLIES
              [ASSEMBLIES ...] [--suffix SUFFIX] [--outdir OUTDIR] [--prefix PREFIX]
              [--makeblastdb MAKEBLASTDB] [--blast BLAST] [--identity IDENTITY]
              [--algorithm ALGORITHM] [--cdhit CDHIT] [--cdhit_args CDHIT_ARGS] [--skip] [--clean]

Searching allele clusters in contigs

optional arguments:
  -h, --help            show this help message and exit
  --clusters CLUSTERS   A two-column tab-delimited file defining allele clusters
  --allele_db ALLELE_DB
                        A FASTA file for allele sequences
  --assemblies ASSEMBLIES [ASSEMBLIES ...]
                        FASTA files for assembled contigs
  --suffix SUFFIX       Characters to be chopped off from the end of assembly filename in order to
                        get a sample name
  --outdir OUTDIR       Output directory
  --prefix PREFIX       A prefix for output files
  --makeblastdb MAKEBLASTDB
                        Path to call makeblastdb
  --blast BLAST         Path to call BLAST
  --identity IDENTITY   The minimum nucleotide identity used for calling a hit.
  --algorithm ALGORITHM
                        BLAST algorithm (the -task argument)
  --cdhit CDHIT         Path to call CD-HIT-EST
  --cdhit_args CDHIT_ARGS
                        A string of arguments for CD-HIT-EST
  --skip                Flag it to skip existing output files
  --clean               Flag it to delete original BLAST output files
```

Particularly, columns of the TSV file for argument `--clusters` are:

- Cluster ID
- Comma-delimited vector of allele names (sequence IDs in the input FASTA file specified by argument `--allele_db`)



This script requires BioPython to run. Any issue report or suggestion for further development will be welcomed.



**Citation**

Wan, Y., Wick, R. R., Zobel, J., Ingle, D. J., Inouye, M., & Holt, K. E. (2020). GeneMates: an R package for Detecting Horizontal Gene Co-transfer between Bacteria Using Gene-gene Associations Controlled for Population Structure. *BioRxiv*, 2020.02.29.970970. https://doi.org/10.1101/2020.02.29.970970.

# Genotube v3.0
Version : Tokyo Drift

<https://github.com/adrienlemeur/genotube-td/>

Genotube is a genomic workflow for the study of Mycobacterium tuberculosis and affiliated species (only MTBC supported for now). It includes download, trimming, alignement, variant calling, lineage caracterrisation and antibioresistance predition for M. tuberculosis.

## Installation
Genotube requires :
- Nextflow (> 21.0)
- Singularity (or Docker) (any recent version should do the trick)
- A Linux AMD64 architecture

Both Nextflow and Singularity can be installed with Conda or Mamba. Docker must be installed with root rights.


## Documentation

[Genotube wiki](https://github.com/adrienlemeur/genotube/wiki) (coming later)


## Quick start

**Running Genotube** :

Command :
```
nextflow main.nf \\
  --sra data/SRA.txt \\
  --referenceSequence data/mtuberculosis/genomes/reference/H37Rv.fa.gz \\
  --referenceAnnotation data/mtuberculosis/genomes/reference/H37Rv.fa.gz
```

**MAIN ARGUMENTS**
**INPUT**
* --sra : path to a line separated file with NCBI run ID to be download. Default : a small dataset of 25 run ID of strains from different lineage
* --referenceSequence : reference sequence for mapping
* --referenceAnnotation : GFF3 annotation file for the provided reference sequence
* --results : Folder where intermediary and output files are stored. Genotube will use these files to resume the analysis if interrupted.
> NOTE : SRA WILL NOT be redownloaded if either raw fastq, trimmed fastq, raw bam, bam, raw vcf or annotated vcf are available. To force re-download, use --FORCE_DOWNLOAD flag

**IMPORTING FILES**
* --fasta : path to a folder with fasta samples. Raw reads will be simulated from the fasta and input in Genotube.
* --fastq : path to a folder with compressed fastq. Paired fastq are discriminated based on the presence of "_1" and "_2" strings.
* --bam : path to a folder with bam alignement files.
* --vcf : path to a folder with compressed vcf.

* **FORCE AND SKIP OPTIONS**
* --FORCE ['DOWNLOAD', 'TRIM', 'MAP', 'CALL' and/or 'ANN'] : step to force
* --SKIP ['COVERAGE_CHECK', 'PROFILING'] : list of step to not perform
* --REMOVE	['FASTQ', 'BAM', 'VCF'] : type of files not to be kept

**GENE ORIENTED**
* --target_region : bed file specifying specific regions to analyse
> NOTE : Make sure bed file CHROM column matches reference annotation.

**CONTAMINATION CHECK**
* --contam_check ['taxo_class', 'compet_mapping', 'none'] : Method to use to process contamination. Read taxonomic classification with Kraken2, competitive mapping against reference genomes from common contaminants.

**VARIANT CALLING OPTIONS**
* --species ['MTBC', 'other'] : species to consider
* --variant_caller ['gatk', 'samtools' or 'freebayes'] : variant caller to use
* --deletion_region : deletion region caller. It returns a vcf file as output.

**PHYLOGENY**
* --tree_build ['none', 'classic' or 'fastest'] : Phylogenetic tree reconstruciton method (max. likelihood). Classic with RAxML, fastest with fasttree. Default : none
* --outgroup : Sample name to root the tree
* --taxonomy ['barcode' or 'placement'] : strain taxonomical assignation. Default : 'barcode' for MTBC, none for other species.

**PHYLOGENETICAL PLACEMENT**
Do not change unless you want to use you own reference samples / taxonomical ranks
* --referenceMSA : MSA file with full genome sequence of reference samples
> Note : the MSA width must match the size of the reference genome ! (by default : H37Rv strain)
* --referenceTree : Newick reference tree for phylogenetical placement. Must match the MSA file
* --referenceModel : RAxML evolution model of the reference tree
* --referenceTaxo : reference strain taxonomy

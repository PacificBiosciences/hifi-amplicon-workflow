# HiFi Amplicon Workflow
Generalized amplicon workflow for single-gene HiFi amplicon targets

## Workflow Outline
  1. Demultiplex by barcode with [lima](https://lima.how/)
  2. Verify and trim primer sequences with [lima](https://lima.how/)
  3. Convert demuxed HiFi bam files to fastq + index
  4. Cluster with [pbaa](https://github.com/PacificBiosciences/pbAA)
  5. Align cluster consensus & labeled HiFi reads
  6. Call variants per cluster with htslib
  7. Annotate variants with [clinVar](https://www.ncbi.nlm.nih.gov/clinvar/) 
  8. Annotate variants with [bcftools consequence](https://academic.oup.com/bioinformatics/article/33/13/2037/3000373)
  9. Generate variant reports 

# To set up and run pipeline
```
# Note that the workflow is set up to demo GBA/GBAP1.  
# Update the config.yaml and support files as described below for other targets

# create the base conda environment
conda create \
  --channel bioconda \
  --channel conda-forge \
  --prefix ./conda_env \
  python=3.9 snakemake pysam pandas openpyxl mamba lockfile

# activate the base conda environment
conda activate ./conda_env

# clone the github repo into 'workflow' subdirectory
git clone https://github.com/PacificBiosciences/hifi-amplicon-workflow-poc.git workflow

# create a couple directories for reference sequence and analysis logs
mkdir reference cluster_logs

# REFERENCE
# drop your reference.fasta, reference.fasta.fai, barcode.fasta and primer.fasta into the reference directory
# For simplicity, only include chromosomes of your targets ( samtools faidx GRCH38.fasta chr1 > GRCh38_chr1.fasta )
# adjust the paths to those files in workflow/config.yaml

# BARCODES and PRIMERS
# drop your barcode.fasta and primer.fasta into the workflow/data directory
# adjust the paths to those files in workflow/config.yaml
# Set 'barcodePreset' in workflow/config.yaml to match your barcode design (see lima link above for help)
# Set 'demuxCross' to True if your primers have different barcodes and you expect large deletions to amplify cross-primer

# Samples can have more than one expected barcode.  
# Simply list all expected barcode combinations with the *same sample name* in the biosample.csv

# TARGETS
# Define the amplified regions in workflow/config.yaml relative to the reference -- this region will be used for the cluster guide
# GBA/GBAP1 are included as an example.  Remove these if they are not your targets!
# (NOTE -- for very highly homologous [>99%] loci like SMN1 & SMN2 it may be desirable to define only one region)
# See config.yaml for example coordinates.  Outer coordinates of primer sequences are a good option.
# Any whole number (1..?) of named targets is allowed for multiplexed assays -- just add the name: coord pair under 'regions'

# ANNOTATION
# clinvar is included for convenience.
# Consequence gff files can be found here: ftp://ftp.ensembl.org/pub/current_gff3/homo_sapiens
# Replace chr1 gff3 with appropriate file for your target(s) and update workflow/config.yaml

# CLUSTER
# Adjust cluster settings in workflow/config.yaml. See pbaa linka above for help on parameter definitions.

# run the workflow for batch <batch_name>
sbatch workflow/run_snakemake.sh <batch_name> <biosample_csv> <hifi_reads>
 
Results will be generated in a new directory named batches/<batch_id>
```

# Directory structure within basedir

```text
.
├── cluster_logs  # slurm stderr/stdout logs
├── reference
│   ├── reference.fasta
│   └── reference.fasta.fai
├── batches
│   └── <batch_name>  # batch_id
│       ├── demux/  # demultiplexed hifi reads
│       ├── logs/  # per-rule stdout/stderr logs
│       ├── pbaa/  # extracted reference sequences used for clustering
│       ├── biosample.csv # used for demux, may include within-sample cross-products if demuxCross=True
│       ├── <sample_id 1>/  # per-sample results, one for each sample
│       :        ...
│       └── <sample_id n>/  # per-sample results, one for each sample
│           ├── 0152.GRCh38.clinvar_annotated.variant_summary.xlsx # variant summary excel sheet
│           ├── clusterQC.report.tsv # Barcode and Primers -- check here for cross-amplified calls
│           ├── hifi.GRCh38.painted.bam # aligned HiFi subset with cluster labels
│           ├── hifi.GRCh38.painted.bam.bai
│           ├── pbaa_0152_GBA_ecr.fasta # error-corrected reads for target GBA
│           ├── pbaa_0152_GBAP1_ecr.fasta # error-corrected reads for target GBAP1
│           ├── pbaa_failed_cluster_sequences.fasta # failed cluster consensus
│           ├── pbaa_passed_cluster_sequences.fasta # passed cluster consensus
│           ├── pbaa_passed_cluster_sequences.GRCh38.bam # cluster consensus aligned to reference
│           ├── pbaa_passed_cluster_sequences.GRCh38.bam.bai
│           ├── pbaa_read_info.txt # per-read clustering info
│           ├── primertrim/  # trimmed/demux input data.  All HiFi reads for this sample
│           ├── primertrim_subsample/  # subsampled fastq, input to pbaa
│           └── vcf/ # vcf/tsv reports for each cluster identified in pbaa clustering
└── workflow  # clone of this repo
    ├── <workflow_files>
    │       :        ...
    ├── <workflow_files>
    └── data
        ├── clinvar_fixnames.vcf.gz
        ├── clinvar_fixnames.vcf.gz.tbi
        ├── primers_gba_gbap1.fasta
        ├── Sequel_96_barcodes_v1.fasta
        └── Homo_sapiens.GRCh38.108.chromosome.1.gff3.gz

```

# Release History
* 0.1.0 - Initial Release (11/16/2022)
* 0.1.1 - Add alignment & color by cluster of HiFi reads
* 0.1.2 - Add annotation, reporting, barcode design option
* 0.1.3 - Add consequence calling of variants, collated excel file

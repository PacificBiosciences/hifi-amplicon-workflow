# HiFi Amplicon Workflow
Generalized amplicon workflow for single-gene HiFi amplicon targets

# To set up and run pipeline
```
# create the base conda environment
conda create \
  --channel bioconda \
  --channel conda-forge \
  --prefix ./conda_env \
  python=3 snakemake=6.15.5 pysam mamba lockfile

# activate the base conda environment
conda activate ./conda_env

# clone the github repo into 'workflow' subdirectory
git clone http://bitbucket.pacificbiosciences.com:7990/scm/~jharting/single_gene_amplicon.git workflow

# create a couple directories for reference sequence and analysis logs
$ mkdir reference cluster_logs

# drop your reference.fasta, reference.fasta.fai, barcode.fasta and primer.fasta into the reference directory
# adjust the paths to those files in workflow/config.yaml
# Define the amplified regions in workflow/config.yaml relative to the reference -- this region will be used for the cluster guide
# GBA/GBAP1 are included as an example.  Remove these if they are not your targets!
# (NOTE -- for very highly homologous [>99%] loci like SMN1 & SMN2 it may be desirable to define only one region)

# run the workflow for batch <batch_name>
$ sbatch workflow/run_snakemake.sh <batch_name> <biosample_csv> <hifi_reads>
 
Results will be generated in a new directory named batches/<batch_id>
```

# Release History
* 0.1.0 - Initial Release (11/16/2022)

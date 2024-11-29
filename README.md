# Clonal hematopoiesis in metastatic urologic malignacies
This repo contains the relevant snakemake rules for to carry out FastQ modification, BAM generation and variant calling as explained in our article "Clonal hematopoiesis in metastatic urothelial and kidney cancer". 

## Pre-requisites
Our analysis pipeline is written in Python using the Snakemake workflow management system. Please follow these instructions for setup:
1. All dependencies needed to run the pipeline are provided in the `envs/snakemake.yaml` file. You can create the snakemake environment by running `conda env create -f snakemake.yaml`
2. GATK's base calibration tool that we use requires the following three files (`resources_broad_hg38_v0_Homo_sapiens_assembly38.known_indels.vcf`, `resources_broad_hg38_v0_Mills_and_1000G_gold_standard.indels.hg38.vcf`, `resources_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf`) which can be downloaded from this [Google Cloud](https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0;tab=objects?prefix=&forceOnObjectsSortingFiltering=false&pli=1) link. 
3. Mutect2 requires a panel of normals, which can be obtained from this [link](https://gatk.broadinstitute.org/hc/en-us/articles/360035890631-Panel-of-Normals-PON).

## Running the workflow
1. Edit the `config.yaml` file with the file locations of relevant files and directories.
2. Workflow expects raw FastQ files to be placed into the `results/data/fastq` directory. All files placed here will be processed by the pipeline.
3. The pipeline can be evoked with the `snakemake` command. By issuing `snakemake -n` you can issue a dry run.  
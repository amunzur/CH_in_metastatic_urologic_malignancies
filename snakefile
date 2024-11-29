import pandas as pd
import re
import os

configfile: "config/config.yaml"
exec(open("workflow/snakemake_functions.py").read())

# References
PATH_hg38 = config["PATH_hg38"]
PATH_hg38_dict = config["PATH_hg38_dict"]
PATH_bed = config["PATH_bed"]
PATH_bed_intervals = config["PATH_bed_intervals"]
PATH_baits = config["PATH_baits"]

# For base recalibration
PATH_known_indels = config["PATH_known_indels"]
PATH_gold_std_indels = config["PATH_gold_std_indels"]
PATH_SNP_db = config["PATH_SNP_db"]

# For Mutect2
PATH_PoN = config["PATH_PoN"]

# Data files
DIR_fastq = config["DIR_fastq"]
DIR_bams = config["DIR_bams"]
DIR_recalibrated_base_scores = config["DIR_recalibrated_base_scores"]

# Fastqc reports
DIR_fastqc = config["DIR_fastqc"]

# Metrics outputs
DIR_umi_metrics = config["DIR_umi_metrics"]

# variant calling
DIR_Vardict = config["DIR_Vardict"]
DIR_Mutect = config["DIR_Mutect"]
DIR_freebayes = config["DIR_freebayes"]
DIR_results = config["DIR_results"]
DIR_ANNOVAR = config["DIR_ANNOVAR"]

# figures
DIR_metrics = config["DIR_metrics"]

all_samples = os.listdir(DIR_fastq)
samples = list({file.replace("_1.fq.gz", "").replace("_2.fq.gz", "") for file in all_samples})

cfDNA_samples = list(filter(lambda x: "cfDNA" in x, samples))
WBC_samples = list(filter(lambda x: "WBC" in x, samples))

pair1 = [sample + "_1" for sample in samples]
pair2 = [sample + "_2" for sample in samples]
all_pairs = pair1 + pair2

###################################################
# TARGET FILES FOR RULE ALL
###################################################
RUN_FastQC_merged = expand(DIR_fastqc + "/merged/{wildcard}_fastqc.html", wildcard=all_pairs)
RUN_FastQC_trimmed = expand(DIR_fastqc + "/trimmed/{wildcard}_fastqc.html", wildcard=all_pairs)

TRIM_FastQ = [expand(DIR_fastq + "/trimmed/{wildcard}_1.fq.gz", wildcard=samples), expand(DIR_fastq + "/trimmed/{wildcard}_2.fq.gz", wildcard=samples)]

make_consensus = expand(
    DIR_bams + "/{consensus_type}_fixmate/{wildcard}.bam",
    consensus_type=["SSCS2"],
    wildcard=samples,
)

filter_consensus = expand(
    DIR_bams + "/{consensus_type}_final/{wildcard}.bam.bai",
    consensus_type=["SSCS2"],
    wildcard=samples,
)

call_variants_chip = expand(
    DIR_results + "/variant_calling_chip/{variant_caller}/{consensus_type}/{wildcard}.vcf.gz",
    variant_caller=["freebayes", "Mutect2", "Vardict"],
    consensus_type=["SSCS2"],
    wildcard=samples,
    )

call_variants_somatic = expand(
    DIR_results + "/variant_calling_somatic/{variant_caller}/{consensus_type}/{wildcard}.vcf.gz",
    variant_caller=["Mutect2", "Vardict"],
    consensus_type=["SSCS2"],
    wildcard=cfDNA_samples,
    )

run_annovar_chip = expand(
    DIR_results + "/data/annovar_outputs/{variant_caller}/{consensus_type}/{wildcard}.hg38_multianno.txt",
    variant_caller=["Mutect2", "freebayes", "Vardict"],
    consensus_type=["SSCS2"],
    wildcard=samples,
)

run_annovar_somatic = expand(
    DIR_results + "/data/annovar_outputs_somatic/{variant_caller}/{consensus_type}/{wildcard}.hg38_multianno.vcf",
    variant_caller=["freebayes"],
    consensus_type=["SSCS2"],
    wildcard=cfDNA_samples,
)

vcfToTable_chip = expand(
    DIR_results + "/data/annovar_outputs/{variant_caller}/{consensus_type}/{wildcard}.hg38_multianno.txt",
    variant_caller=["Mutect2", "Vardict", "freebayes"],
    consensus_type=["SSCS1"],
    wildcard=samples,
)

variant_table_somatic = expand(
    DIR_results + "/data/variant_tables/somatic/{variant_caller}/{consensus_type}/{wildcard}.tsv",
    consensus_type=["SSCS2"],
    wildcard=cfDNA_samples,
    variant_caller=["Mutect2", "Vardict"],
)

variant_table_chip = expand(
    DIR_results + "/data/variant_tables/chip/{variant_caller}/{consensus_type}/{wildcard}.tsv",
    consensus_type=["SSCS2"],
    wildcard=samples,
    variant_caller=["Vardict", "Mutect2", "freebayes"],
)

run_insert_size = expand(
    DIR_insertsize_metrics + "/{consensus_type}/{wildcard}.txt",
    consensus_type=["SSCS2"],
    wildcard=VIP_normals_cfDNA,
)

run_depth = expand(
    DIR_depth_metrics + "/{consensus_type}/{wildcard}.txt",
    consensus_type=["SSCS2"],
    wildcard=samples,
)

RUN_mpileup = expand(
    DIR_mpileup + "/{consensus_type}/{wildcard}.mpileup",
    consensus_type=["SSCS2"],
    wildcard=samples,
)

collect_HS_metrics = expand(
    DIR_HS_metrics + "/{consensus_type}/{wildcard}.HS_metrics",
    consensus_type=["SSCS2"],
    wildcard=samples,
)

rule all:
    input:
        RUN_FastQC_merged,
        RUN_FastQC_trimmed,
        RUN_mpileup,
        run_insert_size, 
        make_consensus,
        filter_consensus,
        call_variants_chip,
        run_annovar_chip,
        call_variants_somatic,
        run_annovar_somatic,
        run_depth, 
        collect_HS_metrics,
        variant_table_somatic,
        variant_table_chip,

##### Modules #####
include: "rules/process_fastq.smk"
include: "rules/make_consensus_bams.smk"
include: "rules/filter_consensus_bams.smk"
include: "rules/run_metrics.smk"
include: "rules/variant_calling.smk"
include: "rules/variant_calling_somatic.smk"
include: "rules/annotation.smk"
include: "rules/index_bams.smk"
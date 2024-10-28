rule run_VarDict_somatic:
    input:
        cfDNA=DIR_bams + "/{consensus_type}_final/{wildcard}.bam",
        wbc=lambda wildcards: get_wbc(wildcards.wildcard, "{consensus_type}", "final")[0],
    output:
        temp(
            DIR_results
            + "/variant_calling_somatic_raw/Vardict/{consensus_type}/{wildcard}.vcf"
        ),
    conda:
        "../envs/vardict_env.yaml"
    params:
        PATH_hg38=PATH_hg38,
        PATH_bed=PATH_bed,
        THRESHOLD_VarFreq=0.005,
        sample_name_cfDNA="{wildcard}",
        sample_name_wbc=lambda wildcards: get_wbc_name(wildcards.wildcard),
    threads: 12
    shell:
        """
        VarDict \
            -G {params.PATH_hg38} \
            -b '{input.cfDNA}|{input.wbc}' \
            -f {params.THRESHOLD_VarFreq} \
            -N {params.sample_name_cfDNA} \
            -k 1 \
            -c 1 \
            -S 2 \
            -E 3 \
            -g 4 \
            -r 4 \
            --nosv \
            -th {threads} \
            {params.PATH_bed} | \
            ~/VarDictJava/build/install/VarDict/bin/testsomatic.R | \
            ~/VarDictJava/build/install/VarDict/bin/var2vcf_paired.pl \
        -N '{params.sample_name_cfDNA}|{params.sample_name_wbc}' -f {params.THRESHOLD_VarFreq} > {output}
        """

rule run_mutect2_somatic:
    input:
        cfDNA=DIR_bams + "/{consensus_type}_final/{wildcard}.bam",
        wbc=lambda wildcards: get_wbc(wildcards.wildcard, "{consensus_type}", "final")[0],
    output:
        vcf=temp(DIR_results + "/variant_calling_somatic_raw/Mutect2/{consensus_type}/{wildcard}.vcf.gz"),
        stats=DIR_results + "/variant_calling_somatic_raw/Mutect2/{consensus_type}/{wildcard}.vcf.gz.stats",
    conda:
        "../envs/chip_variantcalling.yaml"
    params:
        PATH_hg38=PATH_hg38,
        PATH_bed=PATH_bed,
        sample_name_wbc=lambda wildcards: get_wbc_name(wildcards.wildcard),
    threads: 12
    shell:
        "gatk Mutect2 \
        --reference {params.PATH_hg38} \
        --input {input.cfDNA} \
        --input {input.wbc} \
        --normal-sample {params.sample_name_wbc} \
        --output {output.vcf} \
        --force-active true \
        --initial-tumor-lod 0 \
        --tumor-lod-to-emit 0 \
        --intervals {params.PATH_bed}"

from os.path import basename
rule run_freebayes_somatic:
    input:
        cfDNA=DIR_bams + "/{consensus_type}_final/{wildcard}.bam",
        wbc=lambda wildcards: get_wbc(wildcards.wildcard, "{consensus_type}", "final")[0],
    output:
        DIR_results + "/variant_calling_somatic_raw/freebayes/{consensus_type}/{wildcard}.vcf"
    conda:
        "../envs/chip_variantcalling.yaml"
    params:
        PATH_hg38=PATH_hg38,
        PATH_bed=PATH_bed,
        cfDNA_name=lambda wildcards: basename(wildcards.wildcard),
        WBC_name=lambda wildcards: get_wbc_name(wildcards.wildcard)
    shell:
        "freebayes \
            -f {params.PATH_hg38} \
            -t {params.PATH_bed} \
            --pooled-discrete \
            --genotype-qualities \
            --haplotype-length 0 \
            --min-alternate-count 4 \
            --min-alternate-fraction 0.001 \
            {input.cfDNA} {input.wbc} |\
             vcfsamplediff somatic,germline {params.cfDNA_name} {params.WBC_name} - |\
             sed 's/somatic,germline/somatic_germline/g' > {output}"

rule zip_vcf_files_freebayes_somatic:
    input:
        DIR_results + "/variant_calling_somatic_raw/freebayes/{consensus_type}/{wildcard}.vcf",
    output:
        temp(DIR_results + "/variant_calling_somatic_raw/freebayes/{consensus_type}/{wildcard}.vcf.gz"),
    conda:
        "../envs/bcftools.yaml"
    shell:
        "bgzip -c {input} > {output}"

rule zip_vcf_files_vardict_somatic:
    input:
        DIR_results + "/variant_calling_somatic_raw/Vardict/{consensus_type}/{wildcard}.vcf",
    output:
        temp(DIR_results + "/variant_calling_somatic_raw/Vardict/{consensus_type}/{wildcard}.vcf.gz"),
    conda:
        "../envs/bcftools.yaml"
    shell:
        "bgzip -c {input} > {output}"

rule index_vcf_files_somatic:
    input:
        DIR_results + "/variant_calling_somatic_raw/{variant_caller}/{consensus_type}/{wildcard}.vcf.gz",
    output:
        temp(DIR_results + "/variant_calling_somatic_raw/{variant_caller}/{consensus_type}/{wildcard}.vcf.gz.tbi"),
    conda:
        "../envs/bcftools.yaml"
    shell:
        "tabix -p vcf {input}"


rule sort_vcf_somatic:
    input:
        vcf=DIR_results + "/variant_calling_somatic_raw/{variant_caller}/{consensus_type}/{wildcard}.vcf.gz",
        index=DIR_results + "/variant_calling_somatic_raw/{variant_caller}/{consensus_type}/{wildcard}.vcf.gz.tbi",
    output:
        temp(DIR_results + "/variant_calling_somatic_sorted/{variant_caller}/{consensus_type}/{wildcard}.vcf.gz"),
    conda:
        "../envs/bcftools.yaml"
    shell:
        "bcftools sort {input.vcf} -Oz -o {output}"


rule index_sorted_vcf_files_somatic:
    input:
        DIR_results + "/variant_calling_somatic_sorted/{variant_caller}/{consensus_type}/{wildcard}.vcf.gz",
    output:
        temp(DIR_results + "/variant_calling_somatic_sorted/{variant_caller}/{consensus_type}/{wildcard}.vcf.gz.tbi"),
    conda:
        "../envs/bcftools.yaml"
    shell:
        "tabix -p vcf {input}"

rule normalize_variants_somatic:
    input:
        vcf=DIR_results + "/variant_calling_somatic_sorted/{variant_caller}/{consensus_type}/{wildcard}.vcf.gz",
        index=DIR_results + "/variant_calling_somatic_sorted/{variant_caller}/{consensus_type}/{wildcard}.vcf.gz.tbi"
    output:
        temp(DIR_results + "/variant_calling_somatic_normalized/{variant_caller}/{consensus_type}/{wildcard}.vcf.gz"),
    params:
        PATH_hg38_dict="/groups/wyattgrp/reference/hg38/hg38.fa",
    conda:
        "../envs/bcftools.yaml"
    shell:
        "bcftools norm \
        {input.vcf} \
        -m-both \
        -f {params} \
        -o {output}"

rule decompose_blocksubstitutions_somatic:
    input:
        DIR_results + "/variant_calling_somatic_normalized/{variant_caller}/{consensus_type}/{wildcard}.vcf.gz",
    output:
        DIR_results + "/variant_calling_somatic/{variant_caller}/{consensus_type}/{wildcard}.vcf.gz",
    conda:
        "../envs/chip_variantcalling.yaml"
    shell:
        "vt decompose_blocksub {input} -o {output}"
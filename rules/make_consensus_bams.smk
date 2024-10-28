rule FastqToBam:
    input:
        R1=DIR_trimmed_fastq + "/{wildcard}_1.fq.gz",
        R2=DIR_trimmed_fastq + "/{wildcard}_2.fq.gz",
    output:
        temp(DIR_bams + "/uBAM/{wildcard}.bam"),
    params:
        sample="{wildcard}",
    run:
        shell(
            "fgbio FastqToBam -Xmx20G \
        --input {input.R1} {input.R2} \
        --output {output} \
        --read-structure 3M2S+T 3M2S+T \
        --umi-tag RX \
        --sample {params.sample} \
        --library lib1 \
        --platform illumina \
        --sort true"
        )

# Go back to the Fastq format
rule BamtoFastq:
    input:
        DIR_bams + "/uBAM/{wildcard}.bam",
    params:
        PATH_hg38=PATH_hg38,
    output:
        temp(DIR_fastq + "/consensus_fq/{wildcard}.fastq"),
    threads: 12
    run:
        shell(
            "picard SamToFastq -Xmx20G I={input} F={output} VALIDATION_STRINGENCY=SILENT INCLUDE_NON_PF_READS=true INCLUDE_NON_PRIMARY_ALIGNMENTS=true INTERLEAVE=true"
        )

# Generate mapped bam
rule mapBAM:
    input:
        DIR_fastq + "/consensus_fq/{wildcard}.fastq",
    params:
        PATH_hg38=PATH_hg38,
    output:
        temp(DIR_bams + "/mBAM_raw/{wildcard}.bam"),
    threads: 12 # keep the same for deterministic
    run:
        shell("bwa mem {params.PATH_hg38} {input} -p -Y -t {threads} > {output}")

rule MergeBamAlignment:
    input:
        mBAM=DIR_bams + "/mBAM_raw/{wildcard}.bam",
        uBAM=DIR_bams + "/uBAM/{wildcard}.bam",
    params:
        PATH_hg38=PATH_hg38,
    output:
        temp(DIR_bams + "/mBAM/{wildcard}.bam"),
    threads: 12
    run:
        shell(
            "picard MergeBamAlignment -Xmx20G \
            UNMAPPED={input.uBAM} \
            ALIGNED={input.mBAM} \
            O={output} \
            CLIP_OVERLAPPING_READS=false \
            CLIP_ADAPTERS=false \
            INCLUDE_SECONDARY_ALIGNMENTS=true \
            ALIGNER_PROPER_PAIR_FLAGS=true \
            EXPECTED_ORIENTATIONS=FR \
            MAX_INSERTIONS_OR_DELETIONS=-1 \
            VALIDATION_STRINGENCY=SILENT \
            REFERENCE_SEQUENCE={params.PATH_hg38} \
            CREATE_INDEX=true \
            SORT_ORDER=coordinate" # needs to be sorted for ABRA2
        )

rule indel_realignment:
    input:
        MAPPED_bam=DIR_bams + "/mBAM/{wildcard}.bam",
        PATH_bed=PATH_bed,
        PATH_hg38=PATH_hg38,
    params:
        min_mapping_quality=20,
        # tumor_normal_pair_input = lambda wildcards: get_relevant_bams_abra2_INPUT(wildcards.wildcard)
    output:
        DIR_bams + "/abra2/{wildcard}.bam",
    threads: 12
    run:
        shell(
            "java -Xms64G -jar abra2.jar \
        --in {input.MAPPED_bam} \
        --out {output} \
        --ref {input.PATH_hg38} \
        --threads {threads} \
        --index \
        --nosort \
        --no-edge-ci \
        --sa \
        --targets {input.PATH_bed} \
        --tmpdir tmpdir"
        )

rule fixmate:
    input:
        DIR_bams + "/abra2/{wildcard}.bam",
    output:
        temp(DIR_bams + "/fixmate/{wildcard}.bam"),
    threads: 12
    run:
        shell(
            "picard -Xmx40g FixMateInformation I={input} O={output} SORT_ORDER=coordinate VALIDATION_STRINGENCY=SILENT"
        )

rule recalibrate_bases:
    input:
        DIR_bams + "/fixmate/{wildcard}.bam",
    output:
        DIR_recalibrated_base_scores + "/{wildcard}_table",
    params:
        PATH_hg38=PATH_hg38,
        PATH_known_indels=PATH_known_indels,
        PATH_gold_std_indels=PATH_gold_std_indels,
        PATH_SNP_db=PATH_SNP_db,
    threads: 12
    run:
        shell(
            "gatk BaseRecalibrator -I {input} -R {params.PATH_hg38} --known-sites {params.PATH_known_indels} --known-sites {params.PATH_gold_std_indels} --known-sites {params.PATH_SNP_db} -O {output}"
        )


rule apply_base_scores:
    input:
        fixmate_BAM=DIR_bams + "/fixmate/{wildcard}.bam",
        base_scores=DIR_recalibrated_base_scores + "/{wildcard}_table",
    output:
        temp(DIR_bams + "/uncollapsed_BAM/{wildcard}.bam"),
    params: 
        PATH_hg38=PATH_hg38
    threads: 12
    run:
        shell(
            "gatk ApplyBQSR --reference {params.PATH_hg38} --input {input.fixmate_BAM} --output {output} --bqsr-recal-file {input.base_scores}"
        )


# Identify reads or read pairs that originate from the same source molecule based on genomic positions and UMI
rule GroupReadsByUmi:
    input:
        DIR_bams + "/uncollapsed_BAM/{wildcard}.bam",  # output of the fixmate_and_recalibrate_bases rule from the process_bams.smk file
    output:
        bam = temp(DIR_bams + "/grouped_umi_BAM/{wildcard}.bam"),
        family_size_hist = DIR_umi_metrics + "/{wildcard}.umi_metrics"
    params:
        PATH_hg38=PATH_hg38,
        min_map_quality = 20
    run:
        shell(
            "fgbio GroupReadsByUmi -Xmx20G --input={input} --output={output.bam} --strategy=paired --edits=1 --min-map-q={params.min_map_quality} --family-size-histogram={output.family_size_hist}"
        )

rule CallMolecularConsensusReads:
    input:
        DIR_bams + "/grouped_umi_BAM/{wildcard}.bam",  # output of the fixmate_and_recalibrate_bases rule from the process_bams.smk file
    output:
        DIR_bams + "/{consensus_type}_uBAM/{wildcard}.bam",
    threads: 12
    params:
        min_reads=1,
        min_input_base_quality=20,
        error_rate_post_umi=40,
        error_rate_pre_umi=45,
    run:
        shell(
            "fgbio CallMolecularConsensusReads -Xmx20G --input={input} --output={output} --threads={threads} \
            --read-name-prefix=singlex \
            --sort-order=Queryname \
            --consensus-call-overlapping-bases=true \
            --error-rate-pre-umi={params.error_rate_pre_umi} \
            --error-rate-post-umi={params.error_rate_post_umi} \
            --min-input-base-quality={params.min_input_base_quality} \
            --min-reads={params.min_reads}"
        )
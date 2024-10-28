rule BamtoFastq2:
    input:
        DIR_bams + "/{consensus_type}_uBAM/{wildcard}.bam",
    output:
        temp(DIR_fastq + "/{consensus_type}_consensus_fq_2/{wildcard}.fastq"),
    threads: 12
    run:
        shell(
            "picard SamToFastq -Xmx20G \
            I={input}\
            VALIDATION_STRINGENCY=SILENT\
            INCLUDE_NON_PRIMARY_ALIGNMENTS=true \
            INCLUDE_NON_PF_READS=true \
            INTERLEAVE=true \
            F={output}"
        )

rule mapBAM2:
    input:
        DIR_fastq + "/{consensus_type}_consensus_fq_2/{wildcard}.fastq",
    params:
        PATH_hg38=PATH_hg38,
    output:
        temp(DIR_bams + "/{consensus_type}_mBAM_raw_2/{wildcard}.bam"),
    threads: 12
    run:
        shell("bwa mem {params.PATH_hg38} {input} -p -Y -t {threads} > {output}")


rule MergeBamAlignment2:
    input:
        mBAM=DIR_bams + "/{consensus_type}_mBAM_raw_2/{wildcard}.bam",
        uBAM=DIR_bams + "/{consensus_type}_uBAM/{wildcard}.bam",
    params:
        PATH_hg38=PATH_hg38,
    output:
        temp(DIR_bams + "/{consensus_type}_mBAM/{wildcard}.bam"),
    threads: 12
    run:
        shell(
            "picard MergeBamAlignment -Xmx20G UNMAPPED={input.uBAM} ALIGNED={input.mBAM} O={output} R={params.PATH_hg38} \
                CLIP_OVERLAPPING_READS=false \
                CLIP_ADAPTERS=false \
                ALIGNER_PROPER_PAIR_FLAGS=true \
                INCLUDE_SECONDARY_ALIGNMENTS=true \
                PRIMARY_ALIGNMENT_STRATEGY=MostDistant \
                EXPECTED_ORIENTATIONS=FR \
                MAX_INSERTIONS_OR_DELETIONS=-1 \
                VALIDATION_STRINGENCY=SILENT \
                CREATE_INDEX=true \
                SORT_ORDER=unsorted"
        )

rule subset_to_proper_pairs:
    input:
        DIR_bams + "/{consensus_type}_mBAM/{wildcard}.bam"
    output:
        temp(DIR_bams + "/{consensus_type}_proper_pair/{wildcard}.bam"),
    threads: 12
    run:
        shell(
            "sambamba view {input} -F 'proper_pair' -t 12 -f bam -l 0 -o {output}"
        )

# abra2 requires sorted and indexed bams
rule sort_subsetted_bams:
    input:
        DIR_bams + "/{consensus_type}_proper_pair/{wildcard}.bam",
    output:
        temp(DIR_bams + "/{consensus_type}_proper_pair_sorted/{wildcard}.bam"),
    threads: 12
    run:
        shell(
            "samtools sort -o {output} {input}"
        )

rule index_sorted_subsetted_bams:
    input:
        DIR_bams + "/{consensus_type}_proper_pair_sorted/{wildcard}.bam",
    output:
        temp(DIR_bams + "/{consensus_type}_proper_pair_sorted/{wildcard}.bam.bai"),
    threads: 12
    run:
        shell(
            "samtools index {input}"
        )

rule indel_realignment2:
    input:
        MAPPED_bam=DIR_bams + "/{consensus_type}_proper_pair_sorted/{wildcard}.bam",
        MAPPED_bam_index=DIR_bams + "/{consensus_type}_proper_pair_sorted/{wildcard}.bam.bai",
        PATH_bed=PATH_bed,
        PATH_hg38=PATH_hg38,
    output:
        temp(DIR_bams + "/{consensus_type}_abra2/{wildcard}.bam"),
    threads: 12
    run:
        shell(
            "java -Xms64G -jar abra2.jar \
        --in {input.MAPPED_bam} \
        --out {output} \
        --ref {input.PATH_hg38} \
        --threads {threads} \
        --mad 5000 \
        --no-edge-ci \
        --targets {input.PATH_bed} \
        --tmpdir tmpdir"
        )

rule fixmate2:
    input:
        DIR_bams + "/{consensus_type}_abra2/{wildcard}.bam",
    output:
        DIR_bams + "/{consensus_type}_fixmate/{wildcard}.bam"
    threads: 12
    run:
        shell(
            "picard -Xmx40g FixMateInformation I={input} O={output} SORT_ORDER=coordinate VALIDATION_STRINGENCY=SILENT"
        )

rule FilterConsensusReads_SSCS:
    input:
        DIR_bams + "/{consensus_type}_fixmate/{wildcard}.bam",
    params:
        PATH_hg38=PATH_hg38,
        PATH_bed=PATH_bed,
        sample="{wildcard}",
        min_reads=2,
        min_base_quality=30,
        max_read_error_rate=0.025,
        max_no_call_fraction=0.15,
    output:
        temp(DIR_bams + "/{consensus_type}_final_no_rg/{wildcard}.bam"),
    threads: 12
    run:
        shell(
            "samtools sort -n {input} | fgbio FilterConsensusReads -Xmx20G \
            --input=/dev/stdin \
            --output={output} \
            --ref={params.PATH_hg38} \
            --min-reads={params.min_reads} \
            --max-read-error-rate={params.max_read_error_rate} \
            --max-no-call-fraction={params.max_no_call_fraction} \
            --min-base-quality={params.min_base_quality} \
            --reverse-per-base-tags=true \
            --sort-order=coordinate"
        )

rule add_rg_for_freebayes_somatic: 
    input:
        DIR_bams + "/{consensus_type}_final_no_rg/{wildcard}.bam",
    output:
        bam=DIR_bams + "/{consensus_type}_final/{wildcard}.bam",
    params:
        rg_id=lambda wildcards: "cfDNA" if "cfDNA" in wildcards.wildcard else "WBC",
        rg_lib=lambda wildcards: wildcards.wildcard + "_lib",
        rg_pl="ILLUMINA",
        rg_pu="{wildcard}_unit1",
        rg_sm="{wildcard}"
    shell:
        """
        picard -Xmx40g AddOrReplaceReadGroups \
            I={input} \
            O={output.bam} \
            RGID={params.rg_id} \
            RGLB={params.rg_lib} \
            RGPL={params.rg_pl} \
            RGPU={params.rg_pu} \
            RGSM={params.rg_sm} && samtools index {output.bam}
        """

rule index_rg_bams: 
    input:
        DIR_bams + "/{consensus_type}_final/{wildcard}.bam"
    output:
        DIR_bams + "/{consensus_type}_final/{wildcard}.bam.bai"
    shell: 
        "samtools index {input}"
rule index_bams:
    input:
        DIR_results + "/data/bam/{consensus_type}_filtered/{wildcard}.bam"
    output:
        DIR_results + "/data/bam/{consensus_type}_filtered/{wildcard}.bam.bai"
    shell:
        "samtools index {input}"
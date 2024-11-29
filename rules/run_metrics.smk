rule run_depth: 
	input: 
		BAM = DIR_bams + "/{consensus_type}_final/{wildcard}.bam",
		PATH_bed = PATH_bed 
	output:
		DIR_metrics + "depth_metrics/{consensus_type}/{wildcard}.txt"
	threads: 12
	shell:
		"samtools depth -b {input.PATH_bed} {input.BAM} > {output}"

rule run_mpileup: 
	input: 
		BAM = DIR_bams + "/{consensus_type}_final/{wildcard}.bam", 
		BAM_index = DIR_bams + "/{consensus_type}_final/{wildcard}.bam.bai",
		PATH_hg38 = PATH_hg38, 
		PATH_bed = PATH_bed,
	output:
		DIR_mpileup + "/{consensus_type}/{wildcard}.mpileup"
	threads: 12
	shell:
		"samtools mpileup -A -f {input.PATH_hg38} --no-BAQ --positions {input.PATH_bed} --min-MQ 0 --min-BQ 0 {input.BAM} -o {output}"

rule run_insert_size: 
	input: 
		DIR_bams + "/{consensus_type}_final/{wildcard}.bam"
	output:
		metrics = DIR_metrics + "/insert_size_metrics/{consensus_type}/{wildcard}.txt",
		figures = DIR_metrics + "/insert_size_figures/{consensus_type}/{wildcard}.png"
	threads: 12
	shell:
		"picard CollectInsertSizeMetrics \
			I={input} \
			O={output.metrics} \
			H={output.figures} \
			W=2000 \
			M=0.5"

rule run_read_counts: 
	input: 
		DIR_merged_fastq + "/{wildcard}_1.fq.gz"
	params:
		sample_name = "{wildcard}"
	output:
		DIR_metrics + "/fastq_read_counts/merged/{wildcard}.txt",
	shell:
		"paste <(echo {params}) <(echo $(gunzip -c {input}|wc -l)/4|bc) > {output}"

# Using the baits interval file for both the target regions and the baits input.
rule hs_metrics: 
	input: 
		bam = DIR_bams + "/{consensus_type}_final/{wildcard}.bam",
		PATH_baits = PATH_baits,
		PATH_bed_intervals = PATH_bed_intervals,
		PATH_hg38 = PATH_hg38
	output:
		DIR_metrics + "/hs_metrics/{consensus_type}/{wildcard}.HS_metrics",
	shell:
		"picard CollectHsMetrics \
    		I={input.bam} \
    		O={output} \
    		R={input.PATH_hg38} \
    		BAIT_INTERVALS={input.PATH_baits} \
    		TARGET_INTERVALS={input.PATH_baits}"
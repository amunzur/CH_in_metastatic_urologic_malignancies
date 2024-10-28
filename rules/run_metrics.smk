rule run_depth: 
	input: 
		BAM = DIR_bams + "/{consensus_type}_final/{wildcard}.bam",
		PATH_bed = PATH_bed 
	output:
		DIR_depth_metrics + "/{consensus_type}/{wildcard}.txt"
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
		metrics = DIR_insertsize_metrics + "{consensus_type}/{wildcard}.txt",
		figures = DIR_insertsize_figures + "{consensus_type}/{wildcard}.pdf"
	threads: 12
	shell:
		"picard CollectInsertSizeMetrics \
			I={input} \
			O={output.metrics} \
			H={output.figures} \
			M=0.5"

rule run_read_counts: 
	input: 
		DIR_merged_fastq + "/{wildcard}_1.fq.gz"
	params:
		sample_name = "{wildcard}"
	output:
		DIR_readcounts_metrics + "/merged/{wildcard}.txt",
	shell:
		"paste <(echo {params}) <(echo $(gunzip -c {input}|wc -l)/4|bc) > {output}"

# Using the baits interval file for both the target regions and the baits input.
rule hs_metrics: 
	input: 
		SSCS_bam = DIR_bams + "/{consensus_type}_final/{wildcard}.bam",
		PATH_baits = PATH_baits,
		PATH_bed_intervals = PATH_bed_intervals,
		PATH_hg38 = PATH_hg38
	output:
		DIR_HS_metrics + "/{consensus_type}/{wildcard}.HS_metrics",
	shell:
		"picard CollectHsMetrics \
    		I={input.SSCS_bam} \
    		O={output} \
    		R={input.PATH_hg38} \
    		BAIT_INTERVALS={input.PATH_baits} \
    		TARGET_INTERVALS={input.PATH_baits}"

# rule run_tnvstats:
# 	input:
# 		BAM_tumor_path=DIR_bams + "/SSCS1_clipped/{wildcard}.bam",
# 		BAM_normal_path=
# 		BAM_normal_name=
# 				BAM_index=DIR_bams + "/SSCS1_clipped/{wildcard}.bam.bai",

# 	params:
# 	output:
# 	conda:
#         "../envs/pysamstat_v2.yaml"
# 	shell:
# 		'script_dir="/groups/wyattgrp/users/amunzur/gene_panel_pipeline/scripts";
# 		mkdir -p /groups/wyattgrp/users/amunzur/chip_project/tnvstats/kidney_samples; # output dir to place all tnvstats
# 		configfile="/groups/wyattgrp/users/amunzur/gene_panel_pipeline/scripts/config.txt"; # conda envs and stuff

# 		source ${configfile};
# 		source ${conda_profile_path};
# 		conda activate ${conda_pysam_env};

# 		while read name_tumor path_tumor name_normal path_normal
# 		do
# 		    echo $name_tumor
# 		    echo $name_normal
# 		    outputdir="/groups/wyattgrp/users/amunzur/chip_project/tnvstats/kidney_samples/${name_tumor}";

# 		    # tbam=$(readlink -ve ${path_tumor})
# 		    # nbam=$(readlink -ve ${path_normal})
# 		    printf "bash ${script_dir}/make_tnvstat_file.bash ${path_tumor} ${path_normal} ${name_tumor} ${name_normal} ${outputdir} ${configfile}\n";
# 		    sbatch ${script_dir}/make_tnvstat_file.bash ${path_tumor} ${path_normal} ${name_tumor} ${name_normal} ${outputdir} ${configfile};
# 		done < /groups/wyattgrp/users/amunzur/pipeline/results/metrics/tnvstats/kidney_samples/bam_subset.txt'
rule run_fastqc_merged:
    input:
        DIR_fastq + "/merged/{wildcard}.fq.gz",
    output:
        output_zip=DIR_fastqc + "/merged/{wildcard}_fastqc.zip",
        output_html=DIR_fastqc + "/merged/{wildcard}_fastqc.html",
    threads: 5
    shell:
        "fastqc {input} --outdir=`dirname {output.output_zip}`"

rule trim_fastq:
    input:
        R1=DIR_fastq + "/merged/{wildcard}_1.fq.gz",
        R2=DIR_fastq + "/merged/{wildcard}_2.fq.gz",
    output:
        R1=temp(DIR_trimmed_fastq + "/{wildcard}_1.fq.gz"),
        R2=temp(DIR_trimmed_fastq + "/{wildcard}_2.fq.gz"),
        html_report="results/reports/fastp/{wildcard}.html",
        json_report="results/reports/fastp/{wildcard}.json",
    threads: 12
    params:
        minimum_read_length=32,
    run:
        shell(
            "fastp -i {input.R1} -I {input.R2} -o {output.R1} -O {output.R2} \
        --dont_eval_duplication \
        --cut_tail \
        --length_required {params.minimum_read_length} \
        --html {output.html_report} \
        --json {output.json_report} \
        --adapter_sequence AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
        --adapter_sequence_r2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
        --thread {threads}"
        )

rule run_fastqc_trimmed:
    input:
        R1=DIR_fastq + "/trimmed/{wildcard}_1.fq.gz",
        R2=DIR_fastq + "/trimmed/{wildcard}_2.fq.gz",
    output:
        output_zip_R1=DIR_fastqc + "/trimmed/{wildcard}_1_fastqc.zip",
        output_html_R1=DIR_fastqc + "/trimmed/{wildcard}_1_fastqc.html",
        output_zip_R2=DIR_fastqc + "/trimmed/{wildcard}_2_fastqc.zip",
        output_html_R2=DIR_fastqc + "/trimmed/{wildcard}_2_fastqc.html",
    threads: 5
    shell:
        """
        fastqc {input.R1} --outdir=`dirname {output.output_zip_R1}`
        fastqc {input.R2} --outdir=`dirname {output.output_zip_R2}` 
        """
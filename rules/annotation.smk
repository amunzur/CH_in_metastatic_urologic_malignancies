rule run_ANNOVAR_chip:
    input:
        DIR_results + "/variant_calling_chip/{variant_caller}/{consensus_type}/{wildcard}.vcf.gz",
    output:
        DIR_ANNOVAR + "/chip_variants/{variant_caller}/{consensus_type}/{wildcard}.hg38_multianno.txt",
        DIR_ANNOVAR + "/chip_variants/{variant_caller}/{consensus_type}/{wildcard}.hg38_multianno.vcf",
    params:
        placeholder = DIR_results + "/data/annovar_outputs/{variant_caller}/{consensus_type}/{wildcard}",
    shell:
        "perltable_annovar.pl {input} humandb \
        -vcfinput \
        -buildver hg38 \
        -out {params.placeholder} \
        -remove \
        -protocol refGene,cosmic97_coding,avsnp150,clinvar_20221231,gnomad40_exome \
        -operation g,f,f,f,f \
        -nastring ."

rule run_ANNOVAR_ctDNA:
    input:
        DIR_results + "/variant_calling_somatic/{variant_caller}/{consensus_type}/{wildcard}.vcf.gz",
    output:
        DIR_ANNOVAR + "/ctDNA_variants/{variant_caller}/{consensus_type}/{wildcard}.hg38_multianno.vcf",
        DIR_ANNOVAR + "/ctDNA_variants/{variant_caller}/{consensus_type}/{wildcard}.hg38_multianno.txt",
    params:
        DIR_ANNOVAR + "/ctDNA_variants/{variant_caller}/{consensus_type}/{wildcard}",
    shell:
        "perl table_annovar.pl {input} humandb \
        -vcfinput \
        -buildver hg38 \
        -out {params} \
        -remove \
        -protocol refGene,cosmic97_coding,avsnp150,clinvar_20221231,gnomad40_exome \
        -operation g,f,f,f,f \
        -nastring ."

# ctDNA mutations, ildcard is cfDNA.
rule vcfToTable_freebayes_somatic:
    input:
        DIR_ANNOVAR + "/ctDNA_variants/freebayes/{consensus_type}/{wildcard}.hg38_multianno.vcf",
    output:
        DIR_results + "/data/variant_tables/somatic/freebayes/{consensus_type}/{wildcard}.tsv",
    params:
        DIR_ANNOVAR + "/freebayes/{consensus_type}/{wildcard}",
    shell:
        """
        gatk VariantsToTable \
            -V {input} \
            -F CHROM -F POS -F REF -F ALT -F TYPE -F SAF -F SRF -F SAR -F SRR -F SAP -F somatic_germline \
            -F Func.refGene -F Gene.refGene -F ExonicFunc.refGene -F AAChange.refGene -F cosmic97_coding -F avsnp150 -F gnomad40_exome_AF -F CLNALLELEID -F CLNSIG \
            -O {output}
        """

rule vcfToTable_Vardict:
    input:
        DIR_ANNOVAR + "/ctDNA_variants/Vardict/{consensus_type}/{wildcard}.hg38_multianno.vcf",
    output:
        DIR_results + "/data/variant_tables/ctDNA_variants/Vardict/{consensus_type}/{wildcard}.tsv",
    shell:
        """
        gatk VariantsToTable \
        -V {input} \
        -F CHROM -F POS -F REF -F ALT -F TYPE -F FILTER -F STATUS -GF DP -GF VD -GF AF -GF ALD -GF RD -GF SBF -GF ODDRATIO \
        -F Func.refGene -F Gene.refGene -F ExonicFunc.refGene -F AAChange.refGene -F cosmic97_coding -F avsnp150 -F gnomad40_exome_AF -F CLNALLELEID -F CLNSIG \
        -O {output}
        """

rule vcfToTable_Mutect2_somatic:
    input:
        DIR_ANNOVAR + "/ctDNA_variants/Mutect2/{consensus_type}/{wildcard}.hg38_multianno.vcf",
    output:
        DIR_results + "/data/variant_tables/ctDNA_variants/Mutect2/{consensus_type}/{wildcard}.tsv",
    shell:
        """
        gatk VariantsToTable \
            -V {input} \
            -F CHROM -F POS -F REF -F ALT -F TYPE -F EVENTLENGTH -GF SB \
            -F Func.refGene -F Gene.refGene -F ExonicFunc.refGene -F AAChange.refGene -F AF -F cosmic97_coding -F avsnp150 -F gnomad40_exome_AF -F CLNALLELEID -F CLNSIG \
            -O {output}
        """

#CHIP
# wildcard is samples, both cfDNA and gDNA processed individually
rule vcfToTable_freebayes_chip:
    input:
        DIR_ANNOVAR + "/chip_variants/freebayes/{consensus_type}/{wildcard}.hg38_multianno.vcf",
    output:
        DIR_results + "/data/variant_tables/chip_variants/freebayes/{consensus_type}/{wildcard}.tsv",
    params:
        DIR_ANNOVAR + "/freebayes/{consensus_type}/{wildcard}",
    shell:
        """
        gatk VariantsToTable \
            -V {input} \
            -F CHROM -F POS -F REF -F ALT -F TYPE -F SAF -F SRF -F SAR -F SRR -F SAP \
            -F Func.refGene -F Gene.refGene -F ExonicFunc.refGene -F AAChange.refGene -F cosmic97_coding -F avsnp150 -F gnomad40_exome_AF -F CLNALLELEID -F CLNSIG \
            -O {output}
        """

rule vcfToTable_Vardict_chip:
    input:
        DIR_ANNOVAR + "/chip_variants/Vardict/{consensus_type}/{wildcard}.hg38_multianno.vcf",
    output:
        DIR_results + "/data/variant_tables/chip_variants/Vardict/{consensus_type}/{wildcard}.tsv",
    shell:
        """
        gatk VariantsToTable \
            -V {input} \
            -F CHROM -F POS -F REF -F ALT -F TYPE -F FILTER -F STATUS -GF DP -GF VD -GF AF -GF ALD -GF RD -GF BIAS -GF ODDRATIO\
            -F Func.refGene -F Gene.refGene -F ExonicFunc.refGene -F AAChange.refGene -F cosmic97_coding -F avsnp150 -F gnomad40_exome_AF -F CLNALLELEID -F CLNSIG \
            -O {output}
        """

# wildcard is cfDNA
rule vcfToTable_Mutect2_chip:
    input:
        DIR_ANNOVAR + "/chip_variants/Mutect2/{consensus_type}/{wildcard}.hg38_multianno.vcf",
    output:
        DIR_results + "/data/variant_tables/chip_variants/Mutect2/{consensus_type}/{wildcard}.tsv",
    shell:
        """
        gatk VariantsToTable \
            -V {input} \
            -F CHROM -F POS -F REF -F ALT -F TYPE -F EVENTLENGTH -GF SB \
            -F Func.refGene -F Gene.refGene -F ExonicFunc.refGene -F AAChange.refGene -F AF -F cosmic97_coding -F avsnp150 -F gnomad40_exome_AF -F CLNALLELEID -F CLNSIG \
            -O {output}
        """
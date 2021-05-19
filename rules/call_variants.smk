rule gatk_HaplotypeCaller_ERC_GVCF:
    input:    
        bam="reads/recalibrated/{sample}.dedup.recal.bam",
        bai="reads/recalibrated/{sample}.dedup.recal.bai"
#        cram="reads/recalibrated/{sample}.dedup.recal.cram",
#        crai="reads/recalibrated/{sample}.dedup.recal.cram.crai"
    output:
        gvcf="variant_calling/{sample}.g.vcf.gz"
    conda:
       "../envs/gatk.yaml"
    params:
        custom=java_params(tmp_dir=config.get("tmp_dir"), multiply_by=5),
        intervals = config.get("interval_list"),
        genome=resolve_single_filepath(*references_abs_path(), config.get("genome_fasta"))
    log:
        "logs/gatk/HaplotypeCaller/{sample}.genotype_info.log"
    benchmark:
        "benchmarks/gatk/HaplotypeCaller/{sample}.txt"
    threads: 2
    shell:
        "gatk HaplotypeCaller --java-options {params.custom} "
        "-R {params.genome} "
        "-I {input.bam} "
        "-O {output.gvcf} "
        "-ERC GVCF "
        "-L {params.intervals} "
        "-ip 200 "
        "-G StandardAnnotation "
        # "--use-new-qual-calculator "
        ">& {log}"
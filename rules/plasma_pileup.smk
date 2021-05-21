rule gatk_pileup:
    input:    
        bam="reads/recalibrated/{sample}.dedup.recal.bam",
        bai="reads/recalibrated/{sample}.dedup.recal.bai"
#        cram="reads/recalibrated/{sample}.dedup.recal.cram",
#        crai="reads/recalibrated/{sample}.dedup.recal.cram.crai"
    output:
        gvcf="variant_calling/{sample}.pileup.txt"
    conda:
       "../envs/gatk.yaml"
    params:
        custom=java_params(tmp_dir=config.get("tmp_dir"), multiply_by=5),
        intervals = config.get("interval_list"),
        genome=resolve_single_filepath(*references_abs_path(), config.get("genome_fasta"))
    log:
        "logs/gatk/pileup/{sample}.genotype_info.log"
    benchmark:
        "benchmarks/gatk/pileup/{sample}.txt"
    threads: 2
    shell:
        "gatk Pileup --java-options {params.custom} "
        "-R {params.genome} "
        "-I {input.bam} "
        "-O {output.gvcf} "
        "-L {params.intervals} "
        "-ip 200 "
        "--read-filter FragmentLengthReadFilter"
        "--max-fragment-length 1000000 "
        "--min-fragment-length 0 "
        "--minimum-mapping-quality 40" #default 10

        ">& {log}"
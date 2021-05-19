rule vcf_SelectVariants:
    input:
        vcf="variant_calling/all.snp_recalibrated.indel_recalibrated.vcf.gz",
        #lambda wildcards: get_sample_by_set(wildcards, sets)
    output:
        vcf="variant_calling/SelectVariants/{set}.selected.vcf"
    params:
        custom=java_params(tmp_dir=config.get("tmp_dir"), multiply_by=5),
        genome=resolve_single_filepath(*references_abs_path(), config.get("genome_fasta")),
        arguments=_multi_flag(config.get("rules").get("gatk_SelectVariants").get("arguments")),
#        samples_files=_get_samples_set(config.get("rules").get("gatk_SelectVariants").get("samples_files"))
        samples=lambda wildcards: get_sample_by_set(wildcards, sets)
    log:
        "logs/gatk/SelectVariants/{set}.SelectVariants.log"
    benchmark:
        "benchmarks/gatk/SelectVariants/{set}.SelectVariants.txt"
    conda:
       "../envs/gatk.yaml"

    shell:
        "gatk SelectVariants "
        "--java-options {params.custom} "
        "-R {params.genome} "
        "-V {input.vcf} "
        "{params.samples} "
        "-O {output.vcf} "
        "{params.arguments} "
#        "{params.samples_files} "
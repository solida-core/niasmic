

rule filter_mutect:
    input:
        vcf="data/results/{patient}_somatic.vcf.gz",
        cont_tab="data/filters/{patient}_contamination.table"
    output:
        "data/results/{patient}_somatic_filtered.vcf.gz"
    params:
        custom=java_params(tmp_dir=config.get("tmp_dir"), multiply_by=5),
        genome=resolve_single_filepath(*references_abs_path(), config.get("genome_fasta")),
        intervals=config.get("interval_list"),
    log:
        "logs/gatk/Mutect2/{patient}.filter_info.log"
    conda:
       "../envs/gatk.yaml"
    threads: conservative_cpu_count(reserve_cores=2, max_cores=99)
    shell:
        "gatk FilterMutectCalls "
        "--java-options {params.custom} "
        "-V {input.vcf} "
        "--tumor-segmentation "
        "--contamination-table {input.cont_tab} "
        "--ob-priors "
        "-O {output} "
        "-R {params.genome} "
        "-L {params.intervals} "
        "--native-pair-hmm-threads {threads} "
        ">& {log} "




rule gatk_SelectVariants:
    input:
        vcf="data/results/{patient}_somatic_filtered.vcf.gz"
    output:
        vcf="data/results/{patient}_somatic_filtered_selected.vcf.gz"
    params:
        custom=java_params(tmp_dir=config.get("tmp_dir"), multiply_by=5),
        genome=resolve_single_filepath(*references_abs_path(), config.get("genome_fasta")),
        arguments=_multi_flag(config.get("rules").get("gatk_SelectVariants").get("arguments"))
    log:
        "logs/gatk/SelectVariants/{patient}.SelectVariants.log"
    conda:
       "../envs/gatk.yaml"

    benchmark:
        "benchmarks/gatk/SelectVariants/{patient}.SelectVariants.txt"

    threads: conservative_cpu_count(reserve_cores=2, max_cores=99)
    shell:
        "gatk SelectVariants "
        "--java-options {params.custom} "
        "-R {params.genome} "
        "-V {input.vcf} "
        "-O {output.vcf} "
        "{params.arguments} "
        ">& {log} "


rule gatk_Funcotator:
    input:
        vcf="data/results/{patient}_somatic_filtered_selected.vcf.gz"
    output:
        vcf="data/results/{patient}_funcotated.maf"
    params:
        custom=java_params(tmp_dir=config.get("tmp_dir"), multiply_by=5),
        genome=resolve_single_filepath(*references_abs_path(), config.get("genome_fasta")),
        intervals=config.get("interval_list"),
        resource=config.get("funcotator")
    log:
        "logs/gatk/Funcotator/{patient}.funcotator.log"
    conda:
       "../envs/gatk.yaml"
    threads: conservative_cpu_count(reserve_cores=2, max_cores=99)
    shell:
        "gatk Funcotator "
        "--java-options {params.custom} "
        "-R {params.genome} "
        "-L {params.intervals}"
        "-V {input.vcf} "
        "-O {output.vcf} "
        "--data-sources-path {params.resource}"
        "--output-file-format MAF "
        "--ref-version hg19"
        ">& {log} "


#
# gatk VariantEval
# --java-options '-Xms2g -Xmx10g -XX:ParallelGCThreads=5 '
# --eval family_10_santoro.vcf
# -O family_10_santoro_eval.vcf
# -D dbSNP146_chr.vcf \
#    -R /ELS/els9/users/biosciences/references/ucsc/hg19/ucsc.hg19.fasta \
#    -L /ELS/els9/users/biosciences/references/ucsc/hg19/intervals/S04380110_Covered_Target.bed
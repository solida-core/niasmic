

rule filter_mutect_1:
    input:
        vcf="data/results/{patient}_somatic.vcf.gz",
        cont_tab="data/filters/{patient}_contamination.table"
    output:
        "data/results/{patient}_somatic_oncefiltered.vcf.gz"
    params:
        custom=java_params(tmp_dir=config.get("tmp_dir"), multiply_by=5),
    log:
        "logs/gatk/Mutect2/{patient}.filter_1_info.log"
    conda:
       "../envs/gatk.yaml"
    threads: conservative_cpu_count(reserve_cores=2, max_cores=99)
    shell:
        "gatk FilterMutectCalls "
        "--java-options {params.custom} "
        "-V {input.vcf} "
        "--contamination-table {input.cont_tab} "
        "-O {output} "
        "--native-pair-hmm-threads {threads} "
        ">& {log} "



rule filter_mutect_2:
    input:
        vcf="data/results/{patient}_somatic_oncefiltered.vcf.gz",
        artifact="data/filters/{patient}_tumor_artifact.pre_adapter_detail_metrics"
    output:
        "data/results/{patient}_somatic_twicefiltered.vcf.gz"
    params:
        custom=java_params(tmp_dir=config.get("tmp_dir"), multiply_by=5),
    log:
        "logs/gatk/Mutect2/{patient}.filter_2_info.log"
    conda:
       "../envs/gatk.yaml"
    threads: conservative_cpu_count(reserve_cores=2, max_cores=99)
    shell:
        "gatk FilterByOrientationBias "
        "--java-options {params.custom} "
        "--artifact-modes 'G/T' "
        "--artifact-modes 'C/T' "
        "-V {input.vcf} "
        "-P {input.artifact} "
        "-O {output} "
        ">& {log} "




rule gatk_SelectVariants:
    input:
        vcf="data/results/{patient}_somatic_twicefiltered.vcf.gz"
    output:
        vcf="data/results/{patient}_somatic_twicefiltered_selected.vcf.gz"
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



rule pileup_summaries_tumoral:
    input:
        lambda wildcards: get_sample_by_famid(wildcards, patients)
    output:
        "data/filters/{patient}_getpileupsummaries.table"
    params:
        custom=java_params(tmp_dir=config.get("tmp_dir"), multiply_by=5),
        intervals=config.get("interval_list"),
        exac=config.get("exac"),
        java_params=config.get("rules").get("mutect").get("java")
    log:
        "logs/gatk/Mutect2/{patient}_pileupsummaries_T.log"
    conda:
       "../envs/gatk.yaml"
    threads: conservative_cpu_count(reserve_cores=2, max_cores=99)
    shell:
        "gatk GetPileupSummaries "
        "--java-options {params.custom} "
        "-I {input[0]} " # corresponding to tbam
        "-V {params.exac} "
        "-L {params.intervals} "
        "-O {output} "
        ">& {log} "



rule pileup_summaries_normal:
    input:
        lambda wildcards: get_sample_by_famid(wildcards, patients)
    output:
        "data/filters/{patient}_normal_pileups.table"
    params:
        custom=java_params(tmp_dir=config.get("tmp_dir"), multiply_by=5),
        intervals=config.get("interval_list"),
        exac=config.get("exac")
    log:
        "logs/gatk/Mutect2/{patient}_pileupsummaries_C.log"
    conda:
       "../envs/gatk.yaml"
    threads: conservative_cpu_count(reserve_cores=2, max_cores=99)
    shell:
        "gatk GetPileupSummaries "
        "--java-options {params.custom} "
        "-I {input[1]} "
        "-V {params.exac} "
        "-L {params.intervals} "
        "-O {output} "
        ">& {log} "


rule calculate_contamination:
    input:
        tab_t="data/filters/{patient}_getpileupsummaries.table",
        tab_c="data/filters/{patient}_normal_pileups.table"
    output:
        "data/filters/{patient}_contamination.table"
    params:
        custom=java_params(tmp_dir=config.get("tmp_dir"), multiply_by=5),
    log:
        "logs/gatk/Mutect2/{patient}_calculatecontamination.log"
    conda:
       "../envs/gatk.yaml"
    threads: conservative_cpu_count(reserve_cores=2, max_cores=99)
    shell:
        "gatk CalculateContamination "
        "--java-options {params.custom} "
        "-I {input.tab_t} "
        "-O {output} "
        "-matched {input.tab_c} "
        ">& {log} "



rule calculate_seq_artifacts:
    input:
        lambda wildcards: get_sample_by_famid(wildcards, patients)
    output:
        "data/filters/{patient}_tumor_artifact.pre_adapter_detail_metrics"
    params:
        custom=java_params(tmp_dir=config.get("tmp_dir"), multiply_by=5),
        genome=resolve_single_filepath(*references_abs_path(), config.get("genome_fasta")),
        out=("data/filters/{patient}_tumor_artifact.txt").split(".")[0]
    log:
        "logs/gatk/Mutect2/{patient}_seq_artifacts_metrics.log"
    conda:
       "../envs/gatk.yaml"
    threads: conservative_cpu_count(reserve_cores=2, max_cores=99)
    shell:
        "gatk CollectSequencingArtifactMetrics "
        "--java-options {params.custom} "
        "-I {input[0]} " # corresponding to tbam
        "-O {params.out} "
        "-R {params.genome} "
        ">& {log} "





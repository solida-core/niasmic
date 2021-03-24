



rule picard_gc_bias:
    input:
        "reads/recalibrated/{sample}.dedup.recal.bam"
    output:
        chart="qc/picard/{sample}_gc_bias_metrics.pdf",
        summary="qc/picard/{sample}_summary_metrics.txt",
        out="qc/picard/{sample}_gc_bias_metrics.txt"
    params:
        custom=java_params(tmp_dir=config.get("tmp_dir"), multiply_by=5),
        genome=resolve_single_filepath(*references_abs_path(), config.get("genome_fasta")),
        param=config.get("rules").get("picard_gc").get("params"),
        tmp_dir=config.get("tmp_dir")
    log:
        "logs/picard/CollectGcBiasMetrics/{sample}.gcbias.log"
    conda:
       "../envs/picard.yaml"
    threads: conservative_cpu_count(reserve_cores=2, max_cores=99)
    shell:
        "picard "
        "CollectGcBiasMetrics "
        "{params.custom} "
        "I={input} "
        "R={params.genome} "
        "CHART={output.chart} "
        "S={output.summary} "
        "O={output.out} "
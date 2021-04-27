

rule picard_duplicate_estimation:
    input:
        bam="reads/merged/{sample}.bam",
        bai="reads/merged/{sample}.bam.bai"
    output:
        bam="qc/picard/UMIstats/bam/{sample}.bam",
        metrics="qc/picard/UMIstats/{sample}.duplicate_metrics.txt",
        umi="qc/picard/UMIstats/{sample}.umi_metrics.txt"
    params:
        custom=java_params(tmp_dir=config.get("tmp_dir"),multiply_by=5),
        genome=resolve_single_filepath(*references_abs_path(),config.get("genome_fasta")),
    shell:
        "picard UmiAwareMarkDuplicatesWithMateCigar "
        "-I={input.bam} "
        "-O={output.bam} "
        "-M={output.metrics} "
        "-UMI_METRICS={output.umi} "


rule picard_HsMetrics:
    input:
        bam="reads/recalibrated/{sample}.dedup.recal.bam",
#        probes="references/{sample}_probes_header",
#        hsTarget="references/{sample}_hsTarget_header"
    output:
        "qc/picard/hs/{sample}.dedup.recal.hs.txt"
    conda:
        "../envs/picard.yaml"
    params:
        custom=java_params(tmp_dir=config.get("tmp_dir"), multiply_by=5),
        probes=config.get("bait"),
        target=config.get("target")
    benchmark:
        "benchmarks/picard/HsMetrics/{sample}.txt"
    shell:
        "picard {params.custom} CollectHsMetrics "
        "INPUT={input.bam} OUTPUT={output} "
        # "BAIT_INTERVALS={input.probes} TARGET_INTERVALS={input.hsTarget} "
        "BAIT_INTERVALS={params.probes} TARGET_INTERVALS={params.target} "
        "CLIP_OVERLAPPING_READS=true MINIMUM_MAPPING_QUALITY=-1 MINIMUM_BASE_QUALITY=-1 "


rule picard_InsertSizeMetrics:
   input:
      bam="reads/recalibrated/{sample}.dedup.recal.bam"
   output:
       metrics="qc/picard/hs/{sample}.dedup.recal.is.txt",
       histogram="qc/picard/hs/{sample}.dedup.recal.is.pdf"
   conda:
       "../envs/picard.yaml"
   params:
        custom=java_params(tmp_dir=config.get("tmp_dir"), multiply_by=5),
   benchmark:
       "benchmarks/picard/IsMetrics/{sample}.txt"
   shell:
       "picard {params.custom} CollectInsertSizeMetrics "
       "INPUT={input.bam} "
       "OUTPUT={output.metrics} "
       "HISTOGRAM_FILE={output.histogram} "


















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
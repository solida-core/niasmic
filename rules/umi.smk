

rule umi_deduplication:
    input:
        bam="reads/merged/{sample}.bam",
        bai="reads/merged/{sample}.bam.bai"
    output:
        bam="reads/dedup/{sample}.dedup.bam",
        metrics="reads/dedup/{sample}.dedup_stats.txt"
    log:
        "logs/umi_tools/deduplication/{sample}.log"
    conda:
        "../envs/umi_tools.yaml"
    params:
        tmp_dir=java_params(tmp_dir=config.get("tmp_dir"), multiply_by=5)
    shell:
        "umi_tools dedup "
        "-I {input.bam} "
        "--paired "
        "-S {output} "
        "--method=unique " # default: directional
        "-L {log} "
        "--output-stats={output.metrics} "
        "--umi-separator=':' "
#        "--temp-dir={params.tmp_dir} "
        "--unmapped-reads=discard " # si potrebbe usare use (usa solo la r1 se la r2 non mappa


rule fgbio_annotate:
    input:
        bam="reads/merged/{sample}.bam",
        fastq="reads/merged/{sample}.bam.bai"
    output:
        bam="reads/fgbio/{sample}.umi.bam"
    conda:
        "../envs/fgbio.yaml"
    params:
        tmp_dir=java_params(tmp_dir=config.get("tmp_dir"), multiply_by=5)
    shell:
        "fgbio AnnotateBamWithUmis "
        "-i {input.bam} "
        "-f {input.fastq} "
        "-o {output} "

rule extract_umis:
    input:
        bam="reads/merged/{sample}.bam"
    output:
        bam="reads/fgbio/{sample}.umi.bam"
    conda:
        "../envs/fgbio.yaml"
    params:
        tmp_dir=java_params(tmp_dir=config.get("tmp_dir"), multiply_by=5)
    shell:
        "fgbio ExtractUmisFromBam "
        "-i {input.bam} "
        "--read-structure=XXXXXXXXX "
        "-o {output} "

rule fgbio_group:
    input:
        bam="reads/fgbio/{sample}.umi.bam"
    output:
        bam="reads/fgbio/{sample}.umi.grouped.bam"
    conda:
        "../envs/fgbio.yaml"
    params:
        tmp_dir=java_params(tmp_dir=config.get("tmp_dir"), multiply_by=5)
    shell:
        "fgbio AnnotateBamWithUmis "
        "-i {input.bam} "
        "-o {output} "
        "--raw-tag='RX' "
        "--min-map-q=30 "
        "--edits=1 "

rule fgbio_consensus:
    input:
        bam="reads/fgbio/{sample}.umi.grouped.bam"
    output:
        bam="reads/fgbio/{sample}.consensus.bam"
    conda:
        "../envs/fgbio.yaml"
    params:
        tmp_dir=java_params(tmp_dir=config.get("tmp_dir"), multiply_by=5)
    shell:
        "fgbio AnnotateBamWithUmis "
        "-i {input.bam} "
        "-o {output} "
        "--min-reads=???????????"##TODO
        "--tag='MI' "
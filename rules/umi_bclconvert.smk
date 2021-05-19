def get_fastq(wildcards,units):
    # print(wildcards.unit)
    if units.loc[wildcards.unit,["fq2"]].isna().all():
        print("SE")
        # print(units.loc[wildcards.unit,["fq1"]].dropna()[0])
        return expand_filepath(units.loc[wildcards.unit,["fq1"]].dropna()[0])
    else:
        print("PE")
        # print(units.loc[wildcards.unit,["fq1"]].dropna()[0],units.loc[wildcards.unit,["fq2"]].dropna()[0])
        return expand_filepath(units.loc[wildcards.unit,["fq1"]].dropna()[0]),expand_filepath(units.loc[wildcards.unit,["fq2"]].dropna()[0])



rule pre_rename_fastq_pe:
    input:
        lambda wildcards: get_fastq(wildcards,units)
    output:
        r1="reads/untrimmed/{unit}-R1.fastq.gz",
        r2="reads/untrimmed/{unit}-R2.fastq.gz"
    shell:
        "ln -s {input[0]} {output.r1} &&"
        "ln -s {input[1]} {output.r2} "






rule umi_group:
    input:
        "reads/merged/{sample}.bam",
        "reads/merged/{sample}.bam.bai"
    output:
        "reads/group/{sample}_grouped.bam"
    log:
        "logs/umi_tools/group/{sample}.log"
    conda:
        "../envs/umi_tools.yaml"
    params:
        tmp_dir=config.get("tmp_dir")
    shell:
        "umi_tools group "
        "-I {input[0]} "
        "--output-bam "
        "--paired "
        "-S {output} "
        "--chimeric-pairs=discard "
        "-L {log} "
        "--temp-dir={params.tmp_dir} "
        "-umi-separator=':' "


rule umi_dedup:
    input:
        "reads/group/{sample}_grouped.bam",
        "reads/group/{sample}_grouped.bam.bai"
    output:
        "reads/dedup/{sample}.dedup.bam"
    log:
        "logs/umi_tools/dedup/{sample}.log"
    conda:
        "../envs/umi_tools.yaml"
    params:
        tmp_dir=config.get("tmp_dir"),
        stats="qc/umitools/{sample}_dedup_stats"
    shell:
        # "mkdir -p qc/umitools ;"
        "umi_tools dedup "
        "-I {input[0]} "
        "-S {output} "
        "--output-stats {params.stats} "
        "--paired "
        "-L {log} "
        "--temp-dir={params.tmp_dir} "
        "--umi-tag='BX' "
        "--extract-umi-method tag "


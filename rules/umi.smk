def get_fastq(wildcards,units):
    # print(wildcards.unit)
    if units.loc[wildcards.unit,["fq2"]].isna().all():
        print("SE")
        # print(units.loc[wildcards.unit,["fq1"]].dropna()[0])
        return expand_filepath(units.loc[wildcards.unit,["fq1"]].dropna()[0])
    else:
        print("PE")
        # print(units.loc[wildcards.unit,["fq1"]].dropna()[0],units.loc[wildcards.unit,["fq2"]].dropna()[0])
        return expand_filepath(units.loc[wildcards.unit,["fq1"]].dropna()[0]),expand_filepath(units.loc[wildcards.unit,["fq2"]].dropna()[0]),expand_filepath(units.loc[wildcards.unit,["fq3"]].dropna()[0])

rule umi_annotation_R1:
    input:
        lambda wildcards: get_fastq(wildcards,units)
    output:
        r1="reads/untrimmed/{unit}-R1.fastq.gz"
    log:
        "logs/umi_tools/extraction/{unit}_R1.log"
    conda:
        "../envs/umi_tools.yaml"
    params:
        tmp_dir=java_params(tmp_dir=config.get("tmp_dir"), multiply_by=5)
    shell:
        "umi_tools extract "
        "--bc-pattern=NNNNNNNNN "
        "--stdin={input[1]} " # R2 with UMIs
        "--read2-in={input[0]} " # R1
        "--stdout={output.r1} "
        "--read2-stdout "
        "-L {log} "
        "--temp-dir={params.tmp_dir} "


rule umi_annotation_R2:
    input:
        lambda wildcards: get_fastq(wildcards,units)
    output:
        r2="reads/untrimmed/{unit}-R2.fastq.gz"
    log:
        "logs/umi_tools/extraction/{unit}_R2.log"
    conda:
        "../envs/umi_tools.yaml"
    params:
        tmp_dir=java_params(tmp_dir=config.get("tmp_dir"), multiply_by=5)
    shell:
        "umi_tools extract "
        "--bc-pattern=NNNNNNNNN "
        "--stdin={input[1]} " # R2 with UMIs
        "--read2-in={input[2]} " # R3
        "--stdout={output.r2} "
        "--read2-stdout "
        "-L {log} "
        "--temp-dir={params.tmp_dir} "


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
#        tmp_dir=java_params(tmp_dir=config.get("tmp_dir"), multiply_by=5)
        tmp_dir=config.get("tmp_dir")
    shell:
#        "mkdir -p reads/group ;"
        "umi_tools group "
        "-I {input[0]} "
        "--output-bam "
        "--paired "
        "-S {output} "
        "--chimeric-pairs=discard "
        "-L {log} "
        "--temp-dir={params.tmp_dir} "

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
        "mkdir -p qc/umitools ;"
        "umi_tools dedup "
        "-I {input[0]} "
        "-S {output} "
        "--output-stats {params.stats} "
        "-L {log} "
        "--temp-dir={params.tmp_dir} "


rule fgbio_group_reads:
    input:
        bam="reads/merged/{sample}.bam",
        bai="reads/merged/{sample}.bam.bai"
    output:
        bam="reads/fgbio/{sample}.umi.bam"
    conda:
        "../envs/fgbio.yaml"
    params:
        tmp_dir=java_params(tmp_dir=config.get("tmp_dir"), multiply_by=5),
        strategy="Identity", # Identity, Edit, Adjacency
        min_map_qual="30",
    shell:
        "fgbio GroupReadsByUmi "
        "--tmp-dir={params.tmp_dir} "
        "-s {params.strategy} "
        "-i {input.bam} "
        "-o {output} "
        "-m {params.min_map_qual} "


rule fgbio_call_consensus:
    input:
        bam="reads/fgbio/{sample}.umi.bam"
    output:
        bam="reads/fgbio/consensus/{sample}.umi.consensus.bam"
    conda:
        "../envs/fgbio.yaml"
    params:
        tmp_dir=java_params(tmp_dir=config.get("tmp_dir"), multiply_by=5),
        min_reads="3",#The minimum number of reads to produce a consensus base
    shell:
        "fgbio CallMolecularConsensusReads "
        "--tmp-dir={params.tmp_dir} "
        "-i {input.bam} "
        "-o {output} "
        "-M {params.min_reads} "


rule fgbio_filter_consensus:
    input:
        bam="reads/fgbio/consensus/{sample}.umi.consensus.bam"
    output:
        bam="reads/fgbio/dedup/{sample}.umi.dedup.bam"
    conda:
        "../envs/fgbio.yaml"
    params:
        tmp_dir=java_params(tmp_dir=config.get("tmp_dir"), multiply_by=5),
        genome=resolve_single_filepath( * references_abs_path(),config.get("genome_fasta")),
        min_reads="3",#The minimum number of reads to produce a consensus base
    shell:
        "fgbio FilterConsensusReads "
        "--tmp-dir={params.tmp_dir} "
        "-i {input.bam} "
        "-o {output} "
        "-M {params.min_reads} "
        "-r {params.genome} "
        "--max-read-error-rate=0.025 " # default
        "--max-base-error-rate=0.1 " # default


rule samtools_index_dedup_bam:
    input:
        "reads/fgbio/dedup/{sample}.umi.dedup.bam"
    output:
         "reads/fgbio/dedup/{sample}.umi.dedup.bam.bai"
    conda:
        "../envs/samtools.yaml"
    benchmark:
        "benchmarks/samtools/index_deduplicated/{sample}.txt"
    shell:
        "samtools index "
        "{input} "
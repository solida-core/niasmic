rule samtools_sort:
   input:
       "reads/aligned/{unit}_fixmate.cram"
   output:
       temp("reads/sorted/{unit}_sorted.cram")
   conda:
       "../envs/samtools.yaml"
   params:
       tmp_dir=tmp_path(path=config.get("tmp_dir")),
       genome=resolve_single_filepath(*references_abs_path(), config.get("genome_fasta")),
       output_fmt="CRAM"
   benchmark:
       "benchmarks/samtools/sort/{unit}.txt"
   threads: conservative_cpu_count(reserve_cores=2, max_cores=99)
   shell:
       "samtools sort "
       "--threads {threads} "
       "-T {params.tmp_dir} "
       "-O {params.output_fmt} "
       "--reference {params.genome} "
       "-o {output} "
       "{input} "


rule samtools_merge:
    input:
        lambda wildcards: get_units_by_sample(wildcards, samples,
                                              prefix='reads/sorted/',
                                              suffix='_sorted.cram')
    output:
        temp("reads/merged/{sample}.cram")
    conda:
        "../envs/samtools.yaml"
    benchmark:
        "benchmarks/samtools/merge/{sample}.txt"
    params:
        cmd='samtools',
        genome=resolve_single_filepath(*references_abs_path(), config.get("genome_fasta")),
        output_fmt="CRAM"
    threads: conservative_cpu_count(reserve_cores=2, max_cores=99)
    script:
        "scripts/samtools_merge.py"


rule samtools_cram_to_bam:
    input:
        "reads/merged/{sample}.cram"
    output:
        "reads/merged/{sample}.bam"
    conda:
        "../envs/samtools.yaml"
    benchmark:
        "benchmarks/samtools/cram_to_bam/{sample}.txt"
    params:
        genome=resolve_single_filepath(*references_abs_path(),
                                       config.get("genome_fasta")),
        output_fmt="BAM"
    threads: conservative_cpu_count(reserve_cores=2, max_cores=99)
    shell:
        "samtools view -b "
        "--threads {threads} "
        "-T {params.genome} "
        "-o {output} "
        "-O {params.output_fmt} "
        "{input} "

rule samtools_index_bam:
    input:
        "reads/merged/{sample}.bam"
    output:
         "reads/merged/{sample}.bam.bai"
    conda:
        "../envs/samtools.yaml"
    benchmark:
        "benchmarks/samtools/index_1/{sample}.txt"
    shell:
        "samtools index "
        "{input} "


rule samtools_index_cram:
    input:
        "reads/recalibrated/{sample}.dedup.recal.cram"
    output:
         "reads/recalibrated/{sample}.dedup.recal.cram.crai"
    conda:
        "../envs/samtools.yaml"
    benchmark:
        "benchmarks/samtools/index_2/{sample}.txt"
    shell:
        "samtools index "
        "{input} "


rule samtools_index_umigroup:
    input:
        "reads/group/{sample}_grouped.bam"
    output:
         "reads/group/{sample}_grouped.bam.bai"
    conda:
        "../envs/samtools.yaml"
    benchmark:
        "benchmarks/samtools/index_umigrouped/{sample}.txt"
    shell:
        "samtools index "
        "{input} "
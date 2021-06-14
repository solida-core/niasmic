
rule trim_galore_se:
    input:
        "reads/untrimmed/{unit}-R1.fastq.gz"
    output:
        "reads/trimmed/{unit}-R1_trimmed.fq.gz",
        "reads/trimmed/{unit}-R1.fastq.gz_trimming_report.txt"
    params:
        extra=config.get("rules").get("trim_galore_se").get("arguments"),
        outdir="reads/trimmed/"
    log:
        "logs/trim_galore/{unit}.log"
    benchmark:
        "benchmarks/trim_galore/{unit}.txt"
    conda:
        "../envs/trim_galore.yaml"
    threads: (conservative_cpu_count(reserve_cores=2, max_cores=99))/2 if (conservative_cpu_count(reserve_cores=2, max_cores=99)) >2 else 1
    shell:
        "mkdir -p qc/fastqc; "
        "trim_galore "
        "{params.extra} "
        "--cores {threads} "
        "-o {params.outdir} "
        "{input} "
        ">& {log}"


rule post_rename_fastq_se:
    input:
        rules.trim_galore_se.output
    output:
        r1="reads/trimmed/{unit}-R1-trimmed.fastq.gz"
    shell:
        "mv {input[0]} {output.r1} "


def get_trimmed_reads(wildcards,units):
    print(wildcards.unit)
    if units.loc[wildcards.unit,["fq2"]].isna().all():
        # SE
        return None
    # PE
    else:
        return rules.post_rename_fastq_pe.output.r1,rules.post_rename_fastq_pe.output.r2

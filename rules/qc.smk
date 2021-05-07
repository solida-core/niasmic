

rule multiqc:
    input:
#        expand("qc/fastqc/{unit.unit}-R1_fastqc.zip", unit=units.reset_index().itertuples()),
#        expand("qc/fastqc/trimmed_{unit.unit}_fastqc.zip", unit=units.reset_index().itertuples()),
#         expand("reads/trimmed/{unit.unit}-R1.100bp_5prime.fq.gz_trimming_report.txt", unit=units.reset_index().itertuples()),
        expand("qc/picard/UMIstats/{sample.sample}.duplicate_metrics.txt",sample=samples.reset_index().itertuples()),
        expand("qc/picard/hs/{sample.sample}.dedup.recal.hs.txt",sample=samples.reset_index().itertuples()),
        expand("qc/picard/hs/{sample.sample}.dedup.recal.is.txt",sample=samples.reset_index().itertuples()),
        expand("qc/picard/{sample.sample}_gc_bias_metrics.txt",sample=samples.reset_index().itertuples()),
        expand("qc/picard/{sample.sample}_summary_metrics.txt",sample=samples.reset_index().itertuples())
    output:
        "qc/multiqc.html"
    params:
        params=config.get("rules").get("multiqc").get("arguments"),
        outdir="qc",
        outname="multiqc.html",
        fastqc="qc/fastqc/",
        trimming="reads/trimmed/",
        reheader=config.get("reheader")
    conda:
        "../envs/multiqc.yaml"
    log:
        "logs/multiqc/multiqc.log"
    shell:
        "multiqc "
        "{input} "
        "{params.fastqc} "
        "{params.trimming} "
        "{params.params} "
        "-o {params.outdir} "
        "-n {params.outname} "
        "--sample-names {params.reheader} "
        ">& {log}"


# def fastqc_detect(wildcards,units):
#     if units.loc[wildcards.unit,["fq2"]].isna().all():
#         print("SE")
#         #print("reads/untrimmed/se/"+wildcards.unit+"-R1.fq.gz")
#         return "reads/untrimmed/se/"+wildcards.unit+"-R1.fq.gz"
#     else:
#         print("PE")
#         return "reads/untrimmed/"+wildcards.unit+"-R1.fq.gz"



# rule fastqc:
#     input:
#         lambda wildcards: fastqc_detect(wildcards, units)
#     output:
#         html="qc/fastqc/{unit}-R1.html",
#         zip="qc/fastqc/{unit}-R1_fastqc.zip"
#     log:
#         "logs/fastqc/untrimmed/{unit}.log"
#     params:
#         outdir="qc/fastqc"
#     conda:
#         "../envs/fastqc.yaml"
#     shell:
#         "fastqc "
#         "{input} "
#         "--outdir {params.outdir} "
#         "--quiet "
#         ">& {log}"

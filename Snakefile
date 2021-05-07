import pandas as pd
from snakemake.utils import validate, min_version
##### set minimum snakemake version #####
min_version("5.10.0")
##### load config and sample sheets #####


include:
    "rules/functions.py"

samples = pd.read_csv(config["samples"], index_col=["sample"],sep='\t')
units = pd.read_csv(config["units"], index_col=["unit"], dtype=str,sep='\t')
patients = pd.read_csv(config["patients"], index_col=["patient"],sep='\t')


###
localrules: all, pre_rename_fastq_pe, post_rename_fastq_pe

rule all:
    input:
        expand("reads/recalibrated/{sample.sample}.dedup.recal.bam", sample=samples.reset_index().itertuples()),
        # expand("reads/fgbio/{sample.sample}.consensus.bam", sample=samples.reset_index().itertuples()),
#        expand("data/results/{sample.patient}_somatic.vcf.gz",sample=samples.reset_index().itertuples()),
#        expand("data/results/{patient.patient}_somatic_twicefiltered_selected.vcf.gz", patient=patients.reset_index().itertuples()),
#         "qc/bedtools/heatmap_enriched_regions.png",
        # expand("annotation/{patient.patient}/bcftools/selected.annot.lightened.xlsx", patient=patients.reset_index().itertuples()),
        "qc/multiqc.html",
        expand("data/results/{patient.patient}_funcotated.maf", patient=patients.reset_index().itertuples()),
    #    expand("qc/bedtools/{sample.sample}.coverage.tsv", sample=samples.reset_index().itertuples()),
    #     expand("qc/{patient.patient}.multicov.tsv", patient=patients.reset_index().itertuples())



##### load rules #####
include_prefix="rules"
dima_path="dima/"
include:
    include_prefix + "/trimming.smk"
include:
    include_prefix + "/alignment.smk"
include:
    include_prefix + "/samtools.smk"
include:
    include_prefix + "/umi.smk"
include:
    include_prefix + "/bsqr.smk"
include:
    include_prefix + "/mutect_call.smk"
include:
    include_prefix + "/collect_filters.smk"
include:
    include_prefix + "/filter_mutect_call.smk"
include:
    include_prefix + "/annotation.smk"
include:
    include_prefix + "/format_output.smk"
include:
    include_prefix + "/qc.smk"
include:
    include_prefix + "/picard_stats.smk"
include:
    include_prefix + "/coverage.smk"

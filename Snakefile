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
sets = pd.read_csv(config["sets"], index_col=["set"], dtype=str, sep="\t")
plasma_samples = samples[samples["sample_type"]=="Pl"]



###
localrules: all, pre_rename_fastq_pe, post_rename_fastq_pe

rule all:
    input:
        expand("reads/recalibrated/{sample.sample}.dedup.recal.bam",sample=samples.reset_index().itertuples()),
        # expand("reads/fgbio/{sample.sample}.consensus.bam", sample=samples.reset_index().itertuples()),
        #        expand("data/results/{sample.patient}_somatic.vcf.gz",sample=samples.reset_index().itertuples()),
        #        expand("data/results/{patient.patient}_somatic_twicefiltered_selected.vcf.gz", patient=patients.reset_index().itertuples()),
        #         "qc/bedtools/heatmap_enriched_regions.png",
        "variant_calling/all.vcf.gz",
        "variant_calling/all.snp_recalibrated.indel_recalibrated.vcf.gz",
        # expand("annotation/{patient.patient}/bcftools/selected.annot.lightened.xlsx", patient=patients.reset_index().itertuples()),
        "qc/multiqc.html",
        expand("data/results/{patient.patient}_funcotated.maf",patient=patients.reset_index().itertuples()),
        #    expand("qc/bedtools/{sample.sample}.coverage.tsv", sample=samples.reset_index().itertuples()),
        #     expand("qc/{patient.patient}.multicov.tsv", patient=patients.reset_index().itertuples())
        #        "qc/picard/vcf_metrics.txt",
        expand("variant_calling/SelectVariants/{set.set}.selected.vcf",set=sets.reset_index().itertuples()),
#        expand("variant_calling/{sample.sample}.pileup.txt",sample=samples[samples["sample_type"]=="Pl"].reset_index().itertuples()),
#        expand("variant_calling/strelka/{patient.patient}.strelka.somatic.snvs.vcf.gz",patient=patients.reset_index().itertuples()),
#        expand("variant_calling/varscan2/{patient.patient}.varscan.snp.vcf",patient=patients.reset_index().itertuples()),
#        expand("variant_calling/vardict/{patient.patient}.vardict.vcf",patient=patients.reset_index().itertuples()),
        expand("somatic_combiner/{patient.patient}.somatic_combiner.vcf",patient=patients.reset_index().itertuples()),
        expand("somaticseq/{patient.patient}/Consensus.sSNV.vcf",patient=patients.reset_index().itertuples()),
        expand("somatic_combiner/{patient.patient}.somatic_combiner_annotated.flt.txt",patient=patients.reset_index().itertuples()),
        # expand("somaticseq/{patient.patient}.somaticseq_annotated.flt.txt",patient=patients.reset_index().itertuples()),
        expand("somaticseq/{patient.patient}_somaticseq.vcf",patient=patients.reset_index().itertuples()),
##### load rules #####
include_prefix="rules"
dima_path="dima/"
if config.get("pe")=="no":
    if config.get("fastq_numb")==2:
        include:
            include_prefix + "/umi_bclconvert_se.smk"
    else:
        include:
            include_prefix + "/umi_se.smk"
    include:
        include_prefix + "/trimming_se.smk"
    include:
        include_prefix + "/alignment_se.smk"



else:
    # include:
        # include_prefix + "/umi.smk"
    include:
        include_prefix + "/trimming.smk"
    include:
        include_prefix + "/alignment.smk"
    if config.get("fastq_numb")==2:
        include:
            include_prefix + "/umi_bclconvert.smk"
    else:
        include:
            include_prefix + "/umi.smk"
include:
    include_prefix + "/samtools.smk"
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
include:
    include_prefix + "/call_variants.smk"
include:
    include_prefix + "/joint_call.smk"
include:
    include_prefix + "/vsqr.smk"
include:
    include_prefix + "/vcf.smk"
include:
    include_prefix + "/plasma_pileup.smk"
include:
    include_prefix + "/other_callers.smk"
include:
    include_prefix + "/call_combiners.smk"
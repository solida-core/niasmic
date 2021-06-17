rule somatic_combiner:
    input:
        muse=rules.MuSE_sump.output,
        vardict=rules.vardict.output,
        strelka_snvs=rules.strelka.output.snvs,
        strelka_indels=rules.strelka.output.indels,
        mutect=rules.filter_mutect.output,
        varscan_snp=rules.varscan2.output.snp,
        varscan_indel=rules.varscan2.output.indel,
        lofreq_snvs=rules.Lofreq.output.snvs,
        lofreq_indels=rules.Lofreq.output.indels
    output:
        "somatic_combiner/{patient}.somatic_combiner.vcf"
    params:
        somatic_combiner=config.get("somatic_combiner")
    benchmark:
        "benchmarks/muse/{patient}.muse_sump.txt"
    threads: 2
    shell:
        "java -jar {params.somatic_combiner} "
        "--vardict {input.vardict} "
        "--mutect2 {input.mutect} "
        "--strelka-snv {input.strelka_snvs} "
        "--strelka-indel {input.strelka_indels} "
        "--varscan-snv {input.varscan_snp} "
        "--varscan-indel {input.varscan_indel} "
        "--muse {input.muse} "
        "--lofreq-snv {input.lofreq_snvs} "
        "--lofreq-indel {input.lofreq_indels} "
        "--output {output} "
        # "--threshold 0.01"


rule somaticseq:
    input:
        lambda wildcards: get_sample_by_famid(wildcards,patients),
        muse=rules.MuSE_sump.output,
        vardict=rules.vardict.output,
        strelka_snvs=rules.strelka.output.snvs,
        strelka_indels=rules.strelka.output.indels,
        mutect=rules.filter_mutect.output,
        varscan_snp=rules.varscan2.output.snp,
        varscan_indel=rules.varscan2.output.indel,
        lofreq_snvs=rules.Lofreq.output.snvs,
        lofreq_indels=rules.Lofreq.output.indels
    output:
        "somaticseq/{patient}/Consensus.sSNV.vcf",
        "somaticseq/{patient}/Consensus.sINDEL.vcf"
    params:
        somatic_combiner=config.get("somatic_combiner"),
        outdir="somaticseq/{patient}",
        target=config.get("interval_list"),
        genome=resolve_single_filepath(*references_abs_path(), config.get("genome_fasta")),
        dbsnp=config.get("dbsnp_gz")
    conda:
        "../envs/somaticseq.yaml"
    benchmark:
        "benchmarks/somaticseq/{patient}.txt"
    threads: 2
    shell:
        "somaticseq_parallel.py "
        "-outdir {params.outdir} "
        "-ref {params.genome} "
        "-dbsnp {params.dbsnp} "
        # "-cosmic "
        "--inclusion-region {params.target} "
        "paired "
        "--tumor-bam-file {input[0]} "
        "--normal-bam-file {input[1]} "
        "--mutect2-vcf {input.mutect} "
        "--varscan-snv {input.varscan_snp} "
        "--varscan-indel {input.varscan_indel} "
        "--vardict-vcf {input.vardict} "
        "--muse-vcf {input.muse} "
        "--lofreq-snv {input.lofreq_snvs} "
        "--lofreq-indel {input.lofreq_indels} "
        "--strelka-snv {input.strelka_snvs} "
        "--strelka-indel {input.strelka_indels} "


rule somaticseq_fix_vcf:
    input:
        snv="somaticseq/{patient}/Consensus.sSNV.vcf",
        indel="somaticseq/{patient}/Consensus.sINDEL.vcf"
    output:
        snv="somaticseq/{patient}/Consensus.sSNV_header.vcf",
        indel="somaticseq/{patient}/Consensus.sINDEL_header.vcf"
    params:
        dict="/ELS/els9/users/biosciences/references/ucsc/hg19/ucsc.hg19.dict"
    conda:
        "../envs/picard.yaml"
    threads: 2
    shell:
        "picard UpdateVcfSequenceDictionary "
        "I={input.snv} "
        "O={output.snv} "
        "SD={params.dict} ; "
        "picard UpdateVcfSequenceDictionary "
        "I={input.indel} "
        "O={output.indel} "
        "SD={params.dict} "


rule somaticseq_tabix_vcf:
    input:
        snv="somaticseq/{patient}/Consensus.sSNV_header.vcf",
        indel="somaticseq/{patient}/Consensus.sINDEL_header.vcf"
    output:
        snv="somaticseq/{patient}/Consensus.sSNV_header.vcf.gz",
        indel="somaticseq/{patient}/Consensus.sINDEL_header.vcf.gz",
        tbi_snv="somaticseq/{patient}/Consensus.sSNV_header.vcf.gz.tbi",
        tbi_indel="somaticseq/{patient}/Consensus.sINDEL_header.vcf.gz.tbi"
    params:
        dict="/ELS/els9/users/biosciences/references/ucsc/hg19/ucsc.hg19.dict"
    conda:
        "../envs/samtools.yaml"
    threads: 2
    shell:
        "bgzip {input.snv} && "
        "bgzip {input.indel} ; "
        "tabix -p vcf "
        "{output.snv} && "
        "tabix -p vcf "
        "{output.indel} "


rule somaticseq_concat:
    input:
        snv="somaticseq/{patient}/Consensus.sSNV_header.vcf.gz",
        indel="somaticseq/{patient}/Consensus.sINDEL_header.vcf.gz",
        tbi_snv="somaticseq/{patient}/Consensus.sSNV_header.vcf.gz.tbi",
        tbi_indel="somaticseq/{patient}/Consensus.sINDEL_header.vcf.gz.tbi"
    output:
        "somaticseq/{patient}_somaticseq.vcf"
    conda:
        "../envs/bcftools.yaml"
    threads: 2
    shell:
        "bcftools concat "
        "{input.snv} "
        "{input.indel} "
        "-o {output} "
        "-O v -a"
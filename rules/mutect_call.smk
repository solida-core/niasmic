rule somatic_discovery:
    input:
        lambda wildcards: get_sample_by_famid(wildcards, patients)
    output:
        vcf="data/results/{patient}_somatic.vcf.gz",
        bam="data/results/{patient}_tumor_normal.bam",
        fir="data/results/{patient}_tumor_normal_f1r2.tar.gz"
    params:
        custom=java_params(tmp_dir=config.get("tmp_dir"), multiply_by=5),
        genome=resolve_single_filepath(*references_abs_path(), config.get("genome_fasta")),
        intervals=config.get("interval_list"),
        param=config.get("rules").get("mutect").get("params"),
        # pon=config.get("pon"),
        germline_resource=config.get("germline"),
        tbam=lambda wildcards, input: input[0].split("/")[-1].split(".")[0],
        cbam = lambda wildcards, input: input[1].split("/")[-1].split(".")[0]
    log:
        "logs/gatk/Mutect2/{patient}.somatic_call.log"
    conda:
        "../envs/gatk.yaml"
    threads: conservative_cpu_count(reserve_cores=2, max_cores=99)

    shell:
        "gatk "
        "--java-options {params.custom} "
        "Mutect2 "
        "-R {params.genome} "
        "-I {input[0]} "#tumoral bam
        "-I {input[1]} "#control bam
        "-normal {params.cbam} "#control name
#        "-pon {params.pon} "
        "--germline-resource {params.germline_resource} "
        "--af-of-alleles-not-in-resource 0.0000025 "
        "{params.param} "
        "-L {params.intervals} "
        "-O {output.vcf} "
#        "-bamout {output.bam} "
        "--native-pair-hmm-threads {threads} "
        "--max-reads-per-alignment-start 0 "
        "--f1r2-tar-gz {output.fir} "
        ">& {log} "
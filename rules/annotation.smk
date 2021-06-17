
rule kggseq_somatic_combiner:
    input:
        vcf="somatic_combiner/{patient}.somatic_combiner.vcf"
    output:
        "somatic_combiner/{patient}.somatic_combiner_annotated.flt.txt"
    params:
        custom=java_params(tmp_dir=config.get("tmp_dir"), multiply_by=5),
        cmd=config.get("rules").get("kggseq").get("cmd"),
        out_basename='somatic_combiner/{patient}.somatic_combiner_annotated'
    benchmark:
        "benchmarks/kggseq_sc/{patient}.kggseq.txt"
    threads: conservative_cpu_count()
    shell:
        "{params.cmd} -nt {threads} {params.custom} "
        "--no-resource-check --no-lib-check --no-web --no-gz "
        "--no-qc --o-vcf "
        "--db-gene refgene --cosmic-annot --cancer-mut-predict "
        "--vcf-file {input.vcf} "
        "--out {params.out_basename}"


rule kggseq_somaticseq:
    input:
        vcf="somaticseq/{patient}_somaticseq.vcf"
    output:
        "somaticseq/{patient}.somaticseq_annotated.flt.txt"
    params:
        custom=java_params(tmp_dir=config.get("tmp_dir"), multiply_by=5),
        cmd=config.get("rules").get("kggseq").get("cmd"),
        out_basename='somaticseq/{patient}.somaticseq_annotated'
    benchmark:
        "benchmarks/kggseq_ssq/{patient}.kggseq.txt"
    threads: conservative_cpu_count()
    shell:
        "{params.cmd} -nt {threads} {params.custom} "
        "--no-resource-check --no-lib-check --no-web --no-gz "
        "--no-qc --o-vcf "
        "--db-gene refgene --cosmic-annot --cancer-mut-predict "
        "--vcf-file {input.vcf} "
        "--out {params.out_basename}"
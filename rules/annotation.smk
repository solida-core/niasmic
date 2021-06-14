
rule kggseq:
    input:
        vcf="data/results/{patient}_somatic_twicefiltered_selected.vcf.gz"
    output:
        vcf='annotation/{patient}/kggseq/selected.flt.vcf',
        log='annotation/{patient}/kggseq/selected.log',
        txt='annotation/{patient}/kggseq/selected.flt.txt',
        ped='annotation/{patient}/kggseq/selected.ped'
    params:
        custom=java_params(tmp_dir=config.get("tmp_dir"), multiply_by=5),
        cmd=config.get("rules").get("kggseq").get("cmd"),
        arguments=_multi_flag(config.get("rules").get("kggseq").get("arguments")),
        ped_file=config.get("rules").get("kggseq").get("ped_file"),
        out_basename='annotation/{patient}/kggseq/selected'
    benchmark:
        "benchmarks/kggseq/{patient}.kggseq.txt"
    threads: conservative_cpu_count()
    shell:
        "cp {params.ped_file} {output.ped} && "
        "{params.cmd} -nt {threads} {params.custom} "
        "--no-resource-check --no-lib-check --no-web --no-gz "
        "--no-qc --o-vcf "
        "{params.arguments} "
        "--vcf-file {input.vcf} "
        "--ped-file {output.ped} "
        "--out {params.out_basename}"



rule vep:
    input:
        vcf="data/results/{patient}_somatic_twicefiltered_selected.vcf.gz"
    output:
        vcf='annotation/{patient}/vep/{patient}.vep.vcf.gz',
        stats='annotation/{patient}/vep/{patient}.stats.html'
    params:
        arguments=_multi_flag(config.get("rules").get("vep").get("arguments")),
        species=config.get("rules").get("vep").get("species"),
        assembly=config.get("rules").get("vep").get("assembly"),
        cache_dir=config.get("rules").get("vep").get("cache_dir"),
        log='annotation/{patient}/vep/warnings_vcf.log'
    benchmark:
        "benchmarks/vep/{patient}.vep.txt"
    conda:
       "../envs/vep.yaml"
    threads: conservative_cpu_count(reserve_cores=2, max_cores=6)
    shell:
        "vep --cache --format vcf --vcf --offline --pick "
        "-i {input.vcf} "
        "-o {output.vcf} "
        "{params.arguments} "
        "--fork {threads} "
        "--warning_file {params.log} "
        "--stats_file {output.stats} "
        "--species {params.species} "
        "--assembly {params.assembly} "
        "--compress_output gzip "
        "--dir_cache {params.cache_dir}"

rule filter_vep:
    input:
        rules.vep.output
    output:
        temp('annotation/{patient}/vep/{patient}.vep.snpsift.filt.vcf')
    params:
        arguments=_multi_flag(config.get("rules").get("filter_vep").get("arguments"))
    benchmark:
        "benchmarks/vep/{patient}.filter_vep.txt"
    conda:
       "../envs/vep.yaml"
    threads: conservative_cpu_count(reserve_cores=2, max_cores=6)
    shell:
        "filter_vep --format vcf "
        "-i {input} "
        "-o {output} "
        "{params.arguments}"

# Remove fields from VCF
rule bcftools_remove:
    input:
       rules.filter_vep.output
    output:
       temp('annotation/{patient}/vep/{patient}.vep.gtOnly.vcf')
    conda:
        "../envs/bcftools.yaml"
    params:
        fields_to_remove=config.get("rules").get("bcftools_remove").get("fields_to_remove")
    shell:
        "bcftools annotate "
        "-x {params.fields_to_remove} "
        "-o {output} "
        "{input}"

rule vcf_to_tabular:
    input:
       rules.bcftools_remove.output
    output:
       temp('annotation/{patient}/vep/{patient}.vep.gtOnly.tsv')
    params:
       script='../rules/scripts/vcf_to_tabular_futurized.py',
       params=config.get("rules").get("vcf_to_tabular").get("params")
    conda:
        "../envs/future.yaml"
    shell:
        "python {params.script} {params.params} {input} {output}"

# Merge filtered VCF and VEP annotations
rule vep_to_tsv:
    input:
        vep=rules.filter_vep.output,
        tsv=rules.vcf_to_tabular.output
    output:
        'annotation/{patient}/vep/{patient}.annotated.tsv'
    params:
        fields=config.get("rules").get("vep_to_tsv").get("fields")
    benchmark:
        "benchmarks/vep/{patient}.vep_to_tsv.txt"
    conda:
       "../envs/vatools.yaml"
    threads: conservative_cpu_count(reserve_cores=2, max_cores=6)
    shell:
        "vep-annotation-reporter "
        "-t {input.tsv} "
        "-o {output} "
        "{input.vep} "
        "{params.fields}"


rule tsv_to_excel:
    input:
       rules.vep_to_tsv.output
    output:
       'annotation/{patient}/vep/{patient}.annotated.xlsx'
    params:
       script='../rules/scripts/tabular_to_excel.py'
    conda:
        "../envs/excel.yaml"
    shell:
       "python {params.script} "
       "-i {input} "


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
#        "java -Djava.awt.headless=true -jar /ELS/els9/users/biosciences/software/kggseq-1.0/kggseq.jar "
        "{params.cmd} -nt {threads} {params.custom} "
        "--no-resource-check --no-lib-check --no-web --no-gz "
	"--no-qc --o-vcf "
        "--db-gene refgene --cosmic-annot --cancer-mut-predict "
        "--vcf-file {input.vcf} "
        "--out {params.out_basename}"



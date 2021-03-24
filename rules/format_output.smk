rule format_annotation:
    input:
       'annotation/{patient}/kggseq/selected.flt.txt'
    output:
       cname='annotation/{patient}/kggseq/annot.cname',
       header='annotation/{patient}/kggseq/annot.header',
       tab='annotation/{patient}/kggseq/annot.tab'
    params:
        blocks=config.get("rules").get("format_annotation").get("blocks_file")
    script:
        "scripts/format_annotation.py"


rule tabix:
    "Bgzip-compressed and tabix-indexed file with annotations"
    input:
       'annotation/{patient}/kggseq/annot.tab'
    output:
       'annotation/{patient}/kggseq/annot.tab.gz'
    conda:
        "../envs/tabix.yaml"
    params:
       config.get("rules").get("tabix").get("params")
    shell:
        "bgzip {input}; tabix {params} {output}"



rule bcftools_annotate_add:
    input:
       cname='annotation/{patient}/kggseq/annot.cname',
       header='annotation/{patient}/kggseq/annot.header',
       gz='annotation/{patient}/kggseq/annot.tab.gz',
       vcf='annotation/{patient}/kggseq/selected.flt.vcf'
    output:
       'annotation/{patient}/bcftools/selected.annot.vcf'
    conda:
        "../envs/bcftools.yaml"
    params:
        cmd='add'
    script:
        "scripts/bcftools_annotate.py"


rule bcftools_annotate_remove:
    input:
       'annotation/{patient}/bcftools/selected.annot.vcf'
    output:
       'annotation/{patient}/bcftools/selected.annot.lightened.vcf'
    conda:
        "../envs/bcftools.yaml"
    params:
        cmd='remove',
        blocks=config.get("rules").get("format_annotation").get("blocks_file"),
        params=config.get("rules").get("bcftools_annotate_remove").get("params")
    script:
        "scripts/bcftools_annotate.py"


rule vcf_to_tabular_kggseq:
    input:
       'annotation/{patient}/bcftools/selected.annot.lightened.vcf'
    output:
       'annotation/{patient}/bcftools/selected.annot.lightened.tsv'
    params:
       script='../rules/scripts/vcf_to_tabular_futurized.py',
       params= "--do-not-split-sample --print-format "

    conda:
        "../envs/future.yaml"
    shell:
        "python {params.script} {params.params} {input} {output}"

rule tabular_to_excel_full_kggseq:
    input:
       'annotation/{patient}/bcftools/selected.annot.lightened.tsv'
    output:
       'annotation/{patient}/bcftools/selected.annot.lightened.xlsx'
    params:
       script='../rules/scripts/tabular_to_excel.py'
    conda:
        "../envs/excel.yaml"
    shell:
       "python {params.script} "
       "-i {input} "
       "-o {output}"


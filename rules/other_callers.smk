rule strelka:
    input:    
        lambda wildcards: get_sample_by_famid(wildcards, patients)
    output:
        snvs="variant_calling/strelka/{patient}.strelka.somatic.snvs.vcf.gz",
        indels="variant_calling/strelka/{patient}.strelka.somatic.indels.vcf.gz"
    conda:
       "../envs/strelka.yaml"
    params:
        custom=java_params(tmp_dir=config.get("tmp_dir"), multiply_by=5),
        intervals = config.get("interval_list_gz"),
        genome=resolve_single_filepath(*references_abs_path(), config.get("genome_fasta")),
        outdir="strelka/{patient}"
    benchmark:
        "benchmarks/strelka/{patient}.txt"
    threads: 2
    shell:
        "exist_dir({params.outdir}, delete=True) ;"
        "configureStrelkaSomaticWorkflow.py "
        "--normalBam {input[1]} " 
        "--tumorBam {input[0]} "
        "--referenceFasta {params.genome} "
        "--outputCallableRegions --targeted "
        "--callRegions {params.intervals} "
        "--runDir {params.outdir} ; "
        "python {params.outdir}/runWorkflow.py -m local -g 10 ; "
        "cp {params.outdir}/results/variants/somatic.snvs.vcf.gz {output.snvs} && "
        "cp {params.outdir}/results/variants/somatic.indels.vcf.gz {output.indels} "



rule vardict:
    input:
        lambda wildcards: get_sample_by_famid(wildcards, patients)
    output:
        "variant_calling/vardict/{patient}.vardict.vcf"
    conda:
       "../envs/vardict.yaml"
    params:
        custom=java_params(tmp_dir=config.get("tmp_dir"), multiply_by=5),
        intervals = config.get("interval_list"),
        genome=resolve_single_filepath(*references_abs_path(), config.get("genome_fasta")),
        outdir="vardict/{patient}",
        tbam= lambda wildcards,input: input[0].split("/")[-1].split(".")[0],
        cbam= lambda wildcards,input: input[1].split("/")[-1].split(".")[0]
    benchmark:
        "benchmarks/vardict/{patient}.txt"
    threads: 2
    shell:
        "vardict "
        "-G {params.genome} "
        "-f 0.01 -N {params.tbam} "
        "-b '{input[0]} | {input[1]}' "
        "-c 1 -S 2 -E 3 "
        "{params.intervals} "
        # "> DNA21-0414_vs_germinal "
        "| testsomatic.R | var2vcf_paired.pl "
        "-N '{params.tbam}|{params.cbam}' "
        "-f 0.01 > {output} "


rule pre_varscan2_tumoral:
    input:
        lambda wildcards: get_sample_by_famid(wildcards, patients)
    output:
        "variant_calling/samtools_pileup/{patient}.tumoral.pileup"
    conda:
       "../envs/samtools.yaml"
    params:
        custom=java_params(tmp_dir=config.get("tmp_dir"), multiply_by=5),
        intervals = config.get("interval_list"),
        genome=resolve_single_filepath(*references_abs_path(), config.get("genome_fasta")),
        tbam= lambda wildcards,input: input[0].split("/")[-1].split(".")[0],
        cbam= lambda wildcards,input: input[1].split("/")[-1].split(".")[0]
    threads: 2
    shell:
        "samtools "
        "mpileup "
        "-f {params.genome} "
        "-x -C 50 -q 40 -Q 30 "
        "-l {params.intervals} "
        "{input[0]} "
        "-o {output} "


rule pre_varscan2_germinal:
    input:
        lambda wildcards: get_sample_by_famid(wildcards, patients)
    output:
        "variant_calling/samtools_pileup/{patient}.germinal.pileup"
    conda:
       "../envs/samtools.yaml"
    params:
        custom=java_params(tmp_dir=config.get("tmp_dir"), multiply_by=5),
        intervals = config.get("interval_list"),
        genome=resolve_single_filepath(*references_abs_path(), config.get("genome_fasta")),
        tbam= lambda wildcards,input: input[0].split("/")[-1].split(".")[0],
        cbam= lambda wildcards,input: input[1].split("/")[-1].split(".")[0]
    threads: 2
    shell:
        "samtools "
        "mpileup "
        "-f {params.genome} "
        "-x -C 50 -q 40 -Q 30 "
        "-l {params.intervals} "
        "{input[1]} "
        "-o {output} "


rule varscan2:
    input:
        rules.pre_varscan2_tumoral.output,
        rules.pre_varscan2_germinal.output
    output:
        snp="variant_calling/varscan2/{patient}.varscan.snp.vcf",
        indel="variant_calling/varscan2/{patient}.varscan.indel.vcf",
    conda:
       "../envs/varscan2.yaml"
    params:
        custom=java_params(tmp_dir=config.get("tmp_dir"), multiply_by=5),
        intervals = config.get("interval_list"),
        genome=resolve_single_filepath(*references_abs_path(), config.get("genome_fasta"))
    benchmark:
        "benchmarks/varscan/{patient}.txt"
    threads: 2
    shell:
        "varscan somatic "
        "{input[1]} " # normal samtools pileup
        "{input[0]} " # tumoral samtools pileup
        "--output-vcf "
        "--tumor-purity 0.2 "
        "--output-snp {output.snp} "
        "--output-indel {output.indel}"


rule MuSE_call:
    input:
        lambda wildcards: get_sample_by_famid(wildcards, patients)
    output:
        "variant_calling/muse/{patient}.MuSE.txt"
    conda:
       "../envs/muse.yaml"
    params:
        custom=java_params(tmp_dir=config.get("tmp_dir"), multiply_by=5),
        intervals = config.get("interval_list"),
        genome=resolve_single_filepath(*references_abs_path(), config.get("genome_fasta")),
        out="variant_calling/muse/{patient}"
    benchmark:
        "benchmarks/muse/{patient}.muse_call.txt"
    threads: 2
    shell:
        "MuSE call "
        "-f {params.genome} "
        "-l {params.intervals} "
        "{input[0]} " ## tumoral bam (positional)
        "{input[1]} " ## normal bam
        "-O {params.out} "


rule MuSE_sump:
    input:
        rules.MuSE_call.output
    output:
        "variant_calling/muse/{patient}.muse.vcf"
    conda:
       "../envs/muse.yaml"
    params:
        custom=java_params(tmp_dir=config.get("tmp_dir"), multiply_by=5),
        intervals = config.get("interval_list"),
        genome=resolve_single_filepath(*references_abs_path(), config.get("genome_fasta")),
        dbsnp=config.get("dbsnp_gz")
    benchmark:
        "benchmarks/muse/{patient}.muse_sump.txt"
    threads: 2
    shell:
        "MuSE sump "
        "-I {input} " 
        "-D {params.dbsnp} "
        "-E "
        "-O {output} "


rule Lofreq:
    input:
        lambda wildcards: get_sample_by_famid(wildcards, patients)
    output:
        snvs="variant_calling/lofreq/{patient}.somatic_final_minus-dbsnp.snvs.vcf.gz",
        indels="variant_calling/lofreq/{patient}.somatic_final_minus-dbsnp.indels.vcf.gz"
    conda:
       "../envs/lofreq.yaml"
    params:
        custom=java_params(tmp_dir=config.get("tmp_dir"), multiply_by=5),
        intervals = config.get("interval_list"),
        genome=resolve_single_filepath(*references_abs_path(), config.get("genome_fasta")),
        out="variant_calling/lofreq/{patient}.",
        dbsnp=config.get("dbsnp_gz")
    benchmark:
        "benchmarks/lofreq/{patient}.lofreq_somatic.txt"
    threads: 2
    shell:
        "lofreq somatic "
        "-f {params.genome} "
        "-l {params.intervals} "
        "-t {input[0]} " ## tumoral bam 
        "-n {input[1]} " ## normal bam
        "-o {params.out} "
        "-d {params.dbsnp} "


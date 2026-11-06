"""
Module RSL_counter.

Per screen replicate (sample, e.g. RWPE_R1; two time points each):
CrisprCounter RSL counting -> filtered/collapsed log2FC table -> RSL binning ->
SSMD scores (RSLC IRA-SSMD). Requires the library from module crispr_library.
"""


rule create_sample_table:
    output:
        SAMPLE_TABLE,
    params:
        rows=lambda wc: crispr_rows(wc.sample),
    log:
        os.path.join(LOGS, "{cell_line}", "create_sample_table_{sample}.log"),
    conda:
        "../envs/crispr.yaml"
    script:
        "../scripts/make_crisprcounter_sample_table.py"


rule make_config_crisprcounter:
    input:
        samples=SAMPLE_TABLE,
        library=LIBRARY,
        template=config["crispr"]["crisprcounter"]["template"],
    output:
        COUNTER_PROPERTIES,
    params:
        fastqdir=lambda wc: CELL_CFG[wc.cell_line]["fastq_dir"],
        prefix=RSL_COUNTS,
        outdir=lambda wc, output: os.path.dirname(output[0]),
        script=script("makeCounterConfig.py"),
    log:
        os.path.join(LOGS, "{cell_line}", "make_config_crisprcounter_{sample}.log"),
    conda:
        "../envs/crispr.yaml"
    shell:
        """
        python {params.script} --template {input.template} \
            --fastqdir {params.fastqdir} \
            --samples {input.samples} \
            --library {input.library} \
            --prefix {params.prefix} \
            --outdir {params.outdir} > {log} 2>&1
        """


rule count_RSL:
    input:
        properties=COUNTER_PROPERTIES,
        jar=config["crispr"]["crisprcounter"]["jar"],
        fastq=counter_fastqs,
    output:
        RSL_COUNTS,
    params:
        mem=config["crispr"]["crisprcounter"]["java_mem_gb"],
    resources:
        mem_mb=config["crispr"]["crisprcounter"]["java_mem_gb"] * 1024 + 2048,
    log:
        os.path.join(LOGS, "{cell_line}", "crisprcounter_{sample}.log"),
    conda:
        "../envs/crispr.yaml"
    shell:
        "java -Xmx{params.mem}g -jar {input.jar} {input.properties} &> {log}"


rule process_RSL:
    input:
        RSL_COUNTS,
    output:
        PROCESSED_COUNTS,
    params:
        col1=lambda wc: f"{wc.sample}_{STAGES[0]}",
        col2=lambda wc: f"{wc.sample}_{STAGES[1]}",
        script=script("process_RSL_data.py"),
    log:
        os.path.join(LOGS, "{cell_line}", "process_RSL_{sample}.log"),
    conda:
        "../envs/crispr.yaml"
    shell:
        "python {params.script} {input} {output} {params.col1} {params.col2} > {log} 2>&1"


rule bin_RSL:
    input:
        RSL_COUNTS,
    output:
        RSL_COUNTS_BINNED,
    params:
        nr_nt=config["crispr"]["rsl_bin_nt"],
        min_count=config["crispr"]["rsl_bin_min_count"],
        script=script("bin_truncated_rsl_counts.pl"),
    log:
        os.path.join(LOGS, "{cell_line}", "bin_RSL_{sample}.log"),
    conda:
        "../envs/crispr.yaml"
    shell:
        "perl {params.script} {params.nr_nt} {params.min_count} {input} > {output} 2> {log}"


rule calculate_SSMD:
    input:
        RSL_COUNTS_BINNED,
    output:
        scores=SSMD_SCORES,
        normalized=SSMD_NORMALIZED,
    params:
        prefix=lambda wc, output: output.scores.removesuffix("_SSMD_Score.tsv"),
        min_control_count=config["crispr"]["ssmd"]["min_control_count"],
        filter_on_sum=config["crispr"]["ssmd"]["filter_on_sum"],
        negcon=lambda wc: CELL_CFG[wc.cell_line]["neg_con"],
        header=lambda wc: f"RSL.guide,guide.set,Control,{wc.sample}_{STAGES[1]}",
        script=script("IRA-SSMD_modded.R"),
    log:
        os.path.join(LOGS, "{cell_line}", "calculate_SSMD_{sample}.log"),
    conda:
        "../envs/rslc.yaml"
    shell:
        """
        sed '1s/.*/{params.header}/' {input} | \
            Rscript --vanilla {params.script} \
            /dev/stdin {params.prefix} \
            {params.min_control_count} \
            {params.filter_on_sum} \
            {params.negcon} > {log} 2>&1
        """

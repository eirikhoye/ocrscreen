"""
Module JACKS.

Per cell line: RSL bins of each screen replicate are treated as synthetic replicates
and analysed with JACKS (OCR = "gene"). The JACKS OCR scores/p-values are turned into
CRISPRi enhancer-activity BED files (FDR-filtered, unfiltered, narrowPeak) that can be
used as activity input for ABC/STARE. Requires RSL_counter and crispr_library outputs.
"""


rule aggregate_RSL:
    input:
        binned_counts_of_cell,
    output:
        AGGREGATE_RSL,
    log:
        os.path.join(LOGS, "{cell_line}", "aggregate_RSL.log"),
    conda:
        "../envs/crispr.yaml"
    script:
        "../scripts/aggregate_rsl.py"


rule convert_RSL_to_synthetic_replicates:
    input:
        binned=binned_counts_of_cell,
        ref=LIBRARY,
    output:
        rep=RSL_SYNTHETIC_REPLICATES,
        repmap=REPMAP,
    params:
        nr_nt=config["crispr"]["rsl_bin_nt"],
        rows=lambda wc: crispr_rows_of(wc.cell_line),
    log:
        os.path.join(LOGS, "{cell_line}", "convert_RSL_to_synthetic_replicates.log"),
    conda:
        "../envs/crispr.yaml"
    script:
        "../scripts/rsl_synthetic_replicates.py"


rule JACKS:
    input:
        counts=RSL_SYNTHETIC_REPLICATES,
        repmap=REPMAP,
        ref=LIBRARY,
        run_jacks=RUN_JACKS,
    output:
        JACKS_OUTPUTS,
    params:
        ctrl_sample=config["jacks"]["ctrl_sample"],
        sgrna_hdr=config["jacks"]["sgrna_hdr"],
        gene_hdr=config["jacks"]["gene_hdr"],
        outprefix=lambda wc, output: output[0].removesuffix(JACKS_SUFFIXES[0]),
        ctrl_genes=lambda wc: CELL_CFG[wc.cell_line]["neg_con"],
    log:
        os.path.join(LOGS, "{cell_line}", "JACKS.log"),
    conda:
        "../envs/jacks.yaml"
    shell:
        """
        python -s {input.run_jacks} {input.counts} {input.repmap} {input.ref} \
            --ctrl_sample {params.ctrl_sample} \
            --sgrna_hdr {params.sgrna_hdr} \
            --gene_hdr {params.gene_hdr} \
            --outprefix={params.outprefix} \
            --ctrl_genes={params.ctrl_genes} > {log} 2>&1
        """


rule make_enhancer_activity_bed_filtered:
    input:
        results=JACKS_GENE,
        pvals=JACKS_PVAL,
        coord_ref=LIBRARY_COORDS,
    output:
        bed=ENHANCER_BED_FILTERED,
    params:
        fdr=config["jacks"]["fdr_threshold"],
        script=script("make_enhancer_activity_bed.py"),
    log:
        os.path.join(LOGS, "{cell_line}", "make_enhancer_activity_bed_filtered.log"),
    conda:
        "../envs/crispr.yaml"
    shell:
        """
        python {params.script} -r {input.results} -p {input.pvals} -c {input.coord_ref} \
            -n {wildcards.cell_line} -o {output.bed} -t {params.fdr} > {log} 2>&1
        """


rule make_enhancer_activity_bed_unfiltered:
    input:
        results=JACKS_GENE,
        pvals=JACKS_PVAL,
        coord_ref=LIBRARY_COORDS,
    output:
        bed=ENHANCER_BED_UNFILTERED,
    params:
        fdr=config["jacks"]["fdr_threshold"],
        script=script("make_enhancer_activity_bed.py"),
    log:
        os.path.join(LOGS, "{cell_line}", "make_enhancer_activity_bed_unfiltered.log"),
    conda:
        "../envs/crispr.yaml"
    shell:
        """
        python {params.script} -r {input.results} -p {input.pvals} -c {input.coord_ref} \
            -n {wildcards.cell_line} -o {output.bed} -t {params.fdr} --no_fdr_filter > {log} 2>&1
        """


rule make_narrowpeaks_bed:
    input:
        results=JACKS_GENE,
        pvals=JACKS_PVAL,
        coord_ref=LIBRARY_COORDS,
    output:
        bed=ENHANCER_NARROWPEAK,
    params:
        fdr=config["jacks"]["fdr_threshold"],
        script=script("make_enhancer_narrowpeaks.py"),
    log:
        os.path.join(LOGS, "{cell_line}", "make_narrowpeaks_bed.log"),
    conda:
        "../envs/crispr.yaml"
    shell:
        """
        python {params.script} -r {input.results} -p {input.pvals} -c {input.coord_ref} \
            -n {wildcards.cell_line} -o {output.bed} -t {params.fdr} > {log} 2>&1
        """

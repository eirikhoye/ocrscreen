"""
Module GSEA (was crispr_pipeline RUN_GSEA.smk).

ABC interactions of the run config gsea.abc_run are restricted to OCRs with JACKS
FDR < jacks.fdr_threshold; genes are ranked by their summed ABC score and tested with
GSEApy prerank against each Enrichr library in gsea.gene_sets (needs internet access).
Requires JACKS and ABC_STARE outputs.
"""


rule filter_ABC_for_JACKS_fdr:
    input:
        abc=lambda wc: ABC_SCORES.format(cell_line=wc.cell_line, abc_run=config["gsea"]["abc_run"]),
        fdr=JACKS_PVAL,
        bed=ENHANCER_BED_UNFILTERED,
    output:
        ABC_FILTERED_FDR,
    params:
        fdr=config["jacks"]["fdr_threshold"],
    log:
        os.path.join(LOGS, "{cell_line}", "filter_ABC_for_JACKS_fdr.log"),
    conda:
        "../envs/crispr.yaml"
    script:
        "../scripts/filter_abc_by_jacks_fdr.py"


rule rank_genes_for_GSEA:
    input:
        ABC_FILTERED_FDR,
    output:
        ABC_RANK,
    log:
        os.path.join(LOGS, "{cell_line}", "rank_genes_for_GSEA.log"),
    conda:
        "../envs/crispr.yaml"
    script:
        "../scripts/rank_genes_for_gsea.py"


rule run_GSEA:
    input:
        rnk=ABC_RANK,
    output:
        GSEA_RES,
    params:
        outdir=lambda wc, output: os.path.dirname(output[0]),
        permutations=config["gsea"]["permutations"],
        script=script("run_gsea.py"),
    log:
        os.path.join(LOGS, "{cell_line}", "gsea_{gene_set}.log"),
    conda:
        "../envs/gseapy.yaml"
    shell:
        """
        python {params.script} \
            --input {input.rnk} \
            --outdir {params.outdir} \
            --geneset {wildcards.gene_set} \
            --permutations {params.permutations} > {log} 2>&1
        """

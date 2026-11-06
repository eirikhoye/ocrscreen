"""
Module pathway.

OCR-level enrichment: guides are grouped into one "pathway" per OCR (GMT) and tested
with MAGeCK pathway (GSEA mode) on the SSMD ranking of each screen replicate. Also
plots the SSMD/Score and enrichment-score replicate concordance per cell line.
Requires outputs of RSL_counter.
"""


rule make_gmt:
    input:
        counts=PROCESSED_COUNTS,
        coords=lambda wc: CELL_CFG[wc.cell_line]["coords"],
    output:
        GMT_FILE,
    params:
        neg_con=lambda wc: CELL_CFG[wc.cell_line]["neg_con"],
        script=script("convert_RSL_counts_to_gmt.py"),
    log:
        os.path.join(LOGS, "{cell_line}", "make_gmt_{sample}.log"),
    conda:
        "../envs/crispr.yaml"
    shell:
        """
        python {params.script} {input.counts} {input.coords} {wildcards.cell_line} \
            {params.neg_con} {output} > {log} 2>&1
        """


rule mageck_pathway:
    input:
        gene_ranking=SSMD_SCORES,
        gmt=GMT_FILE,
    output:
        PATHWAY_OUT,
    params:
        prefix=lambda wc, output: output[0].removesuffix(".pathway_summary.txt"),
        ranking_col=config["crispr"]["pathway"]["ranking_column"],
        permutations=config["crispr"]["pathway"]["permutations"],
    log:
        os.path.join(LOGS, "{cell_line}", "mageck_pathway_{sample}.log"),
    conda:
        "../envs/crispr.yaml"
    shell:
        """
        mageck pathway \
            --gene-ranking {input.gene_ranking} \
            --gmt-file {input.gmt} \
            --method gsea \
            --single-ranking \
            --ranking-column {params.ranking_col} \
            --output-prefix {params.prefix} \
            --permutation {params.permutations} > {log} 2>&1
        """


rule plot_SSMD:
    input:
        ssmd=replicate_pair_files(SSMD_SCORES),
        es=replicate_pair_files(PATHWAY_OUT),
        pos_con=config["crispr"]["positive_controls"],
    output:
        SSMD_FIGURE,
    params:
        neg_con=lambda wc: CELL_CFG[wc.cell_line]["neg_con"],
        script=script("plot_SSMD.py"),
    log:
        os.path.join(LOGS, "{cell_line}", "plot_SSMD.log"),
    conda:
        "../envs/crispr.yaml"
    shell:
        """
        python {params.script} \
            --ssmd_r1 {input.ssmd[0]} \
            --ssmd_r2 {input.ssmd[1]} \
            --pos_con {input.pos_con} \
            --ES_r1 {input.es[0]} \
            --ES_r2 {input.es[1]} \
            --output {output} \
            --cell_line {wildcards.cell_line} \
            --neg_con {params.neg_con} > {log} 2>&1
        """

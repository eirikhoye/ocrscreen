"""
Module crispri_abc_significant_hits.

Combines the two screen replicates of a cell line (MAGeCK pathway p-values, Fisher's
method + BH), keeps OCRs with combined FDR < 0.05 and p < 0.05 in both replicates,
and intersects them with ABC/STARE contacts (run config crispri_hits.abc_run).
Outputs the significant OCR–gene contacts (promoter/enhancer annotated) and a BED of
the significant OCRs (CREs), which is the input of the TF-enrichment modules.
Requires pathway, crispr_library and ABC_STARE outputs.
"""


rule crispri_abc_significant_hits:
    input:
        pathway=replicate_pair_files(PATHWAY_OUT),
        ssmd=replicate_pair_files(SSMD_SCORES),
        coords=LIBRARY_COORDS,
        abc=lambda wc: ABC_SCORES.format(cell_line=wc.cell_line, abc_run=config["crispri_hits"]["abc_run"]),
    output:
        contacts=SIG_CONTACTS,
        bed=SIG_CRE_BED,
    params:
        coord_col=lambda wc: f"Coord_{wc.cell_line}",
        script=script("crispri_ocr_abc_hits.py"),
    log:
        os.path.join(LOGS, "{cell_line}", "crispri_abc_significant_hits.log"),
    conda:
        "../envs/crispr.yaml"
    shell:
        """
        python {params.script} \
            {input.pathway[0]} {input.pathway[1]} \
            {input.ssmd[0]} {input.ssmd[1]} \
            {input.coords} {input.abc} {output.contacts} \
            --coord-col {params.coord_col} \
            --sig-bed {output.bed} > {log} 2>&1
        """

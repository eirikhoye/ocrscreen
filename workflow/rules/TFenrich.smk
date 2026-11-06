"""
Modules TFenrich and TFenrich_OCR_universe (shared rule file).

UniBind/LOLA enrichment of TF binding sites in the significant CRISPRi OCRs (BED from
module crispri_abc_significant_hits), after liftOver hg19 -> hg38:

  TFenrich               oneSetNoBg: background = UniBind universe
  TFenrich_OCR_universe  oneSetBg:   background = all DNase-seq OCRs of the cell line

Run with --use-singularity (unibind.container) or with a locally deployed
unibind_enrichment (set unibind.container: null).
"""

UNIBIND = config["unibind"]


rule liftover_sig_cre:
    input:
        bed=SIG_CRE_BED,
        chain=CHAIN_HG19_TO_HG38,
    output:
        bed=SIG_CRE_HG38,
        unmapped=SIG_CRE_UNMAPPED,
    log:
        os.path.join(LOGS, "{cell_line}", "liftover_sig_cre.log"),
    conda:
        "../envs/ucsc.yaml"
    shell:
        "liftOver {input.bed} {input.chain} {output.bed} {output.unmapped} > {log} 2>&1"


rule unibind_enrich:
    input:
        lola=UNIBIND["lola_db"],
        universe=UNIBIND["lola_universe"],
        bed=SIG_CRE_HG38,
    output:
        UNIBIND_OUT,
    params:
        outdir=lambda wc, output: os.path.dirname(output[0]),
        enrich=UNIBIND["enrich_script"],
    log:
        os.path.join(LOGS, "{cell_line}", "unibind_enrich.log"),
    container:
        UNIBIND["container"]
    shell:
        """
        bash {params.enrich} oneSetNoBg \
            {input.lola} {input.universe} {input.bed} {params.outdir} > {log} 2>&1
        """


rule create_ocr_universe_bed:
    input:
        narrowpeak=lambda wc: CELL_CFG[wc.cell_line]["dnase_peaks"],
    output:
        bed=OCR_UNIVERSE_HG19,
    log:
        os.path.join(LOGS, "{cell_line}", "create_ocr_universe_bed.log"),
    conda:
        "../envs/download.yaml"
    shell:
        "cut -f1-3 {input.narrowpeak} > {output.bed} 2> {log}"


rule liftover_ocr_universe:
    input:
        bed=OCR_UNIVERSE_HG19,
        chain=CHAIN_HG19_TO_HG38,
    output:
        bed=OCR_UNIVERSE_HG38,
        unmapped=OCR_UNIVERSE_UNMAPPED,
    log:
        os.path.join(LOGS, "{cell_line}", "liftover_ocr_universe.log"),
    conda:
        "../envs/ucsc.yaml"
    shell:
        "liftOver {input.bed} {input.chain} {output.bed} {output.unmapped} > {log} 2>&1"


rule unibind_enrich_ocr_universe:
    input:
        lola=UNIBIND["lola_db"],
        bed=SIG_CRE_HG38,
        universe=OCR_UNIVERSE_HG38,
    output:
        UNIBIND_OUT_OCR,
    params:
        outdir=lambda wc, output: os.path.dirname(output[0]),
        enrich=UNIBIND["enrich_script"],
    log:
        os.path.join(LOGS, "{cell_line}", "unibind_enrich_ocr_universe.log"),
    container:
        UNIBIND["container"]
    shell:
        """
        bash {params.enrich} oneSetBg \
            {input.lola} {input.bed} {input.universe} {params.outdir} > {log} 2>&1
        """

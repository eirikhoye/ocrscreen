"""
Module ABC_STARE.

Activity-by-contact (ABC) enhancer–gene scores and TF–gene affinities with STARE,
using the cell line's own KR-normalised Hi-C contacts (module hic2KR, or a precomputed
contact directory via set_HIC). Each entry of config abc_stare.runs combines an
enhancer-activity source with STARE options:

  dnase      ENCODE/Duke DNase-seq peaks of the cell line (narrowPeak signal, col 7)
  jacks_fdr  CRISPRi enhancer activity of JACKS-significant OCRs (module JACKS)
  jacks_all  CRISPRi enhancer activity of all OCRs (module JACKS)

A stable copy of the scored interactions is written to
results/<CELL>/abc_stare/<run>/<CELL>_<run>_ABC_scoredInteractions.txt.gz.
"""

rule run_stare:
    input:
        bed=abc_activity_bed,
        contacts=lambda wc: contacts_dir(wc.cell_line),
        gtf=abc_gtf,
        fasta=GENOME_FASTA,
        psem=PSEM,
    output:
        metadata=ABC_METADATA,
        scores=ABC_SCORES,
    params:
        outdir=ABC_DIR,
        column=lambda wc: ABC_RUNS[wc.abc_run]["column"],
        tss=lambda wc: ABC_RUNS[wc.abc_run]["tss"],
        bin_size=KR["bin_size"],
        rmtree=script("stare_rmtree.py"),
    threads: config["abc_stare"]["threads"]
    log:
        os.path.join(LOGS, "{cell_line}", "stare_{abc_run}.log"),
    conda:
        "../envs/stare.yaml"
    shell:
        """
        # STARE refuses to write into a pre-existing output folder
        python {params.rmtree} --metadata {output.metadata} --outdir {params.outdir}

        STARE.sh -c {threads} \
            -b {input.bed} \
            -f {input.contacts} \
            -a {input.gtf} \
            -g {input.fasta} \
            -s {input.psem} \
            -o {params.outdir} \
            -n {params.column} \
            -k {params.bin_size} \
            -i {params.tss} > {log} 2>&1

        cp {params.outdir}/ABC_output/{wildcards.abc_run}_ABCpp_scoredInteractions_*.txt.gz {output.scores}
        """


rule abc_to_bedpe:
    input:
        abc_scores=ABC_SCORES,
    output:
        bedpe=ABC_BEDPE,
    log:
        os.path.join(LOGS, "{cell_line}", "abc_to_bedpe_{abc_run}.log"),
    conda:
        "../envs/crispr.yaml"
    script:
        "../scripts/abc_to_bedpe.py"

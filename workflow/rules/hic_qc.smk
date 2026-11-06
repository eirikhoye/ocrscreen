"""
Module hic_qc (optional).

Reproducibility between Hi-C libraries of a cell line: stratum-adjusted correlation
coefficient (SCC; HiCRep, Yang et al. 2017) for every library pair, computed on the
per-library .cool files at hic.qc_resolution. Requires pairtools outputs.
"""


rule hicrep_scc:
    input:
        lambda wc: expand(LIB_COOL, cell_line=wc.cell_line, replicate=hic_replicates_of(wc.cell_line)),
    output:
        HICREP,
    params:
        h=config["hic"]["hicrep"]["h"],
        dbpmax=config["hic"]["hicrep"]["dbpmax"],
    log:
        os.path.join(LOGS, "{cell_line}", "hic", "hicrep.log"),
    conda:
        "../envs/hicrep.yaml"
    script:
        "../scripts/hicrep_pairwise.py"

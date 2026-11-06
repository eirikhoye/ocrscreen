"""Keep ABC interactions whose enhancer (OCR) has JACKS FDR below a threshold.

Snakemake script (rule ``filter_ABC_for_JACKS_fdr``); logic from the original
crispr_pipeline RUN_GSEA.smk. Note: as in the original, the ``fdr`` input is the
JACKS gene p-value table and its two columns are renamed to OCR/FDR.
"""
import pandas as pd

abc_df = pd.read_csv(snakemake.input.abc, sep="\t")
fdr_df = pd.read_csv(snakemake.input.fdr, sep="\t")
fdr_df.columns = ["OCR", "FDR"]
bed_df = pd.read_csv(snakemake.input.bed, sep="\t")

# PeakID as "CHR:START-END" (no 'chr' prefix), as written by STARE
bed_df["PeakID"] = (
    bed_df["#chr"].astype(str).replace("^chr", "", regex=True)
    + ":" + bed_df["start"].astype(str)
    + "-" + bed_df["end"].astype(str)
)

merged_df = pd.merge(bed_df, fdr_df, on="OCR")
abc_merged = pd.merge(abc_df, merged_df, on="PeakID")
abc_merged = abc_merged[abc_merged["FDR"] < snakemake.params.fdr]

abc_merged.to_csv(snakemake.output[0], sep="\t", index=False)

"""Convert STARE ABC scored interactions to BEDPE (enhancer -> TSS links).

Snakemake script (rule ``abc_to_bedpe``); logic from the original ``create_bedpe`` rule.
"""
import pandas as pd

abc_df = pd.read_csv(snakemake.input.abc_scores, sep="\t")
links_df = abc_df[["#chr", "start", "end", "Ensembl ID", "TSS-dist", "ABC-Score"]].copy()

links_df["gene_start"] = (links_df["start"] + links_df["TSS-dist"]).round().astype(int)
links_df["gene_end"] = links_df["gene_start"] + 1  # single bp for gene anchor

bedpe_df = pd.DataFrame({
    "chrom1": links_df["#chr"],
    "start1": links_df["start"],
    "end1": links_df["end"],
    "chrom2": links_df["#chr"],
    "start2": links_df["gene_start"],
    "end2": links_df["gene_end"],
    "name": links_df["Ensembl ID"],
    "score": links_df["ABC-Score"],
})

bedpe_df.to_csv(snakemake.output.bedpe, sep="\t", header=False, index=False)

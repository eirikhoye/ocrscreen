"""Rank genes by their summed ABC score (input .rnk for GSEApy prerank).

Snakemake script (rule ``rank_genes_for_GSEA``); logic from the original RUN_GSEA.smk.
"""
import pandas as pd

abc_df = pd.read_csv(snakemake.input[0], sep="\t")

gene_sums = abc_df.groupby("Gene Name", as_index=False)["ABC-Score"].sum()
gene_sums = gene_sums.rename(columns={"ABC-Score": "ABC-Score_Sum"})
gene_sums_sorted = gene_sums.sort_values(by="ABC-Score_Sum", ascending=False)

ranks = gene_sums_sorted[["Gene Name", "ABC-Score_Sum"]]
ranks.columns = ["gene_id", "score"]

# .rnk: tab-delimited, no header, no index
ranks.to_csv(snakemake.output[0], sep="\t", index=False, header=False)

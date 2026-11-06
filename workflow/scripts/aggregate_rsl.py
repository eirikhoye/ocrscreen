"""Sum binned RSL counts per sgRNA across all screen replicates of a cell line.

Snakemake script (rule ``aggregate_RSL``).
"""
from functools import reduce

import pandas as pd

dfs = [pd.read_csv(f, sep=",") for f in snakemake.input]

# Merge on 'RSL.guide' and 'guide.set'
merged = reduce(
    lambda left, right: pd.merge(left, right, on=["RSL.guide", "guide.set"], how="outer"), dfs
)

# Aggregate bins into total sgRNA counts
value_cols = [col for col in merged.columns if col not in ["RSL.guide", "guide.set"]]
aggregated = merged.groupby("guide.set")[value_cols].sum().reset_index()

aggregated.to_csv(snakemake.output[0], index=False)

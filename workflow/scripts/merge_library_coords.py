"""Add genomic coordinates to the sgRNA library of one cell line.

Snakemake script (rule ``merge_library_coords``). Logic taken unchanged from the
``merge_coords`` rule of the original crispr_pipeline ``add_coords`` module.

Outputs
-------
coords  : library with an extra ``Coord_<CELL>`` column (CSV, header).
library : the library with a header row (sgRNA,sequence,open_reg). This headed
          copy is what CrisprCounter and JACKS use downstream (CrisprCounter skips
          the first line of the design file).
"""
import pandas as pd

columns = snakemake.params.coord_cols  # e.g. ["N20", "Coord_RWPE"]

# Library has no header: sgRNA id, protospacer sequence, targeted open region
ref_df = pd.read_csv(snakemake.input.ref, header=None, names=["sgRNA", "sequence", "open_reg"])

coord_df = pd.read_csv(snakemake.input.coord, sep="\t")
coord_df = coord_df.drop_duplicates(columns[0])
coord_df = coord_df[columns].rename(columns={columns[0]: "sequence"})

merged_df = ref_df.merge(coord_df, on="sequence", how="left")
merged_df[columns[1]] = merged_df[columns[1]].fillna("NA")

merged_df.to_csv(snakemake.output.coords, index=False)
ref_df.to_csv(snakemake.output.library, index=False)

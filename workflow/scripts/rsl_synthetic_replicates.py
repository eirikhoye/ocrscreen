"""Use RSL bins as synthetic replicates for JACKS and write the JACKS replicate map.

Snakemake script (rule ``convert_RSL_to_synthetic_replicates``). Logic taken from
the original crispr_pipeline (JACKS.smk + common.smk::pivot_counts_by_umi).
"""
import os

import pandas as pd


def pivot_counts_by_umi(file_paths, umi_length=2):
    """Extract the RSL bin label of each count row and pivot to guide x (condition_bin)."""
    all_long = []
    for file_path in file_paths:
        df = pd.read_csv(file_path)
        sample_name = os.path.splitext(os.path.basename(file_path))[0]

        df["UMI"] = df["RSL.guide"].str.extract(fr"_([A-Za-z0-9]{{{umi_length}}})$")
        cols = ["guide.set", "UMI"] + [col for col in df.columns if col.startswith(sample_name)]
        df = df[cols]

        df_long = df.melt(id_vars=["guide.set", "UMI"], var_name="Condition", value_name="Count")
        df_long["Condition_UMI"] = df_long["Condition"] + "_" + df_long["UMI"]
        all_long.append(df_long)

    combined_long = pd.concat(all_long, ignore_index=True)
    df_wide = combined_long.pivot(index="guide.set", columns="Condition_UMI", values="Count")
    df_wide = df_wide.fillna(0).astype(int)
    return df_wide.reset_index()


# Library (with header sgRNA,sequence,open_reg)
ref = pd.read_csv(snakemake.input.ref, sep=",")

df = pivot_counts_by_umi(list(snakemake.input.binned), snakemake.params.nr_nt)

merged = df.merge(ref[["sgRNA", "open_reg"]], left_on="guide.set", right_on="sgRNA", how="left")
first_cols = ["guide.set", "open_reg"]
rsl_replicates = [col for col in merged.columns if col not in first_cols + ["sgRNA"]]

merged = merged[first_cols + rsl_replicates].reset_index(drop=True)
merged.to_csv(snakemake.output.rep, index=False)

# Replicate map: <label>_<bin>  ->  <cell_line>_<stage>, control
rows = snakemake.params.rows
with open(snakemake.output.repmap, "w") as out:
    out.write("Replicate\tSample\tControl\n")
    for rsl in rsl_replicates:
        base_name = "_".join(rsl.split("_")[:-1])  # strip the RSL bin suffix
        for info in rows:
            if base_name in info["label"]:
                out.write(f"{rsl}\t{info['cell_line']}_{info['stage']}\t{info['control']}\n")
                break
        else:
            raise ValueError(f"Base replicate '{base_name}' not found in the CRISPR sample sheet")

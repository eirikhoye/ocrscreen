"""Write the per-replicate sample table read by makeCounterConfig.py.

Snakemake script (rule ``create_sample_table``). Columns: file, sample, stage, group.
"""
import csv

with open(snakemake.output[0], "w", newline="") as csvfile:
    writer = csv.writer(csvfile, delimiter="\t")
    writer.writerow(["file", "sample", "stage", "group"])
    for row in snakemake.params.rows:
        writer.writerow([row["fastq"], row["label"], row["stage"], row["group"]])

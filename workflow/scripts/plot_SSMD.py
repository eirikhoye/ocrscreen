import argparse
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import csv
from scipy.stats import pearsonr
from sklearn.linear_model import LinearRegression
import numpy as np


def read_csv_auto_sep(filepath, header=None):
    with open(filepath, 'r') as f:
        sample = f.read(1024)
        dialect = csv.Sniffer().sniff(sample)
        f.seek(0)
        return pd.read_csv(f, delimiter=dialect.delimiter, header=header)


def classify_gene(value, pos_con, cont_label):
    if cont_label in value:
        return "Negative Control"
    if any(gene in value for gene in pos_con):
        return "Positive"
    else:
        return "Open Chromatin"


# Argument parsing
parser = argparse.ArgumentParser(description="Plot SSMD and Score replicate correlation")
parser.add_argument("--ssmd_r1", required=True, help="Path to <CELL> R1 SSMD file")
parser.add_argument("--ssmd_r2", required=True, help="Path to <CELL> R2 SSMD file")
parser.add_argument("--pos_con", required=True, help="Path to list of positive control genes")
parser.add_argument("--ES_r1", required=True, help="Path to <CELL> R1 pathway summary file")
parser.add_argument("--ES_r2", required=True, help="Path to <CELL> R2 pathway summary file")
parser.add_argument("--output", required=True, help="Output plot file (e.g. plot.pdf)")
parser.add_argument("--cell_line", required=True, help="Cell line to show as file name")
parser.add_argument("--neg_con", default="NCR", help="Prefix of negative-control guide sets (e.g. NCR, NCH)")

args = parser.parse_args()

CELL_LINE = args.cell_line

# Load positive control genes
pos_con_data = read_csv_auto_sep(args.pos_con)
pos_con = [i.split('_')[0] for i in pos_con_data[0]]

# Load SSMD data
SSMD_R1 = pd.read_csv(args.ssmd_r1, sep='\t')
SSMD_R2 = pd.read_csv(args.ssmd_r2, sep='\t')

# Merge and classify
SSMD_merged = pd.merge(SSMD_R1, SSMD_R2, on="guide.set", how="outer")
SSMD_merged['type'] = SSMD_merged['guide.set'].apply(lambda row: classify_gene(row, pos_con, args.neg_con))

# Load enrichment score
ES_R1 = pd.read_csv(args.ES_r1, sep='\t')
ES_R2 = pd.read_csv(args.ES_r2, sep='\t')

ES_merged = pd.merge(ES_R1, ES_R2, on='Pathway', suffixes=('_R1', '_R2'))
ES_merged["FDR < 0.1 both"] = (ES_merged["FDR_R1"] < 0.1) & (ES_merged["FDR_R2"] < 0.1)

# Set plot style
sns.set_theme(style="whitegrid", context="talk", font_scale=1.2)
hue_order = ["Open Chromatin", "Positive", "Negative Control"]

# Create figure
fig, axes = plt.subplots(2, 2, figsize=(12, 14))

# SSMD plot
x1 = SSMD_merged[f'{CELL_LINE}_R1_D14.SSMD']
y1 = SSMD_merged[f'{CELL_LINE}_R2_D14.SSMD']
combined = pd.concat([x1, y1], axis=1).dropna()
r_ssmd, p_ssmd = pearsonr(combined.iloc[:, 0], combined.iloc[:, 1])

x_vals = combined.iloc[:, 0].values.reshape(-1, 1)
y_vals = combined.iloc[:, 1].values

model = LinearRegression().fit(x_vals, y_vals)

x_range = np.linspace(combined.iloc[:, 0].min(), combined.iloc[:, 0].max(), 100)
y_pred = model.predict(x_range.reshape(-1, 1))

axes[0, 0].plot(x_range, y_pred, color='black', linestyle='--')

sns.scatterplot(
    ax=axes[0, 0],
    data=SSMD_merged,
    x=f'{CELL_LINE}_R1_D14.SSMD',
    y=f'{CELL_LINE}_R2_D14.SSMD',
    hue='type',
    hue_order=hue_order,
    palette='Set1',
    s=40,
    edgecolor='black',
    linewidth=0.5,
    alpha=0.2,
    legend=False
)

axes[0, 0].set_xlabel(f"SSMD: {CELL_LINE} R1 D14")
axes[0, 0].set_ylabel(f"SSMD: {CELL_LINE} R2 D14")
axes[0, 0].set_title(f"SSMD Replicate Correlation\nr = {r_ssmd:.2f}, p={p_ssmd:.1e}")

# Score plot
x2 = SSMD_merged[f'{CELL_LINE}_R1_D14.Score']
y2 = SSMD_merged[f'{CELL_LINE}_R2_D14.Score']
combined = pd.concat([x2, y2], axis=1).dropna()
r_score, p_score = pearsonr(combined.iloc[:, 0], combined.iloc[:, 1])

x_vals = combined.iloc[:, 0].values.reshape(-1, 1)
y_vals = combined.iloc[:, 1].values
model = LinearRegression().fit(x_vals, y_vals)
x_range = np.linspace(combined.iloc[:, 0].min(), combined.iloc[:, 0].max(), 100)
y_pred = model.predict(x_range.reshape(-1, 1))

axes[0, 1].plot(x_range, y_pred, color='black', linestyle='--')


scatter2 = sns.scatterplot(
    ax=axes[0, 1],
    data=SSMD_merged,
    x=f'{CELL_LINE}_R1_D14.Score',
    y=f'{CELL_LINE}_R2_D14.Score',
    hue='type',
    hue_order=hue_order,
    palette='Set1',
    s=40,
    edgecolor='black',
    linewidth=0.5,
    alpha=0.2,
    legend=True
)
axes[0, 1].set_xlabel(f"Score: {CELL_LINE} R1 D14")
axes[0, 1].set_ylabel(f"Score: {CELL_LINE} R2 D14")
axes[0, 1].set_title(f"Score Replicate Correlation\nr = {r_score:.2f}, p={p_score:.1e}")

# GSEA plot
ax_gsea = plt.subplot2grid((2, 2), (1, 0), colspan=2)

# Plot FDR ≥ 0.1
sns.scatterplot(
    ax=ax_gsea,
    x='ES_R1',
    y='ES_R2',
    data=ES_merged[ES_merged["FDR < 0.1 both"] == False],
    color='grey',
    alpha=0.1,
    legend=False
)

# Hide all visual elements from this plot except the dots
ax_gsea.set_xticks([])
ax_gsea.set_yticks([])
ax_gsea.set_xlabel("")
ax_gsea.set_ylabel("")
ax_gsea.set_title("")
ax_gsea.spines['top'].set_visible(False)
ax_gsea.spines['right'].set_visible(False)
ax_gsea.spines['bottom'].set_visible(False)
ax_gsea.spines['left'].set_visible(False)

# Plot FDR < 0.1
sns.scatterplot(
    ax=ax_gsea,
    x='ES_R1',
    y='ES_R2',
    data=ES_merged[ES_merged["FDR < 0.1 both"] == True],
    color='red',
    alpha=0.8,
    legend=False
)

ax_gsea.set_xlabel(f"Enrichment Score: {CELL_LINE} R1")
ax_gsea.set_ylabel(f"Enrichment Score: {CELL_LINE} R2")
ax_gsea.set_title("GSEA Enrichment Scores R1 vs R2")

# Remove default legend and add clean shared one
#axes[1].get_legend().remove()
#handles, labels = axes[1].get_legend_handles_labels()
#fig.legend(
#    handles,
#    labels,
#    loc='lower center',
#    bbox_to_anchor=(0.5, -0.05),
#    ncol=len(labels),
#    frameon=False
#)

num_sig = ES_merged["FDR < 0.1 both"].sum()

ax_gsea.text(
    0.95, 0.05,
    f"FDR < 0.1 in both: {num_sig}",
    transform=ax_gsea.transAxes,
    ha='right',
    va='bottom',
    fontsize=12,
    color='red'
    )


# Finalize and save
sns.despine(trim=True)
plt.tight_layout(rect=[0, 0.08, 1, 1])
plt.savefig(args.output, dpi=300)
